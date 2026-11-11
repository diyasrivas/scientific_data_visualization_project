"""EAS vs SAS divergence at CPIC variants; compare CPIC level A genes with the rest."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap

sys.path.insert(0, str(Path(__file__).parent))
import style  # noqa: E402

RAW, OUT, FIG = Path("data/raw"), Path("data/processed"), Path("figures")
BG_FILE = OUT / "background_gene_fst.csv"
N_TOP, N_LABEL = 20, 6
LABEL_OFFSETS = [(7, 4, "left", "bottom"), (7, -4, "left", "top"), (-7, 4, "right", "bottom"),
                 (-7, -4, "right", "top"), (0, 14, "center", "bottom"), (0, -14, "center", "top")]
SEED, MIN_AN, MIN_MAF, MIN_VAR = 0, 1000, 0.01, 5
rng = np.random.default_rng(SEED)


def hudson_terms(p1, n1, p2, n2):
    """Per-variant numerator and denominator of Hudson's FST estimator."""
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def load():
    d = pd.read_csv(RAW / "allpop_pharmacogenes.csv")[["gene", "variant_id", "ac_eas", "an_eas", "ac_sas", "an_sas"]]
    d.columns = ["gene", "variant_id", "eas_ac", "eas_an", "sas_ac", "sas_an"]
    lv = pd.read_csv(RAW / "genes.csv").rename(columns={"genesymbol": "gene"})
    d = d.merge(lv, on="gene")
    n_all = len(d)
    d = d[(d.eas_an >= MIN_AN) & (d.sas_an >= MIN_AN)].copy()
    d["p_eas"], d["p_sas"] = d.eas_ac / d.eas_an, d.sas_ac / d.sas_an
    maf = np.maximum(np.minimum(d.p_eas, 1 - d.p_eas), np.minimum(d.p_sas, 1 - d.p_sas))
    d = d[maf >= MIN_MAF].copy()  # common in at least one group
    d["num"], d["den"] = hudson_terms(d.p_eas, d.eas_an, d.p_sas, d.sas_an)
    d["delta"] = (d.p_eas - d.p_sas).abs()
    d["level_A"] = d.level.eq("A")
    d["genesymbol"] = d.gene
    print(f"{n_all} variants (AF>=0.1%) -> {len(d)} common (MAF>=1%) and well covered, "
          f"in {d.gene.nunique()} genes")
    return d


def gene_fst(d):
    g = d.groupby("genesymbol").agg(num=("num", "sum"), den=("den", "sum"),
                                    n_var=("num", "size"), level_A=("level_A", "max"),
                                    mean_delta=("delta", "mean"))
    g["fst"] = g.num / g.den
    g = g[g.n_var >= MIN_VAR]
    g = g.join(drop_one_min(d))
    return g.sort_values("fst", ascending=False)


def drop_one_min(d):
    """Lowest gene FST obtainable by removing one variant: how much of a gene's value is one variant.

    A delete-one jackknife asked the other way round. A gene whose FST survives this is divergent
    across its variants; a gene whose FST collapses owes its rank to a single site.
    """
    out = {}
    for gene, s in d.groupby("genesymbol"):
        num, den = s.num.to_numpy(), s.den.to_numpy()
        tot_n, tot_d = num.sum(), den.sum()
        out[gene] = np.min((tot_n - num) / (tot_d - den))
    return pd.Series(out, name="fst_drop1").rename_axis("genesymbol")


def perm_test(g, n=10000):
    """Is mean gene-level FST larger for level A genes? Permute the gene labels."""
    a, x = g.level_A.to_numpy(), g.fst.to_numpy()
    obs = x[a].mean() - x[~a].mean()
    null = np.empty(n)
    for i in range(n):
        sh = rng.permutation(a)
        null[i] = x[sh].mean() - x[~sh].mean()
    return obs, (np.sum(null >= obs) + 1) / (n + 1)


def save(fig, name):
    fig.savefig(FIG / f"{name}.svg"); fig.savefig(FIG / f"{name}.png", dpi=220); plt.close(fig)


def load_background():
    """FST of random non-pharmacogenes, if that analysis has been run."""
    return pd.read_csv(BG_FILE).fst.to_numpy() if BG_FILE.exists() else None


def ecdf(x):
    x = np.sort(x)
    return x, np.arange(1, len(x) + 1) / len(x)


def title_block(fig, title, subtitle=None, x=0.02):
    """Neutral title (and optional subtitle) aligned to the figure's left edge."""
    fig.text(x, 0.975, title, fontsize=15, fontweight="bold", va="top")
    if subtitle:
        fig.text(x, 0.915, subtitle, fontsize=style.KEY_FONT, color=style.SUBTLE, va="top")


def margin_key(fig, x, y, items, dy=0.055):
    """Coloured-text key in the right margin: (text, colour) rows."""
    for i, (txt, c) in enumerate(items):
        fig.text(x, y - i * dy, "\u25cf " + txt, color=c, fontsize=style.KEY_FONT, fontweight="bold", va="top")


def fig_ranking(g, bg):
    top = g.head(N_TOP).iloc[::-1]
    col = [style.LEVEL_A if b else style.OTHER_CPIC for b in top.level_A]
    fig, ax = plt.subplots(figsize=(10, 6.6))
    fig.subplots_adjust(left=0.19, right=0.76, top=0.87, bottom=0.12)
    if bg is not None:
        q1, med, q3 = np.percentile(bg, [25, 50, 75])
        ax.axvspan(q1, q3, color=style.BACKGROUND, alpha=.5, lw=0, zorder=0)
        ax.axvline(med, color=style.SUBTLE, lw=1.5, zorder=1)
        band, lo = f"Random non-pharmacogenes (n={len(bg)}): median and IQR", q1
    else:
        ax.axvline(g.fst.median(), color=style.SUBTLE, lw=1.5, zorder=0)
        band, lo = f"Median of all {len(g)} CPIC genes", g.fst.median()
    ax.text(lo, 1.015, band, transform=ax.get_xaxis_transform(), ha="left", va="bottom",
            fontsize=style.KEY_FONT, color=style.SUBTLE)  # label sits on the band, above its top edge
    ax.hlines(range(len(top)), 0, top.fst, color=col, lw=3)
    # hollow marker: the same gene with its single most divergent variant removed
    ax.scatter(top.fst_drop1, range(len(top)), facecolor="white", edgecolor=col, linewidth=1.8, s=55, zorder=4)
    ax.scatter(top.fst, range(len(top)), c=col, s=80, zorder=3)
    rug_y = -1.6  # context: where all genes sit, on the same x scale
    ax.vlines(g.fst, rug_y - .35, rug_y + .35, color=style.BACKGROUND, lw=1.2)
    ax.set_yticks(list(range(len(top))) + [rug_y])
    ax.set_yticklabels(list(top.index) + [f"All {len(g)} genes"])
    for lab, flag in zip(ax.get_yticklabels(), list(top.level_A) + [None]):
        lab.set_color(style.LEVEL_A if flag else style.TEXT if flag is False else style.SUBTLE)
    ax.set_ylim(-2.3, N_TOP - .4)
    ax.set_xlim(0, top.fst.max() * 1.05)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Hudson FST, East Asian vs South Asian (dimensionless)")
    ax.set_ylabel(f"Gene (top {N_TOP} of {len(g)} by FST)")
    title_block(fig, "Gene-level East vs South Asian FST, top pharmacogenes")
    margin_key(fig, 0.78, 0.86, [("CPIC level A", style.LEVEL_A), ("other CPIC genes", style.OTHER_CPIC)])
    fig.text(0.78, 0.68, "\u25cb  without the gene's\n     single most\n     divergent variant", color=style.TEXT,
             fontsize=style.KEY_FONT, va="top", linespacing=1.25)
    save(fig, "f1_gene_ranking")


def fig_scatter(d, g):
    fig, ax = plt.subplots(figsize=(10, 6.4))
    fig.subplots_adjust(left=0.09, right=0.62, top=0.86, bottom=0.14)
    cmap = LinearSegmentedColormap.from_list("grey", ["#F0F0F0", "#BDBDBD"])
    ax.hexbin(d.p_sas, d.p_eas, gridsize=40, extent=(0, 1, 0, 1), mincnt=1, cmap=cmap, bins="log", linewidths=0)
    ax.plot([0, 1], [0, 1], color=style.SUBTLE, lw=1.8, ls="--", zorder=2)
    ax.text(.74, .70, "y = x", color=style.SUBTLE, fontsize=style.KEY_FONT, rotation=45, ha="center", va="top")
    fig.canvas.draw()
    taken = []
    rend = fig.canvas.get_renderer()
    labelled = [d[d.genesymbol == gene].sort_values("delta").iloc[-1] for gene in g.head(N_LABEL).index]
    for r in labelled:
        c = style.LEVEL_A if r.level_A else style.OTHER_CPIC
        ax.scatter(r.p_sas, r.p_eas, s=100, color=c, edgecolor="white", linewidth=1.2, zorder=3)
        x, y = ax.transData.transform((r.p_sas, r.p_eas))
        taken.append((x - 9, y - 9, x + 9, y + 9))
    for r in labelled:  # first label offset whose box is free of other labels and points
        for dx, dy, ha, va in LABEL_OFFSETS:
            t = ax.annotate(r.genesymbol, (r.p_sas, r.p_eas), fontsize=style.KEY_FONT, fontweight="bold", zorder=4,
                            color=style.LEVEL_A if r.level_A else style.TEXT, xytext=(dx, dy),
                            textcoords="offset points", ha=ha, va=va)
            bb = t.get_window_extent(rend)
            if all(bb.x1 < q[0] or bb.x0 > q[2] or bb.y1 < q[1] or bb.y0 > q[3] for q in taken):
                taken.append((bb.x0, bb.y0, bb.x1, bb.y1))
                break
            t.remove()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("Alt allele frequency, South Asian (fraction)"); ax.set_ylabel("Alt allele frequency, East Asian (fraction)")
    title_block(fig, "Alt allele frequency per variant, East vs South Asian", f"n = {len(d)} common variants", x=0.02)
    fig.text(0.65, 0.80, "Grey hexagons:\nvariant density\n(darker = more)", color=style.SUBTLE,
             fontsize=style.KEY_FONT, va="top")
    fig.text(0.65, 0.60, f"Points: most divergent\nvariant of each of the\ntop {N_LABEL} genes by FST", color=style.TEXT,
             fontsize=style.KEY_FONT, va="top")
    margin_key(fig, 0.65, 0.38, [("CPIC level A", style.LEVEL_A), ("other CPIC genes", style.OTHER_CPIC)])
    save(fig, "f2_allele_freq_scatter")


def fig_ecdf(g, p, bg):
    """ECDF per group; each label sits at the top end of its own curve."""
    groups = [(g[~g.level_A].fst.to_numpy(), style.OTHER_CPIC, "-", "other CPIC genes", "left", "bottom", 1.0),
              (g[g.level_A].fst.to_numpy(), style.LEVEL_A, "-", "CPIC level A genes", "right", "top", .93)]
    if bg is not None:
        groups.insert(0, (bg, style.SUBTLE, (0, (4, 2)), "random non-pharmacogenes", "right", "bottom", 1.0))
    fig, ax = plt.subplots(figsize=(10, 5.6))
    fig.subplots_adjust(left=0.09, right=0.78, top=0.90, bottom=0.15)
    for x, c, ls, lab, ha, va, ty in groups:
        xs, ys = ecdf(x)
        ax.step(xs, ys, where="post", color=c, ls=ls, lw=3 if ls == "-" else 2.4)
        tx = xs[-1] + (.005 if ha == "left" else -.005)
        ax.text(tx, ty + (.02 if va == "bottom" else 0), f"{lab}\n(n={len(x)})", color=c, fontsize=style.KEY_FONT,
                fontweight="bold", ha=ha, va=va)
    ax.set_ylim(0, 1.2); ax.set_xlim(0, 0.43)
    ax.set_yticks([0, .2, .4, .6, .8, 1.0])
    ax.set_xlabel("Gene-level Hudson FST, East vs South Asian (dimensionless)")
    ax.set_ylabel("Cumulative fraction of genes")
    title_block(fig, "Distribution of gene-level FST by CPIC level")
    fig.text(0.80, 0.86, "Level A vs other\nCPIC genes (one-sided\npermutation test on\nmean gene FST):\n"
             f"p = {p:.3f}", color=style.TEXT, fontsize=style.KEY_FONT, va="top")
    save(fig, "f3_levelA_vs_rest")


# ---- deck variants (SDV_DECK=1): slide-sized, no embedded title or footnote, larger text
def fig_ranking_deck(g, bg):
    n = 15
    top = g.head(n).iloc[::-1]
    col = [style.LEVEL_A if b else style.OTHER_CPIC for b in top.level_A]
    fig, ax = plt.subplots(figsize=(8.9, 4.55))
    fig.subplots_adjust(left=0.175, right=0.97, top=0.895, bottom=0.17)
    q1, med, q3 = np.percentile(bg, [25, 50, 75])
    ax.axvspan(q1, q3, color=style.BACKGROUND, alpha=.5, lw=0, zorder=0)
    ax.axvline(med, color=style.SUBTLE, lw=1.8, zorder=1)
    ax.text(q1, 1.015, f"Random non-pharmacogenes (n={len(bg)}): median and IQR", transform=ax.get_xaxis_transform(),
            ha="left", va="bottom", fontsize=13, color=style.SUBTLE)
    ax.hlines(range(len(top)), 0, top.fst, color=col, lw=3.5)
    ax.scatter(top.fst_drop1, range(len(top)), facecolor="white", edgecolor=col, linewidth=2, s=70, zorder=4)
    ax.scatter(top.fst, range(len(top)), c=col, s=95, zorder=3)
    rug_y = -1.6
    ax.vlines(g.fst, rug_y - .35, rug_y + .35, color=style.BACKGROUND, lw=1.5)
    ax.set_yticks(list(range(len(top))) + [rug_y])
    ax.set_yticklabels(list(top.index) + [f"All {len(g)} genes"], fontsize=12)
    for lab, flag in zip(ax.get_yticklabels(), list(top.level_A) + [None]):
        lab.set_color(style.LEVEL_A if flag else style.TEXT if flag is False else style.SUBTLE)
    ax.set_ylim(-2.3, n - .4)
    ax.set_xlim(0, top.fst.max() * 1.05)
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.set_xlabel("Hudson FST, East Asian vs South Asian (dimensionless)")
    ax.set_ylabel(f"Gene (top {n} of {len(g)} by FST)")
    margin_key(fig, 0.62, 0.62, [("CPIC level A", style.LEVEL_A), ("other CPIC genes", style.OTHER_CPIC)], dy=0.085)
    fig.text(0.62, 0.42, "\u25cb  without the gene's single\n     most divergent variant", color=style.TEXT,
             fontsize=13, va="top", linespacing=1.25)
    style.deck_save(fig, "f1_gene_ranking")


def fig_scatter_deck(d, g):
    fig, ax = plt.subplots(figsize=(8.9, 4.55))
    fig.subplots_adjust(left=0.095, right=0.505, top=0.955, bottom=0.165)
    cmap = LinearSegmentedColormap.from_list("grey", ["#F0F0F0", "#BDBDBD"])
    ax.hexbin(d.p_sas, d.p_eas, gridsize=40, extent=(0, 1, 0, 1), mincnt=1, cmap=cmap, bins="log", linewidths=0)
    ax.plot([0, 1], [0, 1], color=style.SUBTLE, lw=2, ls="--", zorder=2)
    ax.text(.74, .70, "y = x", color=style.SUBTLE, fontsize=13, rotation=45, ha="center", va="top")
    fig.canvas.draw()
    taken = []
    rend = fig.canvas.get_renderer()
    labelled = [d[d.genesymbol == gene].sort_values("delta").iloc[-1] for gene in g.head(N_LABEL).index]
    for r in labelled:
        c = style.LEVEL_A if r.level_A else style.OTHER_CPIC
        ax.scatter(r.p_sas, r.p_eas, s=120, color=c, edgecolor="white", linewidth=1.2, zorder=3)
        x, y = ax.transData.transform((r.p_sas, r.p_eas))
        taken.append((x - 12, y - 12, x + 12, y + 12))
    for r in labelled:
        for dx, dy, ha, va in LABEL_OFFSETS:
            t = ax.annotate(r.genesymbol, (r.p_sas, r.p_eas), fontsize=13, fontweight="bold", zorder=4,
                            color=style.LEVEL_A if r.level_A else style.TEXT, xytext=(dx, dy),
                            textcoords="offset points", ha=ha, va=va)
            bb = t.get_window_extent(rend)
            if all(bb.x1 < q[0] or bb.x0 > q[2] or bb.y1 < q[1] or bb.y0 > q[3] for q in taken):
                taken.append((bb.x0, bb.y0, bb.x1, bb.y1))
                break
            t.remove()
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_aspect("equal")
    ax.set_xlabel("Alt allele frequency, South Asian (fraction)"); ax.set_ylabel("Alt allele frequency, East Asian (fraction)")
    fig.text(0.54, 0.90, "Grey hexagons:\nvariant density\n(darker = more)", color=style.SUBTLE, fontsize=13, va="top")
    fig.text(0.54, 0.65, f"Points: most divergent\nvariant of each of the\ntop {N_LABEL} genes by FST", color=style.TEXT,
             fontsize=13, va="top")
    margin_key(fig, 0.54, 0.38, [("CPIC level A", style.LEVEL_A), ("other CPIC genes", style.OTHER_CPIC)], dy=0.075)
    style.deck_save(fig, "f2_allele_freq_scatter")


def fig_ecdf_deck(g, p, bg):
    groups = [(bg, style.SUBTLE, (0, (4, 2)), "random non-pharmacogenes", "right", "bottom", 1.0),
              (g[~g.level_A].fst.to_numpy(), style.OTHER_CPIC, "-", "other CPIC genes", "left", "bottom", 1.0),
              (g[g.level_A].fst.to_numpy(), style.LEVEL_A, "-", "CPIC level A genes", "right", "top", .93)]
    fig, ax = plt.subplots(figsize=(8.9, 4.55))
    fig.subplots_adjust(left=0.095, right=0.97, top=0.965, bottom=0.17)
    for x, c, ls, lab, ha, va, ty in groups:
        xs, ys = ecdf(x)
        ax.step(xs, ys, where="post", color=c, ls=ls, lw=3.6 if ls == "-" else 3)
        tx = xs[-1] + (.005 if ha == "left" else -.005)
        ax.text(tx, ty + (.02 if va == "bottom" else 0), f"{lab}\n(n={len(x)})", color=c, fontsize=13,
                fontweight="bold", ha=ha, va=va)
    ax.set_ylim(0, 1.2); ax.set_xlim(0, 0.43)
    ax.set_yticks([0, .2, .4, .6, .8, 1.0])
    ax.set_xlabel("Gene-level Hudson FST, East vs South Asian (dimensionless)")
    ax.set_ylabel("Cumulative fraction of genes")
    ax.text(0.43, 0.30, "Level A vs other CPIC genes:\none-sided permutation test\non mean gene FST, "
            f"p = {p:.3f}", color=style.TEXT, fontsize=13, ha="right", va="top")
    style.deck_save(fig, "f3_levelA_vs_rest")


def figures(d, g, p):
    plt.rcParams.update(style.RC)
    bg = load_background()
    if style.DECK:
        fig_ranking_deck(g, bg)
        fig_scatter_deck(d, g)
        fig_ecdf_deck(g, p, bg)
        return
    fig_ranking(g, bg)
    fig_scatter(d, g)
    fig_ecdf(g, p, bg)


def main():
    OUT.mkdir(parents=True, exist_ok=True); FIG.mkdir(exist_ok=True)
    d = load()
    g = gene_fst(d)
    if not style.DECK:
        g.to_csv(OUT / "gene_divergence.csv")
    obs, p = perm_test(g)
    print(f"genes kept (>= {MIN_VAR} common variants): level A {g.level_A.sum()}, other {(~g.level_A).sum()}")
    print(f"median gene FST: level A {g[g.level_A].fst.median():.4f}, other {g[~g.level_A].fst.median():.4f}")
    print(f"mean diff {obs:.4f}, permutation p={p:.4f}")
    print(g.head(10)[["fst", "fst_drop1", "n_var", "level_A"]])
    figures(d, g, p)


if __name__ == "__main__":
    main()
