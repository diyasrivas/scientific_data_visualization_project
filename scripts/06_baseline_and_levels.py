"""Baseline FST of random non-CPIC genes, CPIC pharmacogenes vs baseline, and FST by CPIC level."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import requests
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).parent))
import style

RAW, OUT, FIG = Path("data/raw"), Path("data/processed"), Path("figures")
SEED, MIN_AN, MIN_MAF, MIN_VAR, N_PERM = 0, 1000, 0.01, 5, 10000  # same filters as 03_analyse.py
CPIC_URL = "https://api.cpicpgx.org/v1/pair_view?select=genesymbol,cpiclevel"
RANK = {"A": 0, "B": 1, "C": 2, "D": 3}
LEVEL_MAP = {"A": "A", "A/B": "A", "B": "B", "B/C": "B", "C": "C", "C/D": "C", "D": "D"}
rng = np.random.default_rng(SEED)


def hudson_terms(p1, n1, p2, n2):
    """Per-variant numerator and denominator of Hudson's FST estimator (as in 03_analyse.py)."""
    num = (p1 - p2) ** 2 - p1 * (1 - p1) / (n1 - 1) - p2 * (1 - p2) / (n2 - 1)
    den = p1 * (1 - p2) + p2 * (1 - p1)
    return num, den


def background_fst():
    d = pd.read_csv(RAW / "allpop_background.csv")
    d = d[(d.an_eas >= MIN_AN) & (d.an_sas >= MIN_AN)].copy()
    d["p_eas"], d["p_sas"] = d.ac_eas / d.an_eas, d.ac_sas / d.an_sas
    maf = np.maximum(np.minimum(d.p_eas, 1 - d.p_eas), np.minimum(d.p_sas, 1 - d.p_sas))
    d = d[maf >= MIN_MAF].copy()
    d["num"], d["den"] = hudson_terms(d.p_eas, d.an_eas, d.p_sas, d.an_sas)
    g = d.groupby("gene").agg(num=("num", "sum"), den=("den", "sum"), n_var=("num", "size"))
    g["fst"] = g.num / g.den
    g = g[g.n_var >= MIN_VAR]
    return g[["fst", "n_var"]].reset_index()


def perm_diff(x, lab, strata=None, n=N_PERM):
    """Permutation test on mean(x[lab]) - mean(x[~lab]); labels shuffled within strata if given."""
    obs = x[lab].mean() - x[~lab].mean()
    groups = [np.arange(len(x))] if strata is None else [np.where(strata == s)[0] for s in np.unique(strata)]
    null = np.empty(n)
    for i in range(n):
        sh = lab.copy()
        for ix in groups:
            sh[ix] = rng.permutation(lab[ix])
        null[i] = x[sh].mean() - x[~sh].mean()
    return obs, (np.sum(null >= obs) + 1) / (n + 1)


def compare(bg, ph):
    x = np.concatenate([ph.fst.to_numpy(), bg.fst.to_numpy()])
    nv = np.concatenate([ph.n_var.to_numpy(), bg.n_var.to_numpy()])
    lab = np.r_[np.ones(len(ph), bool), np.zeros(len(bg), bool)]
    q = (lambda s: f"{s.quantile(.25):.4f}-{s.quantile(.75):.4f}")
    print(f"\nBackground: {len(bg)} genes, median FST {bg.fst.median():.4f} (IQR {q(bg.fst)}), "
          f"median n_var {bg.n_var.median():.0f}")
    print(f"CPIC:       {len(ph)} genes, median FST {ph.fst.median():.4f} (IQR {q(ph.fst)}), "
          f"median n_var {ph.n_var.median():.0f}")
    obs, p = perm_diff(x, lab)
    print(f"mean diff (CPIC - background) {obs:.4f}, one-sided permutation p={p:.4f}")
    print(f"Spearman FST vs n_var, pooled: rho={spearmanr(nv, x)[0]:.2f}; "
          f"background only: rho={spearmanr(bg.n_var, bg.fst)[0]:.2f}")
    print(f"n_var is {'higher' if ph.n_var.median() > bg.n_var.median() else 'lower'} in CPIC genes")
    edges = np.unique(np.quantile(nv, [0, 1 / 3, 2 / 3, 1]))
    bins = np.clip(np.digitize(nv, edges[1:-1], right=True), 0, len(edges) - 2)
    for b in range(len(edges) - 1):
        m = bins == b
        a, c = x[m & lab], x[m & ~lab]
        print(f"  n_var tertile {b + 1} ({nv[m].min()}-{nv[m].max()}): CPIC n={len(a)} median {np.median(a) if len(a) else np.nan:.4f}; "
              f"background n={len(c)} median {np.median(c) if len(c) else np.nan:.4f}")
    ok = np.array([lab[bins == b].any() and (~lab[bins == b]).any() for b in bins])
    obs_s, p_s = perm_diff(x[ok], lab[ok], strata=bins[ok])
    print(f"n_var-stratified permutation (within tertiles): mean diff {obs_s:.4f}, p={p_s:.4f}")
    lo, hi = max(ph.n_var.min(), bg.n_var.min()), min(ph.n_var.max(), bg.n_var.max())
    m = (nv >= lo) & (nv <= hi)
    obs_o, p_o = perm_diff(x[m], lab[m])
    print(f"restricted to overlapping n_var {lo}-{hi}: CPIC n={(m & lab).sum()}, bg n={(m & ~lab).sum()}, "
          f"mean diff {obs_o:.4f}, p={p_o:.4f}")


def cpic_levels():
    r = requests.get(CPIC_URL, timeout=60)
    r.raise_for_status()
    df = pd.DataFrame(r.json()).dropna(subset=["cpiclevel"])
    df["lvl"] = df.cpiclevel.map(LEVEL_MAP)  # 'Retired' and anything else maps to NaN and is ignored
    df = df.dropna(subset=["lvl"])
    df["rank"] = df.lvl.map(RANK)
    best = df.groupby("genesymbol")["rank"].min()
    return best.map({v: k for k, v in RANK.items()}).rename("level")


def by_level(ph):
    lv = cpic_levels()
    g = ph.rename(columns={"genesymbol": "gene"}).merge(lv, left_on="gene", right_index=True, how="left")
    n_miss = g.level.isna().sum()
    g = g.dropna(subset=["level"])
    g["level_rank"] = g.level.map(RANK)
    g[["gene", "level", "fst", "n_var"]].sort_values(["level", "fst"], ascending=[True, False]).to_csv(
        OUT / "gene_levels.csv", index=False)
    print(f"\nFST by CPIC level ({len(g)} genes; {n_miss} without a non-retired level):")
    for k in RANK:
        s = g[g.level == k].fst
        print(f"  {k}: n={len(s)}, median {s.median():.4f}, IQR {s.quantile(.25):.4f}-{s.quantile(.75):.4f}")
    rho, p = spearmanr(g.level_rank, g.fst)
    print(f"Spearman FST vs level rank (A=0..D=3): rho={rho:.2f}, p={p:.3f}")
    null = np.array([spearmanr(g.level_rank, rng.permutation(g.fst))[0] for _ in range(N_PERM)])
    print(f"permutation two-sided p={(np.sum(np.abs(null) >= abs(rho)) + 1) / (N_PERM + 1):.4f}")
    return g


def figure_deck(g, bg):
    """Slide-sized variant (SDV_DECK=1): no embedded title, larger text."""
    plt.rcParams.update(style.RC)
    jitter = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(8.9, 4.55))
    fig.subplots_adjust(left=0.125, right=0.775, top=0.965, bottom=0.30)
    q1, med, q3 = bg.fst.quantile([.25, .5, .75])
    ax.axhspan(q1, q3, color=style.BACKGROUND, alpha=.5, zorder=0, lw=0)
    ax.axhline(med, color=style.SUBTLE, lw=1.8, ls="--", zorder=1)
    levels = list(RANK)
    for i, k in enumerate(levels):
        s = g[g.level == k].fst.to_numpy()
        col = style.LEVEL_A if k == "A" else style.OTHER_CPIC
        ax.scatter(i + jitter.uniform(-.18, .18, len(s)), s, s=85, color=col, alpha=.85, zorder=3, lw=0)
        ax.hlines(np.median(s), i - .32, i + .32, color=style.TEXT, lw=4, zorder=4)
    ax.set_xlim(-.6, len(levels) - .4)
    ax.text(1.03, med, f"Random non-CPIC\ngenes (n={len(bg)}):\nmedian {med:.3f},\nshaded = IQR", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=13, color=style.SUBTLE)
    ax.text(1.03, .78, "Black bar =\nmedian per level", transform=ax.transAxes, ha="left", va="center",
            fontsize=13, color=style.TEXT)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"Level {k}\nn = {(g.level == k).sum()}\nmedian {g[g.level == k].fst.median():.3f}" for k in levels],
                       fontsize=12.5)
    ax.set_xlabel("CPIC clinical-annotation level (best level per gene)", labelpad=8)
    ax.set_ylabel("Gene-level Hudson FST")
    style.deck_save(fig, "f5_fst_by_level")
    plt.close(fig)


def figure(g, bg):
    if style.DECK:
        return figure_deck(g, bg)
    plt.rcParams.update(style.RC)
    jitter = np.random.default_rng(SEED)
    fig, ax = plt.subplots(figsize=(10, 5.6))
    fig.subplots_adjust(left=0.10, right=0.76, top=0.88, bottom=0.27)
    q1, med, q3 = bg.fst.quantile([.25, .5, .75])
    ax.axhspan(q1, q3, color=style.BACKGROUND, alpha=.5, zorder=0, lw=0)
    ax.axhline(med, color=style.SUBTLE, lw=1.5, ls="--", zorder=1)
    levels = list(RANK)
    for i, k in enumerate(levels):
        s = g[g.level == k].fst.to_numpy()
        col = style.LEVEL_A if k == "A" else style.OTHER_CPIC
        ax.scatter(i + jitter.uniform(-.18, .18, len(s)), s, s=60, color=col, alpha=.8, zorder=3, lw=0)
        ax.hlines(np.median(s), i - .32, i + .32, color=style.TEXT, lw=3, zorder=4)
    ax.set_xlim(-.6, len(levels) - .4)
    ax.text(1.03, med, f"Random non-CPIC\ngenes (n={len(bg)}):\nmedian {med:.3f},\nshaded = IQR", transform=ax.get_yaxis_transform(),
            ha="left", va="center", fontsize=style.KEY_FONT, color=style.SUBTLE)  # label aligned to the band it names
    ax.text(1.03, .72, "Black bar = median\nof each level", transform=ax.transAxes, ha="left", va="center",
            fontsize=style.KEY_FONT, color=style.TEXT)
    ax.set_xticks(range(len(levels)))
    ax.set_xticklabels([f"Level {k}\nn = {(g.level == k).sum()}\nmedian {g[g.level == k].fst.median():.3f}" for k in levels])
    ax.set_xlabel("CPIC clinical-annotation level (best level per gene)", labelpad=12)
    ax.set_ylabel("Gene-level Hudson FST (dimensionless)")
    fig.text(0.02, 0.975, "Gene-level East vs South Asian FST by CPIC level", fontsize=15, fontweight="bold", va="top")
    fig.savefig(FIG / "f5_fst_by_level.svg")
    fig.savefig(FIG / "f5_fst_by_level.png", dpi=220)
    plt.close(fig)


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    FIG.mkdir(exist_ok=True)
    bg = background_fst()
    bg.to_csv(OUT / "background_gene_fst.csv", index=False)
    ph = pd.read_csv(OUT / "gene_divergence.csv")
    compare(bg, ph)
    g = by_level(ph)
    figure(g, bg)


if __name__ == "__main__":
    main()
