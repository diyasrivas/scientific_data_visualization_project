"""Per-variant East vs South Asian test, Benjamini-Hochberg correction, and the volcano figure F7.

Figure specification (question first, then one line per encoding):

    Which common variants differ enough between the two groups to matter, and does statistical
    significance help us choose them?
    x: allele-frequency difference, East Asian minus South Asian (signed, zero meaningful).
    y: -log10 of the BH-adjusted p-value, on a log scale because it spans 0 to 12,000.
       Nothing is clipped; the log scale is what makes the saturation visible instead of hiding it.
    colour: which group carries the variant more often; grey inside the small-effect band.
      Colour is redundant with the sign of x, so it adds robustness without a conjunction.

Chi-square is computed on allele counts, so the two chromosomes of one person are treated as
independent observations. That inflates significance, which is one more reason the figure ranks
variants by effect rather than by p.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.special import log_ndtr

sys.path.insert(0, str(Path(__file__).parent))
import style  # noqa: E402

RAW, OUT, FIG = Path("data/raw"), Path("data/processed"), Path("figures")
MIN_AN, MIN_MAF = 1000, 0.01     # same filters as 03_analyse.py, so the two figures describe one variant set
ALPHA = 0.05
BIG_EFFECT = 0.20                # |allele-frequency difference| below this is "too small to matter here"
Y_FLOOR = 0.03                   # bottom of the log y axis; the least significant variant is at 0.01
N_LABEL = 7                      # genes labelled, one variant each (the most divergent of that gene)
LABEL_DY = (0, -15, 15, -30, 30, -46, 46, -62, 62)   # offsets tried in order when labels collide
LN10 = np.log(10)


def load():
    d = pd.read_csv(RAW / "allpop_pharmacogenes.csv")[["gene", "variant_id", "ac_eas", "an_eas", "ac_sas", "an_sas"]]
    d = d[(d.an_eas >= MIN_AN) & (d.an_sas >= MIN_AN)].copy()
    d["p_eas"], d["p_sas"] = d.ac_eas / d.an_eas, d.ac_sas / d.an_sas
    maf = np.maximum(np.minimum(d.p_eas, 1 - d.p_eas), np.minimum(d.p_sas, 1 - d.p_sas))
    d = d[maf >= MIN_MAF].copy()
    d["delta"] = d.p_eas - d.p_sas
    return d


def log10_p(d):
    """log10 of the two-sided p-value from a 2x2 chi-square on allele counts.

    Computed through the normal tail, because chi2.logsf underflows to -inf well inside the range
    gnomAD's sample sizes produce (the most extreme variant here is below 1e-1000).
    """
    a, b = d.ac_eas.to_numpy(float), (d.an_eas - d.ac_eas).to_numpy(float)
    c, e = d.ac_sas.to_numpy(float), (d.an_sas - d.ac_sas).to_numpy(float)
    n = a + b + c + e
    chi2 = n * (a * e - b * c) ** 2 / ((a + b) * (c + e) * (a + c) * (b + e))
    return (np.log(2) + log_ndtr(-np.sqrt(chi2))) / LN10


def bh(log10p):
    """Benjamini-Hochberg adjusted p-values, returned as log10(q).

    Done in log space because most of these p-values underflow float64. BH is the running minimum
    of p_(j) * m / j taken from the largest p downwards, which in logs is a running minimum of
    log10 p_(j) + log10(m / j).
    """
    m = len(log10p)
    order = np.argsort(log10p)[::-1]               # largest p-value first
    ranks = np.arange(m, 0, -1)                    # j = m, m-1, ... 1
    running = np.minimum.accumulate(log10p[order] + np.log10(m / ranks))
    out = np.empty(m)
    out[order] = np.minimum(running, 0.0)          # q <= 1
    return out


def analyse():
    d = load()
    d["log10_p"] = log10_p(d)
    d["log10_q"] = bh(d.log10_p.to_numpy())
    d["significant"] = d.log10_q < np.log10(ALPHA)
    d["large_effect"] = d.delta.abs() >= BIG_EFFECT
    d.sort_values("delta", key=abs, ascending=False).to_csv(OUT / "variant_tests.csv", index=False)
    n, n_sig, n_both = len(d), int(d.significant.sum()), int((d.significant & d.large_effect).sum())
    print(f"{n} common variants tested; {n_sig} with q < {ALPHA} ({100 * n_sig / n:.1f}%); "
          f"{n_both} of those also |difference| >= {BIG_EFFECT}")
    print(f"smallest q: 10^{d.log10_q.min():.0f}; median -log10(q): {-d.log10_q.median():.0f}")
    return d


def _labels(d):
    """The most divergent variant of each of the N_LABEL most divergent genes."""
    seen, out = set(), []
    for _, r in d.reindex(d.delta.abs().sort_values(ascending=False).index).iterrows():
        if r.gene not in seen:
            seen.add(r.gene)
            out.append(r)
        if len(out) == N_LABEL:
            return out
    return out


def draw(d, deck):
    """One figure, two sizes. Deck mode drops the embedded title and footnote and enlarges text."""
    plt.rcParams.update(style.RC)
    fig, ax = plt.subplots(figsize=(8.9, 4.55) if deck else (10, 6.2))
    fig.subplots_adjust(left=.135 if deck else .10, right=.985, top=.845 if deck else .80,
                        bottom=.175 if deck else .20)
    key = 13 if deck else style.KEY_FONT

    y = np.maximum(-d.log10_q, Y_FLOOR)
    colour = np.where(d.large_effect, np.where(d.delta > 0, style.POP["eas"], style.POP["sas"]), style.BACKGROUND)
    sig_y = -np.log10(ALPHA)

    n, n_sig, n_both = len(d), int(d.significant.sum()), int((d.significant & d.large_effect).sum())
    ax.set_yscale("log")
    ax.set_xlim(-.80, .95); ax.set_ylim(Y_FLOOR, 10 ** 4.9)   # before any label is placed, so the
    ax.set_xticks(np.arange(-.75, .76, .25))                  # collision boxes are in final coordinates
    ax.set_yticks([0.1, 1, 10, 1e2, 1e3, 1e4])
    ax.set_yticklabels(["0.1", "1", "10", "100", "1,000", "10,000"])
    ax.set_xlabel("Allele-frequency difference, East Asian minus South Asian")
    ax.set_ylabel("−log₁₀(q), log scale")
    ax.axvspan(-BIG_EFFECT, BIG_EFFECT, color="#F4F4F4", lw=0, zorder=0)
    ax.axhline(sig_y, color=style.SUBTLE, lw=1.2, ls=":", zorder=1)
    ax.scatter(d.delta, y, s=15 if deck else 13, c=colour, alpha=.8, linewidths=0, zorder=2)

    ax.text(.93, sig_y * .72, "q = 0.05: everything above this line is “significant”", ha="right",
            va="top", fontsize=key - 1, color=style.SUBTLE, zorder=3)
    ax.text(0, Y_FLOOR * 1.6, "difference too small\nto matter here", ha="center", va="bottom",
            fontsize=key - 1, color=style.SUBTLE, linespacing=1.15, zorder=3)

    heads = [ax.text(x, 1.015, f"◀ more common in South Asians" if p == "sas" else
                     "more common in East Asians ▶", transform=ax.transAxes, ha=ha, va="bottom",
                     fontsize=key, fontweight="bold", color=style.POP[p])
             for x, ha, p in ((0.0, "left", "sas"), (1.0, "right", "eas"))]  # above the axes, clear of the data

    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    labelled = _labels(d)
    taken = [tuple(h.get_window_extent(rend).extents) for h in heads]
    for r in labelled:
        x, y0 = ax.transData.transform((r.delta, -r.log10_q))
        taken.append((x - 7, y0 - 7, x + 7, y0 + 7))
    for r in labelled:  # first offset whose box is free of the points and of the labels already placed
        east = r.delta > 0
        c = style.POP["eas" if east else "sas"]
        for dy in LABEL_DY:
            # a leader line only when the label had to move far from its point (proximity needs help)
            arrow = dict(arrowstyle="-", lw=.8, color=c, shrinkA=0, shrinkB=4) if abs(dy) > 20 else None
            t = ax.annotate(r.gene, (r.delta, -r.log10_q), xytext=(9 if east else -9, dy),
                            textcoords="offset points", ha="left" if east else "right", va="center",
                            fontsize=key - 1, fontweight="bold", color=c, arrowprops=arrow)
            bb = t.get_window_extent(rend)
            if all(bb.x1 < q[0] or bb.x0 > q[2] or bb.y1 < q[1] or bb.y0 > q[3] for q in taken):
                taken.append((bb.x0, bb.y0, bb.x1, bb.y1))
                break
            t.remove()

    ax.text(-.76, 12, f"{n_sig:,} of {n:,} variants have q < {ALPHA}.\nOnly {n_both} differ by {BIG_EFFECT:g} or more.",
            ha="left", va="center", fontsize=key, color=style.TEXT, fontweight="bold", linespacing=1.5)

    if deck:
        style.deck_save(fig, "f7_volcano")
        return
    fig.text(.02, .975, "Significance saturates; effect size does not", fontsize=15, fontweight="bold", va="top")
    fig.text(.02, .930, f"Every common variant in {d.gene.nunique()} CPIC genes, East Asian vs South Asian\n"
             "(gnomAD v4; chi-square on allele counts, Benjamini-Hochberg over all tests)",
             fontsize=style.KEY_FONT, color=style.SUBTLE, va="top", linespacing=1.3)
    fig.text(.02, .035, f"{len(d):,} tests. Allele counts treat a person's two chromosomes as independent, which "
             "inflates significance.\nThe grey band is a judgement about what is worth attention, not a result.",
             fontsize=style.KEY_FONT, color=style.SUBTLE, linespacing=1.3)
    fig.savefig(FIG / "f7_volcano.svg"); fig.savefig(FIG / "f7_volcano.png", dpi=220)
    plt.close(fig)


if __name__ == "__main__":
    draw(analyse(), deck=style.DECK)
