"""Population structure: PCA of allele frequencies, pharmacogenes vs random-gene control.

Usage: python scripts/07_population_pca.py [pharmacogene_csv] [background_csv]

Matrix: 8 gnomAD groups x variants (allele frequency = ac/an). Variants need an >= MIN_AN in
every group and no missing data. Each variant is centred across groups and divided by
sqrt(pbar(1-pbar)), pbar = mean frequency over groups (Patterson et al. 2006 style). Scaling
stops the few common, high-variance variants dominating, and puts variants on a common
footing so distances can be compared between variant sets of different size and composition.
Centring is across groups; with 8 groups there are at most 7 components, so PCA here can
only show the dominant axes of group-to-group frequency difference. It cannot show
sub-structure within groups, cluster membership, or anything about individuals. Variants are
not independent (linkage disequilibrium), so variant counts overstate the information.

Control: (a) Procrustes disparity and correlation of pairwise group distances between the
pharmacogene and background layouts; (b) both sets are subsampled to the same K_DRAW variants
N_DRAWS times (seed 0; the background has fewer variants than the pharmacogene set, so it
cannot be matched upward) to ask whether pharmacogenes separate groups more or
less than the same number of random variants (RMS pairwise group distance per variant, and
share of variance on PC1/PC2).
"""
import sys
import textwrap
from itertools import combinations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
from style import DECK, POP, POP_LABEL, POP_OTHER, RC, SUBTLE, deck_save

RAW, OUT, FIG = Path("data/raw"), Path("data/processed"), Path("figures")
POPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]
MIN_AN, N_DRAWS, SEED, K_DRAW = 1000, 100, 0, 1000
# candidate label offsets in points (dx, dy, ha, va), tried in order until one is free
CAND = [(7, 5, "left", "bottom"), (7, -5, "left", "top"), (-7, 5, "right", "bottom"), (-7, -5, "right", "top"),
        (0, 12, "center", "bottom"), (0, -12, "center", "top"), (10, 0, "left", "center"), (-10, 0, "right", "center"),
        (16, 16, "left", "bottom"), (16, -16, "left", "top"), (-16, 16, "right", "bottom"), (-16, -16, "right", "top"),
        (28, 22, "left", "bottom"), (28, -22, "left", "top"), (-28, 22, "right", "bottom"), (-28, -22, "right", "top"),
        (40, 0, "left", "center"), (-40, 0, "right", "center"), (0, 30, "center", "bottom"), (0, -30, "center", "top")]


def load_matrix(path):
    """Groups x variants frequency matrix after the AN and missing-data filters."""
    d = pd.read_csv(path)
    ac = d[["ac_" + p for p in POPS]].to_numpy(float)
    an = d[["an_" + p for p in POPS]].to_numpy(float)
    ok = (an >= MIN_AN).all(axis=1) & ~np.isnan(ac / an).any(axis=1)
    return (ac[ok] / an[ok]).T, d.loc[ok, "gene"].nunique()


def scale(freq):
    """Centre each variant across groups and divide by sqrt(pbar(1-pbar))."""
    pbar = freq.mean(axis=0)
    keep = (pbar > 0) & (pbar < 1)
    return (freq[:, keep] - pbar[keep]) / np.sqrt(pbar[keep] * (1 - pbar[keep]))


def pca(x):
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    var = s ** 2 / (s ** 2).sum()
    return (u * s)[:, :7], var[:7]


def pdist(x):
    return np.array([np.linalg.norm(x[i] - x[j]) for i, j in combinations(range(len(x)), 2)])


def procrustes(ref, mov):
    """Rotate/reflect/scale mov onto ref (both centred). Returns aligned mov and disparity."""
    a, b = ref - ref.mean(0), mov - mov.mean(0)
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    u, s, vt = np.linalg.svd(b.T @ a)
    aligned = b @ (u @ vt) * (na / nb)
    return aligned, 1 - (s.sum() / (na * nb)) ** 2


def rms_dist(x):
    return np.sqrt((pdist(x) ** 2).mean() / x.shape[1])  # per-variant, so sizes compare


pharma_f, bg_f = (sys.argv[1] if len(sys.argv) > 1 else RAW / "allpop_pharmacogenes.csv",
                  sys.argv[2] if len(sys.argv) > 2 else RAW / "allpop_background.csv")
fp, ngp = load_matrix(pharma_f)
fb, ngb = load_matrix(bg_f)
xp, xb = scale(fp), scale(fb)
sp, vp = pca(xp)
sb, vb = pca(xb)
print(f"pharmacogenes: {xp.shape[1]} variants, {ngp} genes; background: {xb.shape[1]} variants, {ngb} genes")
print("pharmacogene variance explained:", np.round(vp * 100, 1))
print("background variance explained:  ", np.round(vb * 100, 1))

# layout agreement (first two PCs for Procrustes, all groups' pairwise distances in the full space)
bg_aligned, disp = procrustes(sp[:, :2], sb[:, :2])
dist_r = np.corrcoef(pdist(sp), pdist(sb))[0, 1]
print(f"Procrustes disparity (PC1-2) {disp:.3f}; pairwise-distance correlation {dist_r:.3f}")

if DECK:  # slide-sized variant: no embedded title or footnote, larger text; writes no CSV
    plt.rcParams.update(RC)
    REF = POP["nfe"]
    fig, axes = plt.subplots(1, 2, figsize=(8.9, 4.55), sharex=True, sharey=True, gridspec_kw={"wspace": 0.34})
    panels = [(axes[0], sp[:, :2], f"Pharmacogenes\n{xp.shape[1]:,} variants, {ngp} genes"),
              (axes[1], bg_aligned, f"Random non-pharmacogenes\n{xb.shape[1]:,} variants, {ngb} genes")]
    allc = np.vstack([sp[:, :2], bg_aligned])
    pad = 0.12 * np.ptp(allc, axis=0)
    axes[0].set_xlim(allc[:, 0].min() - pad[0], allc[:, 0].max() + pad[0] + 4.5)
    axes[0].set_ylim(allc[:, 1].min() - pad[1] - 2, allc[:, 1].max() + pad[1] + 6)
    for ax, c, title in panels:
        ax.set_title(title, loc="left", fontsize=13)
        ax.set_aspect("equal", adjustable="box")
        for i, p in enumerate(POPS):
            if p in ("eas", "sas"):
                ax.scatter(*c[i], s=240, color=POP[p], edgecolor="white", linewidth=1.5, zorder=4)
            elif p == "nfe":
                ax.scatter(*c[i], s=100, facecolor="white", edgecolor=REF, linewidth=2.6, zorder=3)
            else:
                ax.scatter(*c[i], s=80, color=POP_OTHER, edgecolor="white", zorder=2)
    fig.subplots_adjust(left=0.085, right=0.985, top=0.83, bottom=0.17)
    # direct labels for the four key groups only, at fixed offsets (no leader lines); other groups stay grey dots
    KEY = {"afr": (-2, 14, "left"), "eas": (0, 16, "center"), "sas": (0, -22, "center"), "nfe": (0, 16, "center")}
    for ax, c, _ in panels:
        for i, p in enumerate(POPS):
            if p in KEY:
                dx, dy, ha = KEY[p]
                ax.annotate(POP_LABEL[p], c[i], xytext=(dx, dy), textcoords="offset points", ha=ha, va="center",
                            fontsize=13, color=POP.get(p, SUBTLE), fontweight="bold" if p in POP else "normal")
    axes[0].set_ylabel(f"PC2 ({vp[1] * 100:.0f}% of variance)")
    axes[0].set_xlabel(f"PC1 ({vp[0] * 100:.0f}% of variance)")
    axes[1].set_ylabel(f"Axis 2 (own PC2: {vb[1] * 100:.0f}%)")
    axes[1].set_xlabel(f"Axis 1 (own PC1: {vb[0] * 100:.0f}%)")
    axes[1].tick_params(labelleft=True)
    deck_save(fig, "f6_population_pca")
    raise SystemExit(0)

# subsample BOTH sets to the same variant count K (the background has fewer variants than the
# pharmacogene set, so matching upward is impossible); compare the two distributions of draws
k = min(K_DRAW, xp.shape[1], xb.shape[1])
rng = np.random.default_rng(SEED)


def draw_metrics(freq, ref):
    rows = []
    for i in range(N_DRAWS):
        xs = scale(freq[:, rng.choice(freq.shape[1], size=k, replace=False)])
        ss, vs = pca(xs)
        rows.append({"rms_dist": rms_dist(xs), "pc1": vs[0], "pc12": vs[:2].sum(),
                     "procrustes_vs_pharma_full": procrustes(ref[:, :2], ss[:, :2])[1],
                     "dist_corr_vs_pharma_full": np.corrcoef(pdist(ref), pdist(ss))[0, 1]})
    return pd.DataFrame(rows)


dp, db = draw_metrics(fp, sp), draw_metrics(fb, sp)
summary = []
for m in dp.columns:
    summary.append({"metric": m, "pharma_draw_mean": dp[m].mean(), "pharma_p2.5": dp[m].quantile(.025),
                    "pharma_p97.5": dp[m].quantile(.975), "random_draw_mean": db[m].mean(),
                    "random_p2.5": db[m].quantile(.025), "random_p97.5": db[m].quantile(.975),
                    "frac_random_draws_above_pharma_mean": (db[m] > dp[m].mean()).mean()})
summary = pd.DataFrame(summary)
extra = {"procrustes_full_background": disp, "dist_corr_full_background": dist_r,
         "n_variants_pharmacogene": xp.shape[1], "n_variants_background": xb.shape[1],
         "n_variants_per_draw": k, "n_genes_pharmacogene": ngp, "n_genes_background": ngb}
summary = pd.concat([summary, pd.DataFrame({"metric": list(extra), "pharma_draw_mean": list(extra.values())})])
print(summary.round(4).to_string(index=False))

OUT.mkdir(exist_ok=True)
pcs = [f"PC{i + 1}" for i in range(7)]
scores = pd.concat([pd.DataFrame(sp, columns=pcs).assign(pop=POPS, set="pharmacogenes"),
                    pd.DataFrame(sb, columns=pcs).assign(pop=POPS, set="background")])
scores.to_csv(OUT / "pca_scores.csv", index=False)
pd.DataFrame({"pc": pcs, "pharmacogenes": vp, "background": vb}).to_csv(OUT / "pca_variance.csv", index=False)
summary.to_csv(OUT / "pca_control_summary.csv", index=False)

# figure
plt.rcParams.update(RC)
REF = POP["nfe"]  # European: small hollow ring, a reference rather than a focus group
fig, axes = plt.subplots(1, 2, figsize=(10, 5.9), sharex=True, sharey=True, gridspec_kw={"wspace": 0.12})
panels = [(axes[0], sp[:, :2], f"Pharmacogenes\n{xp.shape[1]:,} variants, {ngp} genes"),
          (axes[1], bg_aligned, f"Random non-pharmacogenes\n{xb.shape[1]:,} variants, {ngb} genes")]
allc = np.vstack([sp[:, :2], bg_aligned])
pad = 0.12 * np.ptp(allc, axis=0)
axes[0].set_xlim(allc[:, 0].min() - pad[0], allc[:, 0].max() + pad[0] + 2.5)
axes[0].set_ylim(allc[:, 1].min() - pad[1], allc[:, 1].max() + pad[1])
for ax, c, title in panels:
    ax.set_title(title, loc="left", fontsize=13)
    ax.set_aspect("equal", adjustable="box")
    for i, p in enumerate(POPS):
        if p in ("eas", "sas"):
            ax.scatter(*c[i], s=200, color=POP[p], edgecolor="white", linewidth=1.5, zorder=4)
        elif p == "nfe":
            ax.scatter(*c[i], s=80, facecolor="white", edgecolor=REF, linewidth=2.2, zorder=3)
        else:
            ax.scatter(*c[i], s=60, color=POP_OTHER, edgecolor="white", zorder=2)
fig.subplots_adjust(left=0.08, right=0.98, top=0.80, bottom=0.27)
fig.canvas.draw()
rend = fig.canvas.get_renderer()
for ax, c, _ in panels:  # greedy label placement avoiding other labels and points
    boxes = [ax.transData.transform(pt) for pt in c]
    taken = [(x - 10, y - 10, x + 10, y + 10) for x, y in boxes]
    for i in np.argsort(-np.abs(c[:, 0] - c[:, 0].mean())):  # outlying groups first, crowded ones pick what is left
        p = POPS[i]
        for dx, dy, ha, va in CAND:
            t = ax.annotate(POP_LABEL[p], c[i], xytext=(dx, dy), textcoords="offset points", ha=ha, va=va,
                            fontsize=12, color=POP.get(p, SUBTLE), fontweight="bold" if p in POP else "normal")
            bb = t.get_window_extent(rend)
            inside = ax.get_window_extent()
            free = all(bb.x1 < q[0] or bb.x0 > q[2] or bb.y1 < q[1] or bb.y0 > q[3] for q in taken)
            if free and bb.x0 > inside.x0 and bb.x1 < inside.x1 and bb.y0 > inside.y0 and bb.y1 < inside.y1:
                taken.append((bb.x0, bb.y0, bb.x1, bb.y1))
                if abs(dx) > 8 or abs(dy) > 8:  # leader line only where the label sits away from its point
                    ax.annotate("", c[i], xytext=(dx, dy), textcoords="offset points",
                                arrowprops={"arrowstyle": "-", "color": "#999999", "lw": 0.8, "shrinkA": 0, "shrinkB": 5})
                break
            t.remove()
        else:
            print("label placement failed for", p)
            ax.annotate(POP_LABEL[p], c[i], xytext=CAND[0][:2], textcoords="offset points", fontsize=12)
axes[0].set_ylabel(f"PC2 ({vp[1] * 100:.0f}% of variance)")
axes[0].set_xlabel(f"PC1 ({vp[0] * 100:.0f}% of variance)")
axes[1].set_ylabel(f"Axis 2 (own PC2: {vb[1] * 100:.0f}%)")
axes[1].set_xlabel(f"Axis 1 (own PC1: {vb[0] * 100:.0f}% of variance)")
axes[1].tick_params(labelleft=True)
fig.text(0.02, 0.975, "PCA of common-variant allele frequencies across gnomAD v4 groups", fontsize=15, fontweight="bold", va="top")
note = ("Each point is one gnomAD group; distance reflects how differently groups' allele frequencies behave across the "
        "plotted variants, not divergence time or clinical difference. Right: random non-pharmacogenes, aligned to the left "
        f"layout (Procrustes disparity {disp:.2f}; distance correlation {dist_r:.2f}). Variants are linked; "
        "8 groups give at most 7 components.")
fig.text(0.02, 0.012, textwrap.fill(note, 112), fontsize=12, color=SUBTLE, va="bottom")
fig.savefig(FIG / "f6_population_pca.svg")
fig.savefig(FIG / "f6_population_pca.png", dpi=220)
