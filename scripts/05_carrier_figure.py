"""Figure 4: share of each ancestry group carrying a clinically relevant variant, by gene-drug pair."""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

sys.path.insert(0, str(Path(__file__).parent))
import style  # noqa: E402

OUT, FIG = Path("data/processed"), Path("figures")
POPS = ("eas", "sas", "nfe")
GROUP_GAP = 0.9           # extra vertical space between the EAS-higher and SAS-higher blocks
BRACKET_X, BRACKET_LABEL_X, GAP_X = -0.47, -0.53, 1.12  # axes-fraction positions of the margin furniture

plt.rcParams.update(style.RC)
q = pd.read_csv(OUT / "carrier_rates.csv")
q = q[q["pop"].isin(POPS) & (q.gene != "CYP3A5")].copy()  # CYP3A5 *3 is the majority allele; "carrier" misleads
q["row"] = q.gene + "  (" + q.drug + ")"
wide = q.pivot(index="row", columns="pop", values="carrier") * 100
wide["gap"] = wide.eas - wide.sas
eas_hi = wide[wide.gap >= 0].sort_values("gap")   # largest gap at the top of each block
sas_hi = wide[wide.gap < 0].sort_values("gap", ascending=False)
wide = pd.concat([sas_hi, eas_hi])
n_sas = len(sas_hi)
y = [k + (GROUP_GAP if k >= n_sas else 0) for k in range(len(wide))]
direction = [style.POP["eas"] if g >= 0 else style.POP["sas"] for g in wide.gap]

if style.DECK:  # slide-sized variant: no embedded title, subtitle or footnote; larger text
    fig, ax = plt.subplots(figsize=(8.9, 4.55))
    fig.subplots_adjust(left=0.345, right=0.835, top=0.86, bottom=0.17)
    ax.hlines(y, wide[["eas", "sas"]].min(axis=1), wide[["eas", "sas"]].max(axis=1), color=direction, lw=7, zorder=1)
    ax.scatter(wide.eas, y, s=210, color=style.POP["eas"], edgecolor="white", linewidth=1.5, zorder=3)
    ax.scatter(wide.sas, y, s=210, color=style.POP["sas"], edgecolor="white", linewidth=1.5, zorder=3)
    ax.scatter(wide.nfe, y, s=80, facecolor="white", edgecolor=style.POP["nfe"], linewidth=2.4, zorder=4)
    top = y[-1]
    for pop, lab, dx, ha in (("sas", "South Asian", 4, "right"), ("nfe", "European\n(reference)", 22, "right"),
                             ("eas", "East Asian", 6, "right")):
        ax.annotate(lab, (wide[pop].iloc[-1], top), xytext=(dx, 12), textcoords="offset points", ha=ha,
                    fontsize=13, fontweight="bold", color=style.POP[pop])
    ax.set_yticks(y); ax.set_yticklabels(wide.index, fontsize=13)
    ax.set_ylim(-.6, top + 1.9)
    ax.set_xlim(-2, 102); ax.set_xticks(range(0, 101, 20))
    ax.set_xlabel("People carrying at least one variant (%)")
    ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
    ax.grid(axis="x", color="#EEEEEE", lw=.8, zorder=0)
    yt = ax.get_yaxis_transform()
    bx, bl = -0.585, -0.655  # left-margin brackets, in axes fractions (axes are about 4.4 in wide)
    for lo, hi, lab, col in ((0, n_sas - 1, "South Asian\nhigher", style.POP["sas"]),
                             (n_sas, len(y) - 1, "East Asian\nhigher", style.POP["eas"])):
        y0, y1 = y[lo] - .32, y[hi] + .32
        ax.plot([bx + .015, bx, bx, bx + .015], [y0, y0, y1, y1], color=col, lw=2.2, transform=yt, clip_on=False)
        ax.text(bl, (y0 + y1) / 2, lab, rotation=90, ha="center", va="center", color=col, linespacing=1.1,
                fontsize=13, fontweight="bold", transform=yt, clip_on=False)
    ax.text(1.065, top + 1.25, "Gap, East Asian\nminus South Asian\n(percentage points)", ha="center", va="bottom",
            fontsize=12.5, color=style.SUBTLE, transform=yt, clip_on=False)
    for yi, g, c in zip(y, wide.gap, direction):
        ax.text(1.065, yi, f"{g:+.0f}".replace("-", "\u2212"), ha="center", va="center", fontsize=15, fontweight="bold",
                color=c, transform=yt, clip_on=False)
    style.deck_save(fig, "f4_carrier_rates")
    raise SystemExit(0)

fig, ax = plt.subplots(figsize=(10, 5.8))
fig.subplots_adjust(left=0.33, right=0.84, top=0.80, bottom=0.24)
ax.hlines(y, wide[["eas", "sas"]].min(axis=1), wide[["eas", "sas"]].max(axis=1), color=direction, lw=6, zorder=1)
ax.scatter(wide.eas, y, s=170, color=style.POP["eas"], edgecolor="white", linewidth=1.5, zorder=3)
ax.scatter(wide.sas, y, s=170, color=style.POP["sas"], edgecolor="white", linewidth=1.5, zorder=3)
ax.scatter(wide.nfe, y, s=60, facecolor="white", edgecolor=style.POP["nfe"], linewidth=2, zorder=4)  # reference group

top = y[-1]
for pop, lab, dx, ha in (("sas", "South Asian", 0, "center"), ("nfe", "European\n(reference)", 0, "center"),
                         ("eas", "East Asian", 8, "right")):
    ax.annotate(lab, (wide[pop].iloc[-1], top), xytext=(dx, 14), textcoords="offset points", ha=ha,
                fontsize=style.KEY_FONT, fontweight="bold", color=style.POP[pop])
ax.set_yticks(y); ax.set_yticklabels(wide.index)
ax.set_ylim(-.6, top + 1.2)
ax.set_xlim(-2, 102); ax.set_xticks(range(0, 101, 20))
ax.set_xlabel("People carrying at least one variant (%)")
ax.spines["left"].set_visible(False); ax.tick_params(axis="y", length=0)
ax.grid(axis="x", color="#EEEEEE", lw=.8, zorder=0)

yt = ax.get_yaxis_transform()
for lo, hi, lab, col in ((0, n_sas - 1, "South Asian\nhigher", style.POP["sas"]),
                         (n_sas, len(y) - 1, "East Asian\nhigher", style.POP["eas"])):  # left-margin block brackets
    y0, y1 = y[lo] - .32, y[hi] + .32
    ax.plot([BRACKET_X + .012, BRACKET_X, BRACKET_X, BRACKET_X + .012], [y0, y0, y1, y1], color=col, lw=2,
            transform=yt, clip_on=False)
    ax.text(BRACKET_LABEL_X, (y0 + y1) / 2, lab, rotation=90, ha="center", va="center", color=col, linespacing=1.1,
            fontsize=style.KEY_FONT, fontweight="bold", transform=yt, clip_on=False)

ax.text(GAP_X, top + .75, "Gap, East Asian\nminus South Asian\n(percentage points)", ha="center", va="bottom",
        fontsize=style.KEY_FONT, color=style.SUBTLE, transform=yt, clip_on=False)
for yi, g, c in zip(y, wide.gap, direction):  # right-margin column, aligned to each row
    ax.text(GAP_X, yi, f"{g:+.0f}".replace("-", "−"), ha="center", va="center", fontsize=13, fontweight="bold",
            color=c, transform=yt, clip_on=False)

fig.text(0.02, 0.965, "Carrier frequency of pharmacogenetic variants by ancestry group", fontsize=15,
         fontweight="bold", va="top")
fig.text(0.02, 0.905, "% of people with at least one copy of the listed variants, per gene-drug pair (CPIC alleles; gnomAD v4)",
         fontsize=style.KEY_FONT, color=style.SUBTLE, va="top")
fig.text(0.02, 0.02, "Carrier = at least one copy of the listed core variants (Hardy-Weinberg).\nVKORC1: dose-lowering variant, not loss of function. "
         "UGT1A1 *28 via tag SNP. No Malay data.", fontsize=style.KEY_FONT, color=style.SUBTLE)
fig.savefig(FIG / "f4_carrier_rates.svg"); fig.savefig(FIG / "f4_carrier_rates.png", dpi=220)
