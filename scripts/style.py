"""Shared visual vocabulary: one colour, one meaning, across every figure.

Colours are Okabe-Ito based (colour-blind safe). Level A has its own hue that no
population uses, so that colour never carries two meanings.
"""
LEVEL_A = "#6A3D9A"        # CPIC level A genes (purple)
OTHER_CPIC = "#8C8C8C"     # other CPIC genes (mid grey)
BACKGROUND = "#D0D0D0"     # random non-pharmacogene background (light grey)
POP = {"eas": "#0072B2", "sas": "#D55E00", "nfe": "#009E73"}  # East Asian, South Asian, European
POP_OTHER = "#B0B0B0"      # all other gnomAD groups when shown for context
POP_LABEL = {"eas": "East Asian", "sas": "South Asian", "nfe": "European",
             "afr": "African", "amr": "Admixed American", "asj": "Ashkenazi Jewish",
             "fin": "Finnish", "mid": "Middle Eastern"}
TEXT, SUBTLE = "#222222", "#666666"
# sized for projection on 16:9 slides: nothing below 12 pt at final size
RC = {"font.family": "sans-serif", "font.size": 12, "axes.labelsize": 14, "axes.titlesize": 15,
      "xtick.labelsize": 12, "ytick.labelsize": 12, "legend.fontsize": 12,
      "axes.spines.top": False, "axes.spines.right": False, "axes.edgecolor": "#444444",
      "axes.titleweight": "bold", "lines.linewidth": 2.5}
KEY_FONT = 12              # annotation and key text


# --- optional deck mode (SDV_DECK=1): slide-sized variants for presentation/assets/figures/.
# Default behaviour and default outputs are unchanged. In deck mode the figure code draws without
# embedded title or footnote, at the size it appears on the slide, with larger text.
import os
from pathlib import Path

DECK = os.environ.get("SDV_DECK") == "1"
DECK_DIR = Path("presentation/assets/figures")
if DECK:
    KEY_FONT = 13
    RC.update({"font.size": 13, "axes.labelsize": 14, "xtick.labelsize": 13, "ytick.labelsize": 13,
               "legend.fontsize": 13, "axes.titlesize": 14, "lines.linewidth": 3})


def deck_save(fig, name):
    """Write a deck variant (PNG, 300 dpi) and close the figure."""
    import matplotlib.pyplot as plt
    DECK_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(DECK_DIR / f"{name}.png", dpi=300)
    plt.close(fig)
