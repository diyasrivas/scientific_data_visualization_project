"""Pull the CPIC gene list and each gene's level into data/raw/genes.csv.

Level is "A" if the gene has at least one gene-drug pair of CPIC level A (including A/B),
otherwise "other". Retired pairs are ignored.
"""
from pathlib import Path

import pandas as pd
import requests

RAW = Path("data/raw")
CPIC = "https://api.cpicpgx.org/v1"


def gene_levels():
    pairs = pd.DataFrame(requests.get(f"{CPIC}/pair_view?select=genesymbol,cpiclevel", timeout=60).json())
    pairs = pairs[pairs.cpiclevel != "Retired"]
    lv = pairs.groupby("genesymbol")["cpiclevel"].apply(lambda s: "A" if s.str.startswith("A").any() else "other")
    return lv.rename("level").reset_index()


def main():
    RAW.mkdir(parents=True, exist_ok=True)
    genes = gene_levels().sort_values(["level", "genesymbol"])
    genes.to_csv(RAW / "genes.csv", index=False)
    print(f"{len(genes)} CPIC genes, {(genes.level == 'A').sum()} level A -> {RAW / 'genes.csv'}")


if __name__ == "__main__":
    main()
