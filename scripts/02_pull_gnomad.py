"""Pull gnomAD v4 allele counts for ALL main populations, per gene, for two gene sets.

  background    N_BG random autosomal protein-coding genes that are not CPIC genes (HGNC list)
  pharmacogenes the CPIC genes in data/raw/genes.csv

Usage: python3 scripts/02_pull_gnomad.py background|pharmacogenes
Output: data/raw/allpop_<set>.csv, one row per variant, columns ac_<pop>, an_<pop>.
Only variants with minor allele frequency >= MIN_MAF in at least one population are kept.
Cached per gene; re-running skips genes already in the output file.
"""
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

RAW = Path("data/raw")
GNOMAD = "https://gnomad.broadinstitute.org/api"
HGNC = "https://storage.googleapis.com/public-download-files/hgnc/tsv/tsv/locus_types/gene_with_protein_product.txt"
POPS = ["afr", "amr", "asj", "eas", "fin", "mid", "nfe", "sas"]
N_BG, SEED, MIN_MAF = 120, 1, 0.01
QUERY = ('{ gene(gene_symbol:"%s", reference_genome: GRCh38){ variants(dataset: gnomad_r4){ '
         'variant_id exome{populations{id ac an}} genome{populations{id ac an}} } } }')


def gene_list(which):
    if which == "pharmacogenes":
        return list(pd.read_csv(RAW / "genes.csv").genesymbol)
    hgnc = RAW / "hgnc_protein_coding.tsv"
    if not hgnc.exists():  # public HGNC protein-coding gene table, about 9 MB
        hgnc.write_bytes(requests.get(HGNC, timeout=120).content)
    h = pd.read_csv(hgnc, sep="\t", low_memory=False)
    h = h[(h.status == "Approved") & h.location.fillna("").str.match(r"^(\d|1\d|2[0-2])[pq]")]
    cpic = set(pd.read_csv(RAW / "genes.csv").genesymbol)
    syms = sorted(set(h.symbol) - cpic)
    rng = np.random.default_rng(SEED)
    return list(rng.choice(syms, size=N_BG, replace=False))


def fetch_gene(sym):
    for attempt in range(6):
        try:
            r = requests.post(GNOMAD, json={"query": QUERY % sym}, timeout=240)
        except requests.RequestException:
            time.sleep(15); continue
        if r.status_code == 200:
            gene = (r.json().get("data") or {}).get("gene")
            return gene["variants"] if gene else None
        time.sleep(20 * (attempt + 1))  # 429 rate limit
    return None


def counts(variant, pop):
    ac = an = 0
    for src in ("exome", "genome"):
        for p in (variant.get(src) or {}).get("populations", []):
            if p["id"] == pop:
                ac += p["ac"]; an += p["an"]
    return ac, an


def main(which):
    out = RAW / f"allpop_{which}.csv"
    done = pd.read_csv(out) if out.exists() else pd.DataFrame()
    have = set(done["gene"]) if len(done) else set()
    cols = ["gene", "variant_id"] + [f"{k}_{p}" for p in POPS for k in ("ac", "an")]
    rows = []
    for sym in gene_list(which):
        if sym in have:
            continue
        vs = fetch_gene(sym)
        if vs is None:
            print(f"{sym}: not available, skipped", flush=True); time.sleep(5); continue
        n = 0
        for v in vs:
            c = {p: counts(v, p) for p in POPS}
            afs = [c[p][0] / c[p][1] for p in POPS if c[p][1]]
            if afs and any(0.01 <= a <= 0.99 or 0.01 <= 1 - a <= 0.99 for a in afs) and max(min(a, 1 - a) for a in afs) >= MIN_MAF:
                rows.append([sym, v["variant_id"]] + [x for p in POPS for x in c[p]]); n += 1
        print(f"{sym}: {len(vs)} variants, {n} kept", flush=True)
        pd.concat([done, pd.DataFrame(rows, columns=cols)]).to_csv(out, index=False)
        time.sleep(5)
    print("saved", out, flush=True)


if __name__ == "__main__":
    main(sys.argv[1])
