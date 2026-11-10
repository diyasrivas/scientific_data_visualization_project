"""Carrier rates of CPIC reduced-function alleles by gnomAD v4 population.

Alleles: CPIC star/named alleles defined by exactly one SNV with a clinical function label
of No function or Decreased function (loss of function, LoF). Allele frequency q per
population is the sum over a gene's single-SNV LoF alleles; carrier rate = 1-(1-q)^2
(Hardy-Weinberg, alleles treated as non-overlapping).
"""
import time
from pathlib import Path

import pandas as pd
import requests

RAW, OUT = Path("data/raw"), Path("data/processed")
CPIC = "https://api.cpicpgx.org/v1"
GNOMAD = "https://gnomad.broadinstitute.org/api"
POPS = ["eas", "sas", "nfe", "afr", "amr"]
GENES = {  # gene: drug shown on the figure
    "CYP2C19": "clopidogrel", "CYP2C9": "warfarin", "SLCO1B1": "simvastatin",
    "NUDT15": "thiopurines", "DPYD": "fluorouracil", "UGT1A1": "irinotecan", "CYP3A5": "tacrolimus",
}
# Curated core alleles: well-established, on different haplotypes, so their frequencies can be added.
# Automatic inclusion of every CPIC allele double counts alleles in LD (e.g. CYP2C19 *35 with *2).
CORE = {"CYP2C19": {"*2", "*3"}, "CYP2C9": {"*2", "*3"}, "SLCO1B1": {"*5"}, "NUDT15": {"*3"},
        "UGT1A1": {"*6", "*80+*28"}, "CYP3A5": {"*3"}}
DPYD_CORE = ("c.1905+1G>A", "c.2846A>T", "c.1679T>G", "c.1129-5923C>G")
LOF = {"No function", "Decreased function"}
BATCH = 15


def cpic(path):
    r = requests.get(f"{CPIC}/{path}", timeout=60); r.raise_for_status(); return r.json()


def single_snv_alleles():
    chrom = {g["symbol"]: g["chr"].replace("chr", "") for g in cpic("gene?select=symbol,chr") if g["chr"]}
    rows = []
    for gene in GENES:
        if gene not in chrom:
            print(f"{gene}: no chromosome in CPIC, skipped"); continue
        alleles = pd.DataFrame(cpic(f"allele?genesymbol=eq.{gene}&select=name,clinicalfunctionalstatus,definitionid"))
        if alleles.empty:
            continue
        ids = ",".join(map(str, alleles.definitionid.dropna().astype(int)))
        vals = pd.DataFrame(cpic(f"allele_location_value?alleledefinitionid=in.({ids})&select=alleledefinitionid,locationid,variantallele"))
        locs = pd.DataFrame(cpic(f"sequence_location?genesymbol=eq.{gene}&select=id,chromosomelocation,position,dbsnpid"))
        vals = vals[vals.variantallele != "R"]  # "R" = reference base at that position
        share = vals.groupby("locationid").alleledefinitionid.nunique() / len(alleles)
        vals = vals[vals.locationid.map(share) <= 0.5]  # drop background positions shared by most alleles
        d = alleles.merge(vals, left_on="definitionid", right_on="alleledefinitionid").merge(locs, left_on="locationid", right_on="id")
        m = d.chromosomelocation.str.extract(r"g\.(\d+)([ACGT])>([ACGT])$")
        d["ref"], d["alt"] = m[1], m[2]
        d = d.dropna(subset=["ref"])  # ignore non-SNV positions (indels, repeats): not readable from SNV frequencies
        d = d[d.variantallele == d.alt]  # allele carries the alt base
        d = d[d.groupby("definitionid").locationid.transform("nunique") == 1]  # exactly one defining SNV
        d["gene"] = gene
        d["variant_id"] = chrom[gene] + "-" + d.position.astype(int).astype(str) + "-" + d.ref + "-" + d.alt
        rows.append(d[["gene", "name", "clinicalfunctionalstatus", "dbsnpid", "variant_id"]])
    return pd.concat(rows).drop_duplicates(["gene", "name"])  # variant_id may repeat across alleles; deduplicated when summing


def vkorc1_row():
    """VKORC1 -1639G>A (rs9923231): the warfarin dose-lowering variant. Not a loss-of-function allele."""
    loc = pd.DataFrame(cpic("sequence_location?genesymbol=eq.VKORC1&dbsnpid=eq.rs9923231&select=chromosomelocation,position,dbsnpid"))
    m = loc.chromosomelocation.str.extract(r"g\.(\d+)([ACGT])>([ACGT])$").iloc[0]
    vid = f"16-{int(loc.position.iloc[0])}-{m[1]}-{m[2]}"
    return pd.DataFrame([{"gene": "VKORC1", "name": "-1639G>A", "clinicalfunctionalstatus": "Lower warfarin dose",
                          "dbsnpid": "rs9923231", "variant_id": vid}])


def fetch(ids):
    pops = "populations{id ac an}"
    q = "{" + " ".join(f'v{i}: variant(variantId:"{v}", dataset: gnomad_r4){{ exome{{{pops}}} genome{{{pops}}} }}'
                       for i, v in enumerate(ids)) + "}"
    for attempt in range(8):
        r = requests.post(GNOMAD, json={"query": q}, timeout=120)
        if r.status_code == 200:
            return r.json().get("data") or {}
        time.sleep(20 * (attempt + 1))
    r.raise_for_status()


def pop_af(rec, pop):
    ac = an = 0
    for s in ("exome", "genome"):
        for p in (rec.get(s) or {}).get("populations", []):
            if p["id"] == pop:
                ac += p["ac"]; an += p["an"]
    return ac, an


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    al = single_snv_alleles()
    al = al[al.clinicalfunctionalstatus.isin(LOF)]
    keep = al.apply(lambda r: r["name"] in CORE.get(r["gene"], set())
                    or (r["gene"] == "DPYD" and r["name"] in DPYD_CORE), axis=1)
    al = al[keep]
    al = pd.concat([al, vkorc1_row()], ignore_index=True)
    print(f"{len(al)} single-SNV LoF alleles across {al.gene.nunique()} genes")
    ids, rows = list(al.variant_id.unique()), []
    for i in range(0, len(ids), BATCH):
        b = ids[i:i + BATCH]
        data = fetch(b)
        for k, v in enumerate(b):
            rec = data.get(f"v{k}")
            for pop in POPS:
                ac, an = pop_af(rec, pop) if rec else (0, 0)
                rows.append((v, pop, ac, an, rec is not None))
        time.sleep(5)
    f = pd.DataFrame(rows, columns=["variant_id", "pop", "ac", "an", "in_gnomad"])
    f["af"] = (f.ac / f.an).where(f.an > 0, 0.0)
    al = al.merge(f, on="variant_id")
    al.to_csv(OUT / "lof_allele_freqs.csv", index=False)
    uv = al.drop_duplicates(["gene", "variant_id", "pop"])  # alleles sharing one SNV count once
    q = uv.groupby(["gene", "pop"]).agg(q=("af", "sum"), n_variants=("variant_id", "nunique"),
                                        n_found=("in_gnomad", "sum")).reset_index()
    q["carrier"] = 1 - (1 - q.q.clip(upper=1)) ** 2
    q["two_copies"] = q.q ** 2
    q["drug"] = q.gene.map({**GENES, "VKORC1": "warfarin"})
    q.to_csv(OUT / "carrier_rates.csv", index=False)
    print(q.pivot(index="gene", columns="pop", values="carrier").round(3)[POPS])
    print("distinct variants per gene:", al.groupby("gene").variant_id.nunique().to_dict())
    print("not found in gnomAD:", al[~al.in_gnomad].drop_duplicates("variant_id")[["gene", "name"]].values.tolist())


if __name__ == "__main__":
    main()
