# SNP table for blocks LO..HI from blockN_variants.tsv (row order = eigen row): usage python make_ref.py <ld_dir> <LO-HI> <out.tsv>
import sys, os, pandas as pd
ld, rng_, out = sys.argv[1], sys.argv[2], sys.argv[3]; lo, hi = map(int, rng_.split("-"))
parts = []
for b in range(lo, hi + 1):
    f = os.path.join(ld, f"block{b}_variants.tsv")
    if not os.path.exists(f): continue
    v = pd.read_csv(f, sep="\t"); v["block"] = b; v["ld_index"] = range(len(v)); parts.append(v)
t = pd.concat(parts, ignore_index=True)
t = t.rename(columns={"Chrom": "chr", "PhysPos": "bp37", "A1Freq": "AF_PanUKB"})
t["alt"], t["ref"], t["freq"] = t["A1"], t["A2"], t["AF_PanUKB"]
t[["ID", "chr", "bp37", "ref", "alt", "AF_PanUKB", "block", "ld_index", "A1", "A2", "freq"]].to_csv(out, sep="\t", index=False)
print(len(t), "SNPs in", len(parts), "blocks")
