# Ideal-case truth in eigen space: K non-overlapping tissue annotations, causal proportion PCAUSAL of all SNPs, each
# causal SNP's effect comes from its tissue (or the residual if unannotated). Causal SNPs: PCAUSAL * M, drawn with
# weight ENRICH for tissue SNPs and 1 for the others; beta from the slab mixture of viprs_annot (shares MIX, relative
# variances MIXVAR), scaled to h2.
# b_hat = X (X' beta) + X z / sqrt(n), z ~ N(0, I_k), X = U Lambda^{1/2}; R is never formed.
# usage: python gen_ideal.py <ld_dir> <snp_table: tsv/parquet with ID, block, ld_index, A1, A2, freq> <out_dir> <seed>
#        env: NGWAS (300000), H2 (0.3), PCAUSAL (0.01), K (3), PT (0.1 per tissue), BLOCKS ("1-100" or "" = all),
#             ENRICH (1), MIX ("0.6,0.3,0.1"), MIXVAR ("0.01,0.1,1")
import os, sys, numpy as np, pandas as pd
ld, tab, out, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
n = float(os.environ.get("NGWAS", "300000")); h2 = float(os.environ.get("H2", "0.3"))
pc = float(os.environ.get("PCAUSAL", "0.01")); K = int(os.environ.get("K", "3")); pt = float(os.environ.get("PT", "0.1"))
en = float(os.environ.get("ENRICH", "1"))
mix = np.array([float(x) for x in os.environ.get("MIX", "0.6,0.3,0.1").split(",")]); mix /= mix.sum()
mv = np.array([float(x) for x in os.environ.get("MIXVAR", "0.01,0.1,1").split(",")])
rng = np.random.default_rng(seed); os.makedirs(out, exist_ok=True)

def read_X(path):
    sz = os.path.getsize(path)
    with open(path, "rb") as f:
        m, k = np.frombuffer(f.read(8), "<i4")
        if sz == 20 + 12 * k + 2 * m * k:      # plink4r int16 layout
            f.read(12); lam = np.frombuffer(f.read(8 * k), "<f8"); sc = np.frombuffer(f.read(4 * k), "<f4")
            q = np.frombuffer(f.read(2 * m * k), "<i2").reshape((m, k), order="F")
            return q * (sc * np.sqrt(np.maximum(lam, 0)))
        if sz == 16 + 4 * k + 4 * m * k:       # SBayesRC float layout
            f.read(8); lam = np.frombuffer(f.read(4 * k), "<f4").astype(float)
            U = np.frombuffer(f.read(4 * m * k), "<f4").reshape((m, k), order="F")
            return U * np.sqrt(np.maximum(lam, 0))
        raise ValueError("unknown eigen layout: " + path)

snp = pd.read_parquet(tab) if tab.endswith(".parquet") else pd.read_csv(tab, sep="\t")
bl = os.environ.get("BLOCKS", "")
blocks = sorted(snp.block.unique())
if bl:
    lo, hi = map(int, bl.split("-")); blocks = [b for b in blocks if lo <= b <= hi]
snp = snp[snp.block.isin(blocks)].sort_values(["block", "ld_index"]).reset_index(drop=True)
M = len(snp)
tissue = rng.choice(K + 1, size=M, p=[pt] * K + [1 - K * pt])          # K = unannotated
w = np.where(tissue < K, en, 1.)
causal = np.zeros(M, bool); causal[rng.choice(M, size=max(1, int(round(pc * M))), replace=False, p=w / w.sum())] = True
cls = rng.choice(len(mix), size=M, p=mix)
beta = np.where(causal, rng.standard_normal(M) * np.sqrt(mv[cls]), 0.)
# one pass per block (X is not kept): t_b = X_b' beta_b and the noise draw; scale to h2 afterwards
T, Zs, gv = {}, {}, 0.
for b in blocks:
    X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))
    idx = np.flatnonzero(snp.block.values == b); li = snp.ld_index.values[idx]
    T[b] = X[li].T @ beta[idx]; gv += T[b] @ T[b]
    z = rng.standard_normal(X.shape[1]); Zs[b] = (X @ z)[li]
sc = np.sqrt(h2 / gv); beta *= sc
bh = np.empty(M)
for b in blocks:
    X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))
    idx = np.flatnonzero(snp.block.values == b); li = snp.ld_index.values[idx]
    bh[idx] = X[li] @ (T[b] * sc) + Zs[b] / np.sqrt(n)
C = np.zeros((M, K + 1)); C[np.arange(M), tissue] = beta               # columns t1..tK, residual
a1, a2 = ("A1", "A2") if "A1" in snp else ("alt", "ref")
fr = snp["freq"] if "freq" in snp else snp.get("AF_PanUKB", snp.get("A1Freq"))
pd.DataFrame(dict(SNP=snp.ID, A1=snp[a1], A2=snp[a2], freq=fr, b=bh, se=1 / np.sqrt(n), p=0.5, N=n)).to_csv(f"{out}/trait.ma", sep="\t", index=False)
ann = pd.DataFrame({f"t{k+1}": (tissue == k).astype(float) for k in range(K)}); ann.insert(0, "Index", snp.ld_index.values); ann.insert(0, "Block", snp.block.values); ann.insert(0, "ID", snp.ID.values)
ann.to_csv(f"{out}/annot.tsv", sep="\t", index=False)
np.save(f"{out}/beta_comp.npy", C); snp[["ID", "block", "ld_index"]].to_csv(f"{out}/snps.tsv", sep="\t", index=False)
print(M, "SNPs", len(blocks), "blocks; causal per component", (C != 0).sum(0))
