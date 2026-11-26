# Partition truth by the variance route (the model of viprs_annot): one causal set of proportion PI over all SNPs,
# annotated or not, raised inside small categories to at least MIN_CAUSAL expected causal SNPs; a causal SNP's
# effect is the sum of independent components: the intercept beta_j^(0), t(RESID_DF) with unit variance, on every
# causal SNP, and beta_j^(c) ~ N(0, A_jc sigma2_c) for each category (A_jc = membership
# of category c, the largest of its annotation columns / 100) with the slab class (shares MIX, relative variances
# MIXVAR) common to the categories of the SNP. Each component is rescaled so that its genetic variance beta_c' R beta_c equals h2 * share_c,
# then the total so that beta' R beta = h2 exactly (the cross terms between components make the realised component
# h2 differ slightly from the design; they are printed and saved).
# b_hat = X (X' beta) + X z / sqrt(n), z ~ N(0, I_k), X = U Lambda^{1/2}; R is never formed.
# usage: python gen_var.py <ld_dir> <snp_table: tsv with ID, block, ld_index, A1, A2, freq> <out_dir> <seed>
#   env: NGWAS (300000), H2 (0.3), PI (0.01), BLOCKS ("1-100" or "" = all), MIX ("0.6,0.3,0.1"),
#        MIXVAR ("0.01,0.1,1"), RESID_DF (3), MIN_CAUSAL (100: expected causal SNPs per category at least; 0 = off),
#        ANNOT: tsv / parquet with ID and annotation columns (0-100 or 0-1); without it K (3) random 0/1 columns
#               t1..tK on PT (0.1) of the SNPs each, not overlapping,
#        DESIGN: categories and h2 shares, "name=col1+col2:share,..." (intercept: 1 - the sum), e.g.
#               "coding=High+Moderate:0.5,liver=Hepatocyte_Cholangiocyte:0.3"; default every random column at
#               (1 - RESID) / K, RESID (0.2).
# writes trait.ma, annot.tsv (ID, Block, Index and the annotation columns / 100), beta_comp.npy (M x (C + 1): the
# categories, then intercept), categories.tsv (name, columns, design h2, realised h2 = beta_c' R beta_c), snps.tsv.
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_io import read_X
ld, tab, out, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
env = os.environ.get
n = float(env("NGWAS", "300000")); h2 = float(env("H2", "0.3")); pi = float(env("PI", "0.01"))
mix = np.array([float(x) for x in env("MIX", "0.6,0.3,0.1").split(",")]); mix /= mix.sum()
mv = np.array([float(x) for x in env("MIXVAR", "0.01,0.1,1").split(",")]); kap = mv / (mix * mv).sum()
df = float(env("RESID_DF", "3")); min_causal = float(env("MIN_CAUSAL", "100"))
rng = np.random.default_rng(seed); os.makedirs(out, exist_ok=True)

snp = pd.read_csv(tab, sep="\t")
bl = env("BLOCKS", "")
blocks = sorted(snp.block.unique())
if bl:
    lo, hi = map(int, bl.split("-")); blocks = [b for b in blocks if lo <= b <= hi]
snp = snp[snp.block.isin(blocks)].sort_values(["block", "ld_index"]).reset_index(drop=True)
M = len(snp)

if env("ANNOT"):
    f = env("ANNOT")
    a = pd.read_parquet(f) if f.endswith(".parquet") else pd.read_csv(f, sep="\t")
    cols = [c for c in a.columns if c not in ("ID", "CHR", "BP37", "BP38", "SNP", "A1", "A2")]
    a = a.drop_duplicates("ID").set_index("ID")
    print(f"annotation rows matched by ID: {int(snp.ID.isin(a.index).sum())} of {M} SNPs")
    a = a.reindex(snp.ID)[cols].fillna(0).astype(float)
    if a.values.max() > 1: a = a / 100.
else:
    K = int(env("K", "3")); pt = float(env("PT", "0.1"))
    t = rng.choice(K + 1, size=M, p=[pt] * K + [1 - K * pt])
    cols = [f"t{k+1}" for k in range(K)]
    a = pd.DataFrame({c: (t == k).astype(float) for k, c in enumerate(cols)})
A = a.values
if env("DESIGN"):
    cats = []
    for item in env("DESIGN").split(","):
        nm, rest = item.split("="); cc, sh = rest.split(":"); cats.append((nm, cc.split("+"), float(sh)))
else:
    r = float(env("RESID", "0.2")); cats = [(c, [c], (1 - r) / len(cols)) for c in cols]
for nm, cc, _ in cats:
    miss = [c for c in cc if c not in cols]
    if miss: raise SystemExit(f"DESIGN names unknown columns: {miss}")
share_r = 1 - sum(s for _, _, s in cats)
if share_r < 0: raise SystemExit("DESIGN shares sum above 1")
memb = np.stack([A[:, [cols.index(c) for c in cc]].max(1) for _, cc, _ in cats], 1)   # M x C membership

causal = rng.random(M) < pi
# coverage-stratified causal set: each category gets at least MIN_CAUSAL expected causal SNPs among its members
# (extra causal SNPs drawn there with p = (MIN_CAUSAL - PI m_c) / (m_c (1 - PI)), m_c = members), so that sigma2_c is
# estimable; 0 turns it off
for i in range(memb.shape[1]):
    mem = memb[:, i] > 0; mc = int(mem.sum())
    if mc == 0: continue
    p_add = min(1., max(0., (min_causal - pi * mc) / (mc * (1 - pi))))
    if p_add > 0: causal |= mem & (rng.random(M) < p_add)
cls = rng.choice(len(mix), size=M, p=mix)
sc = np.where(causal, np.sqrt(kap[cls] / pi), 0.)
C = np.zeros((M, len(cats) + 1))
for i in range(len(cats)): C[:, i] = sc * np.sqrt(memb[:, i]) * rng.standard_normal(M)
C[:, -1] = np.where(causal, rng.standard_t(df, M) / np.sqrt(df / (df - 2)), 0.)

def gram(V):
    # V' R V over the blocks (one pass, X not kept), plus the X' V per block for later
    Q = np.zeros((V.shape[1], V.shape[1])); T = {}
    for b in blocks:
        X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))
        idx = np.flatnonzero(snp.block.values == b); li = snp.ld_index.values[idx]
        T[b] = X[li].T @ V[idx]; Q += T[b].T @ T[b]
    return Q, T
Q, _ = gram(C)
target = np.array([s for _, _, s in cats] + [share_r]) * h2
for i in range(C.shape[1]):
    if Q[i, i] > 0: C[:, i] *= np.sqrt(target[i] / Q[i, i])
    elif target[i] > 0: raise SystemExit(f"component {i} has no causal SNP: raise PI or the share")
beta = C.sum(1)
Q, T = gram(np.c_[C, beta])
f = np.sqrt(h2 / Q[-1, -1]); C *= f; beta *= f
own = np.diag(Q)[:-1] * f * f                                         # beta_c' R beta_c: h2 of each component
bh = np.empty(M)
for b in blocks:
    X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))
    idx = np.flatnonzero(snp.block.values == b); li = snp.ld_index.values[idx]
    z = rng.standard_normal(X.shape[1])
    bh[idx] = X[li] @ (T[b][:, -1] * f) + (X @ z)[li] / np.sqrt(n)
names = [nm for nm, _, _ in cats] + ["intercept"]
a1, a2 = ("A1", "A2") if "A1" in snp else ("alt", "ref")
fr = snp["freq"] if "freq" in snp else snp.get("AF_PanUKB", snp.get("A1Freq"))
pd.DataFrame(dict(SNP=snp.ID, A1=snp[a1], A2=snp[a2], freq=fr, b=bh, se=1 / np.sqrt(n), p=0.5, N=n)).to_csv(f"{out}/trait.ma", sep="\t", index=False)
ann = a.reset_index(drop=True).copy(); ann.insert(0, "Index", snp.ld_index.values); ann.insert(0, "Block", snp.block.values)
ann.insert(0, "ID", snp.ID.values); ann.to_csv(f"{out}/annot.tsv", sep="\t", index=False)
np.save(f"{out}/beta_comp.npy", C); snp[["ID", "block", "ld_index"]].to_csv(f"{out}/snps.tsv", sep="\t", index=False)
pd.DataFrame(dict(name=names, columns=["+".join(cc) for _, cc, _ in cats] + [""], design=target,
                  h2=own)).to_csv(f"{out}/categories.tsv", sep="\t", index=False)
print(M, "SNPs", len(blocks), "blocks;", int(causal.sum()), "causal; unannotated SNPs", int((A.max(1) == 0).sum()))
print("causal per category", dict(zip(names[:-1], [int((causal & (memb[:, i] > 0)).sum()) for i in range(memb.shape[1])])))
print("h2 design  ", dict(zip(names, np.round(target, 4))))
print("h2 realised", dict(zip(names, np.round(own, 4))), "(beta_c' R beta_c; total", h2, "exactly)")
