# Partition truth simulated from the variance model of viprs_annot itself: E[beta_j^2] = sigma2_0 + sum_c A_jc sigma2_c,
# sigma2_c = ratio_c sigma2_0 (the ratios of DESIGN, e.g. on the MASH atoms 1, 10, 100), the intercept on every SNP.
# One causal set of proportion PI over all SNPs; a causal SNP's effect is the sum of independent components, the
# intercept ~ N(0, kappa sigma2_0 / PI) and per category c ~ N(0, A_jc kappa sigma2_c / PI) (A_jc = membership of c, the
# largest of its annotation columns / 100), with the slab class kappa (shares MIX, relative variances MIXVAR, E kappa = 1)
# common to the components of the SNP; sigma2_0 is set so that beta' R beta = h2 exactly. Realised component h2
# (beta_c' R beta_c) is printed and saved.
# b_hat = X (X' beta) + X z / sqrt(n), z ~ N(0, I_k), X = U Lambda^{1/2}; R is never formed.
# usage: python gen_var.py <ld_dir> <snp_table: tsv with ID, block, ld_index, A1, A2, freq> <out_dir> <seed>
#   env: NGWAS (300000), H2 (0.3), PI (0.01), BLOCKS ("1-100" or "" = all), MIX ("0.6,0.3,0.1"),
#        MIXVAR ("0.01,0.1,1"),
#        ANNOT: tsv / parquet with ID and annotation columns (0-100 or 0-1); without it random 0/1 columns (K, PT,
#               SEG, DECOY, NOISE: see below),
#        DESIGN: categories and variance ratios sigma2_c / sigma2_0, "name=col1+col2:ratio,...", e.g.
#               "coding=High+Moderate:100,liver=Hepatocyte_Cholangiocyte:10"; default every column (random annotations:
#               the targets t1..tK) at ratio 10.
# writes to out_dir trait.ma, beta_comp.npy (float32, M x (C + 1): the categories, then intercept), v_design.npy (float32,
# the design variance E[beta_j^2] of each SNP, for the model-definition truth of eval_var.py), oracle.tsv (annotation,
# sigma2, pi: the oracle of run_var.R method ORACLE, the variance model h2_j = sigma2_0 + sum_k A_jk sigma2_k on the
# design columns fitted to v_design by the likelihood stationarity sum_j A_jk (v_j - h2_j) / h2_j^2 = 0 (MM steps to
# convergence, every SNP); pi the realised share of causal SNPs), categories.tsv
# (name, columns, ratio, design h2 = sigma2_c sum_j A_jc, realised h2 = beta_c' R beta_c); annot.tsv (ID, Block, Index and the annotation columns
# / 100) and snps.tsv to SHARED (env, default out_dir; random annotations: always out_dir) unless already there.
import os, sys, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from sim_io import read_X
ld, tab, out, seed = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4])
env = os.environ.get
n = float(env("NGWAS", "300000")); h2 = float(env("H2", "0.3")); pi = float(env("PI", "0.01"))
mix = np.array([float(x) for x in env("MIX", "0.6,0.3,0.1").split(",")]); mix /= mix.sum()
mv = np.array([float(x) for x in env("MIXVAR", "0.01,0.1,1").split(",")]); kap = mv / (mix * mv).sum()
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
    # random annotations in runs of SEG neighbouring SNPs (annotations follow genomic regions, so neighbours in LD
    # share them): K targets t1..tK of coverage PT each, not overlapping; DECOY = 1 adds per target a large column
    # big<k> (the target plus 25% of all runs, like an enhancer set) and a partial one part<k> (half of the target
    # runs plus PT / 2 others); NOISE unrelated columns of coverage PT
    K = int(env("K", "3")); pt = float(env("PT", "0.1")); seg = int(env("SEG", "10"))
    run = np.arange(M) // seg; ns = int(run.max()) + 1
    ts = rng.choice(K + 1, size=ns, p=[pt] * K + [1 - K * pt])
    cols = [f"t{k+1}" for k in range(K)]
    a = {c: (ts == k)[run].astype(float) for k, c in enumerate(cols)}
    if int(env("DECOY", "0")):
        for k in range(K):
            a[f"big{k+1}"] = ((ts == k) | (rng.random(ns) < 0.25))[run].astype(float)
            a[f"part{k+1}"] = (((ts == k) & (rng.random(ns) < 0.5)) | (rng.random(ns) < pt / 2))[run].astype(float)
    for i in range(int(env("NOISE", "0"))): a[f"n{i+1}"] = (rng.random(ns) < pt)[run].astype(float)
    cols = list(a); a = pd.DataFrame(a)
A = a.values
if env("DESIGN"):
    cats = []
    for item in env("DESIGN").split(","):
        nm, rest = item.split("="); cc, rt = rest.split(":"); cats.append((nm, cc.split("+"), float(rt)))
else:
    tc = [c for c in cols if c[0] == "t"] if not env("ANNOT") else cols   # random annotations: the targets only
    cats = [(c, [c], 10.) for c in tc]
for nm, cc, _ in cats:
    miss = [c for c in cc if c not in cols]
    if miss: raise SystemExit(f"DESIGN names unknown columns: {miss}")
memb = np.stack([A[:, [cols.index(c) for c in cc]].max(1) for _, cc, _ in cats], 1)   # M x C membership
ratio = np.array([r for _, _, r in cats])

# the variance model itself: causal with probability PI everywhere, slab class common to the components of a SNP,
# component c ~ N(0, memb_c ratio_c kappa / PI), intercept ~ N(0, kappa / PI) on every causal SNP; so
# E[beta_j^2] = sigma2_0 (1 + sum_c memb_jc ratio_c), sigma2_c = ratio_c sigma2_0 (sigma2_0 set by h2 below)
causal = rng.random(M) < pi
cls = rng.choice(len(mix), size=M, p=mix)
sc = np.where(causal, np.sqrt(kap[cls] / pi), 0.)
C = np.zeros((M, len(cats) + 1))
for i in range(len(cats)): C[:, i] = sc * np.sqrt(memb[:, i] * ratio[i]) * rng.standard_normal(M)
C[:, -1] = sc * rng.standard_normal(M)
if not causal.any(): raise SystemExit("no causal SNP: raise PI")

def gram(V):
    # V' R V over the blocks (one pass, X not kept), plus the X' V per block for later
    Q = np.zeros((V.shape[1], V.shape[1])); T = {}
    for b in blocks:
        X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))
        idx = np.flatnonzero(snp.block.values == b); li = snp.ld_index.values[idx]
        T[b] = X[li].T @ V[idx]; Q += T[b].T @ T[b]
    return Q, T
beta = C.sum(1)
Q, T = gram(np.c_[C, beta])
f = np.sqrt(h2 / Q[-1, -1]); C *= f; beta *= f
s20 = f * f                                                           # sigma2_0
vdes = s20 * (1. + memb @ ratio)                                      # design variance E[beta_j^2]
target = np.r_[ratio * s20 * memb.sum(0), s20 * M]                    # design h2 sigma2_c sum_j memb_jc (no LD)
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
# annot.tsv and snps.tsv do not change between seeds with ANNOT: written once to SHARED (default out_dir) if missing
sh = env("SHARED", out) if env("ANNOT") else out
os.makedirs(sh, exist_ok=True)
if not env("ANNOT") or not os.path.exists(f"{sh}/annot.tsv"):
    ann = a.reset_index(drop=True).copy(); ann.insert(0, "Index", snp.ld_index.values); ann.insert(0, "Block", snp.block.values)
    ann.insert(0, "ID", snp.ID.values); tmp = f"{sh}/annot.tsv.{os.getpid()}"   # written aside, then renamed: parallel tasks never read a partial file
    ann.astype({c: "float32" for c in cols}).to_csv(tmp, sep="\t", index=False); os.replace(tmp, f"{sh}/annot.tsv")
if not env("ANNOT") or not os.path.exists(f"{sh}/snps.tsv"):
    tmp = f"{sh}/snps.tsv.{os.getpid()}"
    snp[["ID", "block", "ld_index"]].to_csv(tmp, sep="\t", index=False); os.replace(tmp, f"{sh}/snps.tsv")
np.save(f"{out}/beta_comp.npy", C.astype(np.float32)); np.save(f"{out}/v_design.npy", vdes.astype(np.float32))
# oracle variance model on the design columns (the same MM fit as the model-definition truth of eval_var.py)
oc = list(dict.fromkeys(sum([cc for _, cc, _ in cats], [])))
Ao = np.c_[np.ones(M), A[:, [cols.index(c) for c in oc]]]
s2 = np.full(Ao.shape[1], vdes.sum() / Ao.sum())
for it in range(100000):
    hh = Ao @ s2; N_ = Ao.T @ (vdes / hh ** 2); D_ = Ao.T @ (1 / hh)
    r_ = np.divide(N_, D_, out=np.ones_like(N_), where=D_ > 0); s2 = s2 * np.sqrt(r_)
    if np.max(np.abs(r_[s2 > 0] - 1)) < 1e-8: break
pd.DataFrame(dict(annotation=["intercept"] + oc, sigma2=s2, pi=float(causal.mean()))).to_csv(f"{out}/oracle.tsv", sep="\t", index=False)
pd.DataFrame(dict(name=names, columns=["+".join(cc) for _, cc, _ in cats] + [""], ratio=np.r_[ratio, 1.], design=target,
                  h2=own)).to_csv(f"{out}/categories.tsv", sep="\t", index=False)
print(M, "SNPs", len(blocks), "blocks;", int(causal.sum()), "causal; unannotated SNPs", int((A.max(1) == 0).sum()))
print("causal per category", dict(zip(names[:-1], [int((causal & (memb[:, i] > 0)).sum()) for i in range(memb.shape[1])])))
print("h2 design  ", dict(zip(names, np.round(target, 4))))
print("h2 realised", dict(zip(names, np.round(own, 4))), "(beta_c' R beta_c; total", h2, "exactly)")
