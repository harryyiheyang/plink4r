# Per category (categories.tsv of gen_var.py): prediction R2 of the category's partitioned PRS (the sum of the
# BETA_STD_<col> of its columns) for the true category genetic value, and the Shapley heritability ratio
# _h2 = estimated (h2_shapley of <method>_h2.tsv, summed over the columns of the category: Shapley values add up)
# / true Shapley value beta_c' R beta (beta the true total effect); true values printed as "true shapley"; _vs: the same
# with var_shapley (b_c' R b of the output posterior means instead of the mean over the draws); _or: the
# oracle ratio, the estimated per-SNP share BETA_c / BETA applied to the true beta (w_c beta)' R beta / truth.
# Model-definition truth (per method, for the columns the method was given, read from <method>_h2.tsv and annot.tsv of
# SHARED): sigma2_k of the variance model h2_j = sigma2_0 + sum_k A_jk sigma2_k fitted to the design variance
# v_j = E[beta_j^2] of gen_var.py (v_design.npy) by the likelihood stationarity sum_j A_jk (v_j - h2_j) / h2_j^2 = 0
# (MM steps to convergence; the SNPs not suppressed, suppressed.npy of gen_var.py); the true Shapley value of column k is
# then (w_k beta)' R beta with w_jk = A_jk sigma2_k / h2_j and the true beta (k = 0 the intercept on every SNP).
# _mh2 = estimated h2_shapley / this truth; printed as "<method> model truth" (for ORACLE: its own fixed sigma2_k).
# usage: python eval_var.py <ld_dir> <sim_dir> <method>...   env SHARED: directory of snps.tsv (default sim_dir)
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from sim_io import read_X
ld, sim = sys.argv[1], sys.argv[2]
snp = pd.read_csv(os.path.join(os.environ.get("SHARED", sim), "snps.tsv"), sep="\t")
C = np.load(f"{sim}/beta_comp.npy").astype(float)
cat = pd.read_csv(f"{sim}/categories.tsv", sep="\t", keep_default_na=False)
names = list(cat.name); cols = [c.split("+") if c else ["intercept"] for c in cat["columns"]]
vecs, labels = [C[:, i] for i in range(C.shape[1])] + [C.sum(1)], [("true", nm) for nm in names] + [("true", "total")]
for m in sys.argv[3:]:
    d = pd.read_csv(f"{sim}/{m}.tsv", sep="\t").set_index("ID").reindex(snp.ID).fillna(0)
    vecs.append(d.BETA_STD.values); labels.append((m, "total"))
    for nm, cc in zip(names, cols):
        k = [f"BETA_STD_{c}" for c in cc if f"BETA_STD_{c}" in d]
        if k: vecs.append(d[k].sum(axis=1).values); labels.append((m, nm))
V = np.stack(vecs, 1); Q = np.zeros((V.shape[1], V.shape[1])); Rb = np.zeros(len(snp))
for b, g in snp.groupby("block", sort=True):
    X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))[g.ld_index.values]
    T = X.T @ V[g.index.values]; Q += T.T @ T; Rb[g.index.values] = X @ T[:, len(names)]
L = {lab: i for i, lab in enumerate(labels)}
ann = None
def model_truth(acols, beta, Rb, tol=1e-8, maxit=100000):
    global ann
    if ann is None: ann = pd.read_csv(os.path.join(os.environ.get("SHARED", sim), "annot.tsv"), sep="\t").set_index("ID").reindex(snp.ID).fillna(0)
    keep = ~np.load(f"{sim}/suppressed.npy") if os.path.exists(f"{sim}/suppressed.npy") else np.ones(len(snp), bool)
    A = np.c_[np.ones(len(snp)), ann[acols].values.astype(float)][keep]
    v = np.load(f"{sim}/v_design.npy").astype(float)[keep]
    s2 = np.full(A.shape[1], v.sum() / A.sum())
    for it in range(maxit):
        h = A @ s2; N = A.T @ (v / h ** 2); D = A.T @ (1 / h)
        r = np.divide(N, D, out=np.ones_like(N), where=D > 0); s2 = s2 * np.sqrt(r)
        if np.max(np.abs(r[s2 > 0] - 1)) < tol: break
    h = A @ s2; W = A * s2 / h[:, None]
    return dict(zip(["intercept"] + list(acols), W.T @ (beta * Rb)[keep])), s2
def r2(i, j): return Q[i, j] ** 2 / (Q[i, i] * Q[j, j]) if Q[i, i] * Q[j, j] > 0 else float("nan")
tot = L[("true", "total")]
tsh = {nm: Q[L[("true", nm)], tot] for nm in names}
print("true h2", dict(zip(names, np.round(cat.h2, 4))))
print("true shapley", {k: round(float(v), 4) for k, v in tsh.items()})
for m in sys.argv[3:]:
    out = {"total": round(float(r2(L[(m, "total")], L[("true", "total")])), 3)}
    h = None
    if os.path.exists(f"{sim}/{m}_h2.tsv"):
        ht = pd.read_csv(f"{sim}/{m}_h2.tsv", sep="\t").set_index("annotation"); h, hv = ht.h2_shapley, ht.var_shapley
        ac = [k for k in ht.index if k != "intercept"]; mt, _ = model_truth(ac, V[:, tot], Rb)
        print(m, "model truth", {nm: round(float(sum(mt.get(c, 0.) for c in cc)), 4) for nm, cc in zip(names, cols)},
              "other", round(float(sum(v for k, v in mt.items() if k not in set(sum(cols, [])))), 4))
    for nm, cc in zip(names, cols):
        if (m, nm) not in L: continue
        i, j = L[(m, nm)], L[("true", nm)]
        out[nm] = round(float(r2(i, j)), 3)
        if h is not None:
            out[nm + "_h2"] = round(float(h.reindex(cc).fillna(0).sum() / tsh[nm]), 3)
            out[nm + "_vs"] = round(float(hv.reindex(cc).fillna(0).sum() / tsh[nm]), 3)
            tm = sum(mt.get(c, 0.) for c in cc)
            out[nm + "_mh2"] = round(float(h.reindex(cc).fillna(0).sum() / tm), 3) if tm > 0 else float("nan")
        # oracle: the estimated allocation share w = BETA_k / BETA of each SNP applied to the true beta
        bt_ = V[:, L[(m, "total")]]; w = np.divide(V[:, i], bt_, out=np.zeros_like(bt_), where=np.abs(bt_) > 1e-12)
        out[nm + "_or"] = round(float((w * V[:, tot]) @ Rb / tsh[nm]), 3)
    if h is not None:
        used = set(sum(cols, [])); out["other_h2"] = round(float(h[[k for k in h.index if k not in used]].sum()), 4)
    print(m, out)
