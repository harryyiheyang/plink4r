# Per category (categories.tsv of gen_var.py): prediction R2 of the category's partitioned PRS (the sum of the
# BETA_STD_<col> of its columns) for the true category genetic value, and the Shapley heritability ratio
# _h2 = estimated (h2_shapley of <method>_h2.tsv, summed over the columns of the category: Shapley values add up)
# / true Shapley value beta_c' R beta (beta the true total effect); true values printed as "true shapley"; _or: the
# oracle ratio, the estimated per-SNP share BETA_c / BETA applied to the true beta (w_c beta)' R beta / truth.
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
def r2(i, j): return Q[i, j] ** 2 / (Q[i, i] * Q[j, j]) if Q[i, i] * Q[j, j] > 0 else float("nan")
tot = L[("true", "total")]
tsh = {nm: Q[L[("true", nm)], tot] for nm in names}
print("true h2", dict(zip(names, np.round(cat.h2, 4))))
print("true shapley", {k: round(float(v), 4) for k, v in tsh.items()})
for m in sys.argv[3:]:
    out = {"total": round(float(r2(L[(m, "total")], L[("true", "total")])), 3)}
    h = None
    if os.path.exists(f"{sim}/{m}_h2.tsv"):
        h = pd.read_csv(f"{sim}/{m}_h2.tsv", sep="\t").set_index("annotation").h2_shapley
    for nm, cc in zip(names, cols):
        if (m, nm) not in L: continue
        i, j = L[(m, nm)], L[("true", nm)]
        out[nm] = round(float(r2(i, j)), 3)
        if h is not None: out[nm + "_h2"] = round(float(h.reindex(cc).fillna(0).sum() / tsh[nm]), 3)
        # oracle: the estimated allocation share w = BETA_k / BETA of each SNP applied to the true beta
        bt_ = V[:, L[(m, "total")]]; w = np.divide(V[:, i], bt_, out=np.zeros_like(bt_), where=np.abs(bt_) > 1e-12)
        out[nm + "_or"] = round(float((w * V[:, tot]) @ Rb / tsh[nm]), 3)
    if h is not None:
        used = set(sum(cols, [])); out["other_h2"] = round(float(h[[k for k in h.index if k not in used]].sum()), 4)
    print(m, out)
