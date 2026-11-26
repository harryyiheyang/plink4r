# per-component R2 (LD-weighted r^2) of estimated weights vs truth; usage: python eval_ideal.py <ld_dir> <sim_dir> <est.tsv>...
import sys, os, numpy as np, pandas as pd
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__))); from gen_ideal_io import read_X
ld, sim = sys.argv[1], sys.argv[2]
snp = pd.read_csv(f"{sim}/snps.tsv", sep="\t"); C = np.load(f"{sim}/beta_comp.npy"); K = C.shape[1] - 1
names = [f"t{k+1}" for k in range(K)] + ["residual"]
# all vectors to compare: truth columns, then per file total + components; one pass over blocks (X not kept)
vecs, labels = [C[:, k] for k in range(K + 1)] + [C.sum(1)], [("true", nm) for nm in names] + [("true", "total")]
for f in sys.argv[3:]:
    d = pd.read_csv(f, sep="\t").set_index("ID").reindex(snp.ID).fillna(0)
    vecs.append(d.BETA_STD.values); labels.append((f, "total"))
    for nm in names:
        if f"BETA_STD_{nm}" in d: vecs.append(d[f"BETA_STD_{nm}"].values); labels.append((f, nm))
V = np.stack(vecs, 1); Q = np.zeros((V.shape[1], V.shape[1]))
for b, g in snp.groupby("block", sort=True):
    X = read_X(os.path.join(ld, f"block{b}.eigen.bin"))[g.ld_index.values]
    T = X.T @ V[g.index.values]; Q += T.T @ T                       # Gram matrix under R
L = {lab: i for i, lab in enumerate(labels)}
tv = np.array([Q[L[("true", nm)], L[("true", nm)]] for nm in names])
print("true share", dict(zip(names, np.round(tv / tv.sum(), 3))))
def r2(i, j): return Q[i, j] ** 2 / (Q[i, i] * Q[j, j]) if Q[i, i] * Q[j, j] > 0 else float("nan")
for f in sys.argv[3:]:
    out = {"total": r2(L[(f, "total")], L[("true", "total")])}
    for nm in names:
        if (f, nm) in L:
            i, j = L[(f, nm)], L[("true", nm)]; out[nm] = r2(i, j); out[nm + "_var"] = Q[i, i] / Q[j, j]
    print(os.path.basename(f), {k: round(float(v), 3) for k, v in out.items()})
