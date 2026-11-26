import os, numpy as np
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

