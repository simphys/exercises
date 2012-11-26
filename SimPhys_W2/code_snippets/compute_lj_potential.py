def compute_lj_potential(rij, eps = EPS, sig = SIG):
    q = sig/np.linalg.norm(rij)
    return 4*eps*(q**12-q**6)
