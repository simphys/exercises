def compute_lj_potential(rij, eps = EPS, sig = SIG, rcoff = RCOFF,\
                         vcoff = VCOFF):
    norm = np.linalg.norm(rij)
    if norm > rcoff: return 0
    q = sig/norm
    return 4*eps*(q**12-q**6) - vcoff
