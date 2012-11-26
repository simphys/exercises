def compute_lj_force(rij, eps = EPS, sig = SIG, rcoff = RCOFF):
    norm = np.linalg.norm(rij)
    if norm > rcoff: return np.zeros(3)
    q = sig/norm
    return 4*eps*(12*q**11-6*q**5)*q/norm**2*rij
