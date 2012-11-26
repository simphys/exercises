def compute_forces(x, L = L):
    """Compute and return the forces acting onto the particles,
    depending on the positions x."""
    global epsilon, sigma
    _, N = x.shape
    f = np.zeros_like(x)
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = x[:,j] - x[:,i]
            rij -= np.rint(rij/L)*L #BECAUSE OF PBC
            fij = compute_lj_force(rij)
            f[:,i] -= fij
            f[:,j] += fij
    return f
