d = np.zeros((NParticles,3))
d[:,0] = np.linspace(0.85,2.5,NParticles)
potential = np.array(map(compute_lj_potential, d))
force = np.array(map(compute_lj_force, d))
