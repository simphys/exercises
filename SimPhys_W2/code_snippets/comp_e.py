double compute_energy(double *x, double *v) {
	double E_pot = 0;
	double E_kin = 0;
	double rij[3];

    // sum up potential energies
    for (int i=1; i < N; i++) {
    	for (int j = 0; j < i; j++) {
			// distance vector
    		minimum_image(x, i, j, rij);
            E_pot += compute_lj_potential(rij);
    	}
    }

    // sum up kinetic energy
    for (int i=0; i < N; i++) {
    	double dot = 0;
    	dot +=  v[i]*v[i];
    	dot +=  v[i+N]*v[i+N];
    	dot +=  v[i+2*N]*v[i+2*N];
        E_kin += 0.5 * dot;
    }
    return E_pot + E_kin;
}
