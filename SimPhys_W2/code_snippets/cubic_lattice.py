# total number of particles
N = n*n*n

# particle positions on a cubic lattice
x = np.zeros((3,N))

# compute system size
L = n/density**(1/3)

# - set up n*n*n particles on a cubic lattice
positions = np.arange(0.5, n + 0.5)/density**(1/3)
for a in range(n):
    for b in range(n):
        for c in range(n):
            x[0,(a*n*n)+(b*n)+c] = positions[a]
            x[1,(a*n*n)+(b*n)+c] = positions[b]
            x[2,(a*n*n)+(b*n)+c] = positions[c]
