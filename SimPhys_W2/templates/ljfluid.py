from numpy import *
from matplotlib.pyplot import *

def minimum_image(ri, rj):
    global L
    # compute distance
    rij = rj-ri
    # wrap the distance into [-0.5, 0.5]
    rij -= rint(rij/L)*L
    return rij

def compute_lj_force(rij):
    # TODO
    # .
    # .
    # .
    return fij

def compute_lj_potential(rij):
    # TODO
    # .
    # .
    # .
    return pot

def compute_forces(x, L):
    """Compute and return the forces acting onto the particles,
    depending on the positions x."""
    _, N = x.shape
    f = zeros_like(x)
    rijs = empty_like(x)
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = minimum_image(x[:,i], x[:,j])
            fij = compute_lj_force(rij)
            f[:,i] -= fij
            f[:,j] += fij
    return f

def compute_energy(x, v, L):
    """Compute and return the total energy of the system with the
    particles at positions x."""
    _, N = x.shape
    E_pot = 0.0
    E_kin = 0.0
    # sum up potential energies
    for i in range(1,N):
        for j in range(i):
            # distance vector
            rij = minimum_image(x[:,i], x[:,j])
            E_pot += compute_lj_potential(rij)
    # sum up kinetic energy
    for i in range(N):
        E_kin += 0.5 * dot(v[:,i],v[:,i])
    return E_pot + E_kin
    
def step_vv(x, v, f, dt):
    # update positions
    x += v*dt + 0.5*f * dt*dt
    # half update of the velocity
    v += 0.5*f * dt
        
    # compute new forces
    f = compute_forces(x, L)
    # we assume that m=1 for all particles

    # second half update of the velocity
    v += 0.5*f * dt

    return x, v, f

# CONSTANTS
# density
density = 0.7
# number of particles per side
n = 3
# timestep
dt = 0.01
# length of run
tmax = 1.0
# cutoff
rcut = 2.5
# potential shift
shift = -0.016316891136
# total number of particles
N = n*n*n

# particle positions on a cubic lattice
x = empty((3,N))

# TODO: 
# - compute the size of the system
# - set up n*n*n particles on a cubic lattice
# .
# .
# .


# random particle velocities
v = 2.0*random.random((3,N))-1.0

# RUNNING VARIABLES
t = 0.0

# variables to cumulate data
ts = []
Es = []
traj = []

# open the trajectory file
vtffile = open('ljfluid.vtf', 'w')
# write the structure of the system into the file: 
# N particles ("atoms") with a radius of 0.5
vtffile.write('atom 0:%s radius 0.5\n' % (N-1))
vtffile.write('pbc %s %s %s\n' % (L, L, L))

# main loop
f = compute_forces(x, L)
while t < tmax:
    x, v, f = step_vv(x, v, f, dt)
    t += dt

    E = compute_energy(x, v, L)
    print "t=%s, E=%s" % (t, E)

    ts.append(t)
    Es.append(E)
    
    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))

vtffile.close()
plot(ts, Es)
        
