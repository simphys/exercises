# main loop
while t < tmax:
    xfolded = x.copy()
    xfolded -= np.floor(x/L)*L  #BECAUSE OF PBC
    
    x, v, f = step_vv(xfolded, v, f, dt)
    t += dt
    
    traj.append(x.copy())
    Es.append(compute_energy(xfolded, v))
    
    # write out that a new timestep starts
    vtffile.write('timestep\n')
    # write out the coordinates of the particles
    for i in range(N):
        vtffile.write("%s %s %s\n" % (x[0,i], x[1,i], x[2,i]))

vtffile.close()

traj = np.array(traj)
