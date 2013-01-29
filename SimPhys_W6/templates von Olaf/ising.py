from numpy import *
from matplotlib.pyplot import *
import sys

def compute_energy(sigma):
    n, m = sigma.shape
    shifted_ixs = range(1,n) + [0]
    E = -(sigma*sigma[shifted_ixs,:]).sum()
    shifted_ixs = range(1,m) + [0]
    E += -(sigma*sigma[:,shifted_ixs]).sum()
    return E

def compute_act_error(x):
    x = asarray(x)
    N = len(x)
    xmean = x.mean()
    xvar = x.var()
    acfs = []
    tau_ints = []
    k = 0
    tau_int = 0.5
    while k < 6*tau_int:
        acf = ((x[:N-k]*x[k:]).mean() - xmean*xmean) / xvar #autocorrelatonfunction
        tau_int += acf
        N_eff = N/(2*tau_int)
        err_tau = tau_int*sqrt(12./N_eff)

        acfs.append(acf)
        tau_ints.append(tau_int)
        k += 1

    err_x = sqrt(xvar/N*2.0*tau_int)
    return xmean, xvar, err_x, tau_int, err_tau, N_eff, \
        array(acfs), array(tau_ints)

def monte_carlo_ising(L, T, N):
    V = L*L
    beta = 1.0/T

    # generate random configuration
    sigma = empty((L, L), dtype=int)
    for i in range(L):
        for j in range(L):
            if random.rand() > 0.5: 
                sigma[i,j] = -1
            else:
                sigma[i,j] = 1

    E = compute_energy(sigma)
    current_mag = sigma.sum()

    Es = []
    mags = []

    for sweep in range(N):
        for step in range(V):
            # flip single spin
            i, j = random.randint(0, L, 2)
            sigma[i,j] *= -1

            deltaE = -2*sigma[i,j]*(sigma[(i-1)%L, j] +
                                    sigma[(i+1)%L, j] +
                                    sigma[i, (j-1)%L] +
                                    sigma[i, (j+1)%L])

            if random.rand() < exp(-beta*deltaE):
                # accept move
                E += deltaE
                current_mag += 2*sigma[i,j]
            else:
                # reject move, i.e. restore spin
                sigma[i,j] *= -1

        mag = abs(current_mag)
        Es.append(E/float(V))
        mags.append(mag/float(V))

        print "T = {} {:5}/{:5}\r".format(T, sweep, N),
        sys.stdout.flush()

    E, _, err_E, tau_E, _, _, _, _ = compute_act_error(array(Es))
    mag, _, err_mag, tau_M, _, _, _, _ = compute_act_error(array(mags))
    print "\rT = {} tau_E = {} tau_M = {} E = {}+/-{} M = {}+/-{}"\
        .format(T, tau_E, tau_M, E, err_E, mag, err_mag)

    return E, err_E, mag, err_mag, sigma

L = 4
Ts = arange(1.0, 5.1, 0.1)
num_sweeps = 1000
finalstates = []

print "MC (L={})".format(L)

Ems = []
Emerrs = []
Mms = []
Mmerrs = []
finalstates = []

for T in Ts:
    E, Eerr, M, Merr, sigma = monte_carlo_ising(L, T, num_sweeps)
    Ems.append(E)
    Emerrs.append(Eerr)
    Mms.append(M)
    Mmerrs.append(Merr)
    finalstates.append(sigma)
    
figure()
subplot(211, title='Energy vs. Temperature')
errorbar(Ts, Ems, yerr=Emerrs, fmt='o-')

subplot(212, title='Magnetization vs. Temperature')
errorbar(Ts, Mms, yerr=Mmerrs, fmt='o-')

figure('States')
numplots = len(finalstates)
cols = int(ceil(sqrt(numplots)))
rows = int(ceil(numplots/float(cols)))
for i in range(numplots):
    subplot(rows, cols, i+1, title='T={}'.format(Ts[i]))
    axis('off')
    imshow(finalstates[i], interpolation='nearest')

show()
