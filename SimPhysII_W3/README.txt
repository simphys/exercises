Files
-------------------------------------------------------------------
- spc216.pdb: pre equilibrated structuture of 216 water mulecules

- topol.top: itp-files are the force files (here spc.itp / spce.itp / tip3p.itp)

- grompp.mdp:

integrator = md (velocity verlet)
tinit= 0
dt = 0.002 (size of time step)
nsteps = 250000

Steps to do for simulation
-------------------------------------------------------------------
1) vmd spc216.pdb
2) vim topol.top
3) vim gromp.mdp: change nsteps = 250000
4) create directory with all files for each watermodell
5) groupp -v
6) mdrun -v

Analyzation
-------------------------------------------------------------------
- radio distribution function
- number of H-bonds (hydrogen to oxygen distance is bellow 0.35 nm, angle is less than 30 deg)
- mean squared displacement and diffusion coefficient
< DELTA r(t)^2 > = 6 D t

1) g_rdf -n index
2) g_hbond
3) g_msd -n index

have a look at the resulting files witrh gnu-plot

Gromacs
-------------------------------------------------------------------
download "gromacs-4.5.5.tar.gz" from www.gromacs.org
tar -xzf gromacs-4.5.5.tar.gz
./configure
make

or

use the path /usr/local64/gromacs4.0.4-gcc4.3-static/bin/ ...

vim .basharc:
PATH = $PATH:/usr/local64/gromacs4.0.4-gcc4.3-static/bin oder source /usr/local64/gromacs4.0.4-gcc4.3-static/bin
export PATH
