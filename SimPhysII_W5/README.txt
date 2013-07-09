1) Download Espresso

2) goto src

3) Modify myconfig-default.h: add ...

#define EXTERNAL_FORCES(is allready added)
#define CONSTRAINTS
#define DFD
#define TUNABLE_SL1P
#define LB
#define LENNARD_JONES (is allready added)

4) Then configure with ./configure

5) and finally make