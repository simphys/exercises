################## DESCRIPTION ########################
# "Lattice-Boltzmann simulations of a single particle"#
#######################################################

puts " "
puts "============================================="
puts "           lb_particle.tcl               "
puts "============================================="
puts " "

puts "Program Information: \n[code_info]\n"

###################################################
#                   Parameters                    #
###################################################

# precision of data output
set tcl_precision     6

# initialize random number generator
set t_random_seed     [expr int(rand()*99999)^[pid]]

# System identification:

set box_l 10.0

# Outout for VMD: yes or not
set vmd_output "no"

# Checkpoint write-out after int_steps
set chkpoint             1000

# Checkpointing for warmup
set warm_check           5000

set all_warm_check       1000

set volume               [expr $box_l*$box_l*$box_l]

set box_x                $box_l

set box_y                $box_l

set box_z                $box_l

# reuse of Verlet lists
setmd verlet_reuse

# skin depth
setmd skin               0.4

# box length
setmd box_l              $box_l $box_l $box_l

# volume of system
set volume               [expr $box_x*$box_y*$box_z]

# periodic boundary conditions
setmd periodic           1 1 1

cellsystem domain_decomposition -no_verlet_list

# maximal number of cells: 10^3
setmd max_num_cells      1000

set max_cells            1000

# minimal number of cells: 2^3
setmd min_num_cells      8

set min_cells            8

# cell size
setmd cell_size

# cell grid
setmd cell_grid

################################################
#         Solvent parameters                   #
################################################

set particle_id 0

################################################
#         Interaction parameters               #
################################################

# Lennard-Jones force cap for warm up
set warm_cap                    10000

# Lennard-Jones force cap for integration
set int_cap                     10000

# length of time step
setmd time_step                 0.01

set timestep                    0.01

# Length of integration scheme
set int_steps                   100

# Length of integration run
set int_loops                   1000

# Length of warmup for all
set warm_steps                  100

# Length of warmup run for all
set warm_loops                  1000

# Length of warmup for all
set warm_steps_all              100

# Length of warmup run for all
set warm_loops_all              500

set n_part                      1
#################################################
#        Useful parameters                      #
#################################################

set name  "LBP"
set ident ""

# print out of parameters
puts "...skin:         [setmd skin]"
puts "...PBC:          [setmd periodic]"
puts "...max_cells:    [setmd max_num_cells]"
puts "...min_cells:    [setmd min_num_cells]"
puts "...cell_grid:    [setmd cell_grid]"
puts "...cell size:    [setmd cell_size]"
puts "...verlet_reuse: [setmd verlet_reuse]"
puts "...random seed   $t_random_seed"
puts "...node grid:    [setmd node_grid]"

proc openfile_particle {} {
    global int_steps
    global file_particle
    global timestep
    global time
    global step

    set file_particle [open "particle_diffusion.out" "a"]

    set time [expr $timestep*$int_steps*$step]

    puts $file_particle "$time [part 0 print pos type v]"

    close $file_particle
}

#####################################################
#           Integration of main scheme              #
#####################################################
# Particle

for {set i 0} { $i < $n_part} { incr i } {
    set posx [expr $box_x*[t_random]]
    set posy [expr $box_y*[t_random]]
    set posz [expr 1.0+8.0*[t_random]]

    set vx   [expr 0.1*[t_random]]
    set vy   [expr 0.1*[t_random]]
    set vz   [expr 0.1*[t_random]]

    part $i pos $posx $posy $posz type $particle_id v $vx $vy $vz
}

lbfluid agrid 1.0 dens 3.0 visc 1.0 tau 0.01 friction 1.0
thermostat lb 1.0

puts "Warm up the system with solvent ..."

puts "...Interactions:    [inter]"

    for {set step 0} {$step < $warm_loops_all} {incr step} {
        integrate $warm_steps_all
        puts "$step"
    }
puts "Integration of main scheme"
puts "... with thermostat: [thermostat]"
puts " ...and interactions: [inter]"


for {set step 0} {$step < $int_loops} {incr step} {
    integrate $int_steps
    #puts "[part 0 print pos type v]"
    openfile_particle
}
puts "Finished: Simulation $name$ident"

exit

########################################################
