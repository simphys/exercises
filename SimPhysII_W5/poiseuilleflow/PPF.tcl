#####################################################
#                                                   #
#               Polymer_channel.tcl                 #
#                                                   #
#               Problem 18.06.2007                  #
#                                                   #
#                                                   #
#                                                   #
#                                                   #
#####################################################

################## DESCRIPTION ######################
# "Plane Poiseuille Flow in a microchannel"         #
#####################################################

puts " "
puts "============================================="
puts "           PPF.tcl               "
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

set ext_f 0.1

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

cellsystem domain_decomposition -verlet_list

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

set density_solvent [expr 3.75*10/8]

set n_solvent [expr $density_solvent*$volume]

set solvent_id 0

set left_wall 1

set right_wall 2
################################################
#         Interaction parameters               #
################################################

# Lennard-Jones ################################

# epsilon
set lj_eps               1.0

# sigma
set lj_sig               1.0

# cutoff
set lj_cut               1.0

# shift
set lj_shift             0.0

# offset
set lj_off               0.0

# DPD friction coefficient
set gamma                  5.0

setmd temperature

set temp  1.0

set dpd_cut 1.0
################################################
#        Integration parameters                #
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
set warm_steps_all              100

# Length of warmup run for all
set warm_loops_all              500

#################################################
#        Useful parameters                      #
#################################################

set n_part         $n_solvent


set name  "PPF"
set ident "_ext_force$ext_f"

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


proc openfile_solvent {} {
    global int_steps
    global file_solvent
    global solvent_number
    global n_part
    global timestep
    global time
    global name
    global ident
    global step

    set file_solvent  [open "particle_positions_velocities.out" "a"]

    set time          [expr $timestep*$int_steps*$step]

    for {set solvent_number 0} {$solvent_number < $n_part} {incr solvent_number} {
        puts $file_solvent "$time $solvent_number [part $solvent_number print pos type v]"
    }

    close $file_solvent
}

#####################################################
#           Integration of main scheme              #
#####################################################
# Solvent particles

for {set i 0} { $i < $n_part} { incr i } {
    set posx [expr $box_x*[t_random]]
    set posy [expr $box_y*[t_random]]
    set posz [expr 1.0+8.0*[t_random]]

    set vx   [expr 0.1*[t_random]]
    set vy   [expr 0.1*[t_random]]
    set vz   [expr 0.1*[t_random]]

    part $i pos $posx $posy $posz type $solvent_id v $vx $vy $vz ext_force $ext_f 0.0 0.0
}
galileiTransformParticles

constraint plane cell -10 -10 0.0 type $left_wall
constraint plane cell -10 -10 10.0 type $right_wall


thermostat dpd $temp $gamma $dpd_cut
inter $solvent_id $left_wall lennard-jones 1.0 1.0 1.0 0.0 0.0
inter $solvent_id $right_wall lennard-jones 1.0 1.0 1.0 0.0 0.0
inter $solvent_id $left_wall tunable_slip 1.0 4.0 2.0 0.01 0 0 0
inter $solvent_id $right_wall tunable_slip 1.0 4.0 2.0 0.01 0 0 0



puts "Warm up the system with solvent ..."

puts "...Interactions:    [inter]"

    inter ljforcecap $warm_cap
    for {set step 0} {$step < $warm_loops_all} {incr step} {
        integrate $warm_steps_all
        puts "$step"
    }
puts "Integration of main scheme"
puts "... with thermostat: [thermostat]"
puts " ...and interactions: [inter]"

for {set step 0} {$step < $int_loops} {incr step} {
    puts "$step"
    integrate $int_steps
    openfile_solvent

}
puts "Finished: Simulation $name$ident"

exit

########################################################
