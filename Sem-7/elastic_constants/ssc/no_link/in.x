units			real
dimension		3
processors		* * *
boundary		p p p
atom_style		full
pair_style		lj/cut 10.0
bond_style		harmonic
angle_style 	harmonic
read_data		PVA4_equilibrated.data
variable simname index ssc_not_linked_latest_1ns
variable tmp equal "lx"
variable L0 equal ${tmp}
variable strain equal "(lx - v_L0)/v_L0"
variable p1 equal "v_strain"
variable p2 equal "-pxx/10000*1.01325"
variable p3 equal "-pyy/10000*1.01325"
variable p4 equal "-pzz/10000*1.01325"
variable p5 equal "lx"
variable p6 equal "ly"
variable p7 equal "lz"
variable p8 equal "temp"
fix		1 all npt temp 300 300 50 y 1 1 100 z 1 1 100
fix		2 all deform 1 x erate 1e-6 units box remap x
fix        3 all ave/time 1 100 100 v_p2 v_p3
variable avg_stress equal f_3[1]
fix def1 all print 100 "${p1} ${p2} ${avg_stress}" file ${simname}.txt screen no
thermo_style	custom step temp v_p2 v_avg_stress lx ly lz vol v_strain pxx
thermo          100
timestep	1
run		1000000
unfix 1
unfix 2
unfix def1

print "All done"