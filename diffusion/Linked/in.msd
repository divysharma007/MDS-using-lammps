# Initialization
units			real
dimension		3
processors		* * *
boundary		p p p
atom_style		full
pair_style		lj/cut 10.0
bond_style		harmonic
angle_style 	harmonic
read_data		PVA4_OXY10_GA20_linked_equilibrated.data

#output
group O2 type 9
compute         msd O2 msd
fix             val O2 vector 100 c_msd[4]
variable        fitslope equal slope(f_val)
thermo			100
thermo_style	custom step c_msd[4] v_fitslope
timestep 1

# nvt for 10 ns
fix			1 all nvt temp 300 300 100
run			10000000
unfix			1

#final 
write_data		PVA4_OXY10_GA20_linked_MSD.data

#parameter
neighbor		2.0 bin