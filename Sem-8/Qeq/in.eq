# Initialization
units			real
dimension		3
processors		* * *
boundary		p p p
atom_style		full
pair_style		lj/cut 10.0
bond_style		harmonic
angle_style 	harmonic
dihedral_style harmonic
read_data		PVA_GA_NC_OXY_LARGE_DIHEDRALS_AFTER_LINKING.data

#output

thermo			100
thermo_style	custom step temp pe press vol lx density

fix 1 all qeq/point 1 100 1.0e-6 200 param1.qeq
run 0
unfix 1

#final 
write_data		PVA_GA_NC_OXY_LARGE_DIHEDRALS_INPUT_WITH_CHARGE_AFTER_LINKING.data.data

#parameter
neighbor		2.0 bin
