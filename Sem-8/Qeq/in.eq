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
read_data		PVA_GA_NP_H2O.data

#output

thermo			100
thermo_style	custom step temp pe press vol lx density

fix 1 all qeq/point 1 100 1.0e-6 200 param2.qeq
run 0
unfix 1

#final 
write_data		PVA_GA_NP_H2O_INPUT_WITH_CHARGE.data

#parameter
neighbor		2.0 bin
