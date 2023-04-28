# Initialization
log MSD_link_win.lammps
units			real
dimension		3
processors		* * *
boundary		p p p
atom_style		full
pair_style		lj/cut 10.0
bond_style		harmonic
angle_style 	harmonic
dihedral_style harmonic
read_data		PVA_GA_NP_H2O_EQUILIBRATED_WITH_CHARGE_AFTER_LINKING.data


timestep 1

group  water type 11
compute         msd water msd
fix             val water vector 100 c_msd[4]

# First nano particle 6217:6420
group         NP_1   id   6217:6420
compute     msdNP_1 NP_1 msd
fix        AVE_MSD_NP1 NP_1 ave/time 1 1 100 c_msdNP_1[4] file msddumpfile_msdNP_1


# Second nano particle 6421:6624
group         NP_2   id   6421:6624
compute     msdNP_2 NP_2 msd
fix        AVE_MSD_NP2 NP_2 ave/time 1 1 100 c_msdNP_2[4] file msddumpfile_msdNP_2



#output
thermo			100
thermo_style	custom step c_msd[4] c_msdNP_1[4]  c_msdNP_2[4]
dump			1 water custom 100 MSD_link_win.lammpstrj type id xu yu zu # unwrapped coordinates


# nvt for 10 ns
fix			1 all nvt temp 300 300 100
run			10000000
unfix			1
#final 
write_data		PVA_GA_NP_H2O_EQUILIBRATED_WITH_CHARGE_AFTER_LINKING_MSD.data

#parameter
neighbor		2.0 bin
