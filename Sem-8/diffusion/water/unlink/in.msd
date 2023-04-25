# Initialization
log MSD_unlink_win.lammps
units			real
dimension		3
processors		* * *
boundary		p p p
atom_style		full
pair_style		lj/cut 10.0
bond_style		harmonic
angle_style 	harmonic
dihedral_style harmonic
read_data		PVA_NP_H2O_EQUILIBRATED_WITH_CHARGE.data


timestep 1

group  water type 11 6
compute         msd water msd
fix             val water vector 100 c_msd[4]

# First nano particle 5617:5820
group         NP_1   id   5617:5820
compute     msdNP_1 NP_1 msd
fix        AVE_MSD_NP1 NP_1 ave/time 1 1 100 c_msdNP_1[4] file msddumpfile_msdNP_1


# Second nano particle 5821:6024
group         NP_2   id   5821:6024
compute     msdNP_2 NP_2 msd
fix        AVE_MSD_NP2 NP_2 ave/time 1 1 100 c_msdNP_2[4] file msddumpfile_msdNP_2



#output
thermo			100
thermo_style	custom step c_msd[4] c_msdNP_1[4]  c_msdNP_2[4]
dump			1 water custom 100 MSD_unlink_win.lammpstrj type id xu yu zu # unwrapped coordinates


# nvt for 10 ns
fix			1 all nvt temp 300 300 100
run			10000000
unfix			1
#final 
write_data		PVA_NP_H2O_EQUILIBRATED_WITH_CHARGE_MSD.data

#parameter
neighbor		2.0 bin
