#uniaxial tensile test of graphene

##---------------INITIALIZATION-------------------------------

units			real
dimension		3
processors		* * *
boundary		p p p
atom_style		full
pair_style		lj/cut 10.0
bond_style		harmonic
angle_style 	harmonic
read_data		after_cross_linking_correct.data



timestep 	1
variable   ts equal 1
variable   temp equal lx
variable   Lx equal ${temp}
variable   Vol equal "vol*1e-30"

##---------------COMPUTES-------------------------------------

compute peratom all stress/atom NULL
compute fx all reduce sum c_peratom[1] c_peratom[2] c_peratom[3] 
 

##---------------RELAXATION--------------------------------------

reset_timestep 0

##---------------DEFORMATION--------------------------------------
#fix		1 all npt temp 300 300 1 y 0 0 1 z 0 0 1 drag 2
fix        2 all ave/time 1 100 100 c_fx[1] c_fx[2] c_fx[3]
variable   srate equal 1.0e11
variable   srate1 equal "v_srate / 1.0e15"
variable   conv equal 101325e-9
fix		4 all deform 1 x erate ${srate1} units box remap x
##---------------THERMO-OUTPUTS--------------------------------------
variable strain equal "(lx-v_Lx)/v_Lx"
variable sigmavolume equal "c_fx[1]*v_conv/vol"
variable sigmaxx equal "-pxx*v_conv"
variable avg_sigmaxx equal "f_2[1]*v_conv/vol"


thermo 	100
thermo_style custom step v_strain v_sigmavolume v_sigmaxx v_avg_sigmaxx lx vol pxx


run            200000
