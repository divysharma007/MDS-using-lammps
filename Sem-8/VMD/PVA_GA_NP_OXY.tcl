#Loading equilibrated lammps data file to form cross-links
topo readlammpsdata PVA_GA_NC_OXY_LARGE_DIHEDRALS_INPUT_WITH_CHARGE_EQUILIBRATED_BEFORE_LINKING.data  

#system is loaded in unwrapped condition, thus wrapping it a cube
pbc wrap

# C1 carbon of GA                             
set sel14 [atomselect top "type 14"]    
# O of PVA 
set sel11 [atomselect top "type 11"]        
set cutoff 5.0


# Array of indexes of type 14 atoms
set idx14 [$sel14 get index]   
# Array of indexes of type 11 atoms    
set idx11 [$sel11 get index]  


# Size (length) of idx14 array      
set num14 [$sel14 num]          
# Size (length) of idx11 array   
set num11 [$sel11 num]             


for {set i 0} {$i < $num14} {incr i} {

    # Getting index of atom present at the ith position in the idx14 array
    set idx1 [lindex $idx14 $i]   
    # Getting atom corresponding to idx1                
    set atom14 [atomselect top "index $idx1"]       
    set OH1 -1

    # Getting coordinates of atom14 (C1 of GA)
    set pos1 [$atom14 get {x y z}]  

    # Getting index list, A list of atoms connected to this atom                
    set bonds [$atom14 getbonds]                   
    set b [lindex $bonds 0]

    # Carbon connected to three atoms, 1 double bonded O, 1 C2, 1 H
    set bonded_atom1 [lindex $b 0]                  
    set bonded_atom2 [lindex $b 1]
    set bonded_atom3 [lindex $b 2]

    # Getting atoms from the indexes
    set BA1 [atomselect top "index $bonded_atom1"] 
    set BA2 [atomselect top "index $bonded_atom2"]
    set BA3 [atomselect top "index $bonded_atom3"]

    # Getting type of the atoms connected to atom14
    set BA1_type [$BA1 get type]                   
    set BA2_type [$BA2 get type]
    set BA3_type [$BA3 get type]

    # Finding O connected to atom14 (C1 of GA)
    if {$BA1_type == "17"} {                       
        set O $bonded_atom1
    } elseif {$BA2_type == "17"} {
        set O $bonded_atom2
    } else {
        set O $bonded_atom3
    }

        # Searching for O of PVA (O1/ type 11)               
        for {set j 0} {$j < $num11} {incr j} {   

            # Getting index of atom present at the ith position in the idx11 array               
            set idx2 [lindex $idx11 $j]     
            # Getting atom corresponding to idx2                   
            set atom11 [atomselect top "index $idx2"]  

            # Getting types of atom11 (O) and atom14 to see their cross-linking state        
            set type11 [$atom11 get type]                      
            set type14 [$atom14 get type]

            # type11==11 implies fresh O, type14==14 implies fresh C
            if {$type11 == "11" && $type14=="14"} {  

                # Getting coordinates of atom11 (O of PVA)          
                set pos2 [$atom11 get {x y z}]  

                # Getting distance between atom11 and atom14     
                set dist [measure bond [list $idx1 $idx2]]     
                # if the distance is less than the cutoff, we begin bond-formation procedure

                if {[expr $dist <$cutoff]} {                    
                    puts "first crosslinking dist- $dist"
                    set bonds [$atom11 getbonds]
                    set b [lindex $bonds 0]
                    set bonded_atom1 [lindex $b 0]
                    set bonded_atom2 [lindex $b 1]
                    set BA1 [atomselect top "index $bonded_atom1"]
                    set BA2 [atomselect top "index $bonded_atom2"]
                    set BA1_type [$BA1 get type]
                    set BA2_type [$BA2 get type] 
                     # H connected to this atom11 is found
                    if {$BA1_type == "8"} {                    
                        set OH1 $bonded_atom1
                    } else {
                          set OH1 $bonded_atom2
                    }
    puts " single- H- $OH1 O- $O"

    # H connected to atom11(O of PVA) is connected to O of GA
    topo addbond $OH1 $O   

    # Bond between atom11 and this H is deleted         
    topo delbond $OH1 $idx2  

    # Bond between atom11 and atom14  
    topo addbond $idx1 $idx2 

    # Changing atom type of atom14 to 19 which represents single cross-linked C1 in GA      
    $atom14 set type 19   

    # Changing atom type of atom11 to 18 which represents O of PVA connected to C1 in GA          
    $atom11 set type 18    

    set atom8 [atomselect top "index $OH1"]
    set atom17 [atomselect top "index $O"]

    # Changing type of H connected to atom11 to 21 as it is now connected to O of GA
    $atom8 set type 21  

    # Changing type of O connected to atom14 to 22 representing single bonded O        
    $atom17 set type 22            
                }
            } elseif {$type11 == "11" && $type14=="19"} { 

                # single-crosslinked C1 of GA and a fresh O of PVA       
                # Similiar steps formed as done in single cross-linking procedure

                set pos2 [$atom11 get {x y z}]                    
                set dist [measure bond [list $idx1 $idx2]]
                if {[expr $dist <$cutoff]} {
                    puts "Second crosslinking dist- $dist"
                    puts "prev OH1- $OH1 O- $O "
                    set bonds [$atom11 getbonds]
                    set b [lindex $bonds 0]
                    # puts " D  for b in "
                    set bonded_atom1 [lindex $b 0]
                    set bonded_atom2 [lindex $b 1]
                    set BA1 [atomselect top "index $bonded_atom1"]
                    set BA2 [atomselect top "index $bonded_atom2"]
                    set BA1_type [$BA1 get type]
                    set BA2_type [$BA2 get type]   
                    set OH2 -1
                    if {$BA1_type == "8"} {
        set OH2 $bonded_atom1
    } else {
        set OH2 $bonded_atom2
    }
    puts "cur OH2- $OH2"

    # Bond between O of PVA and H connected to it is broken
    topo delbond $OH2 $idx2   
    # Bond is created between atom14 and atom11, forming Second cross-link
    topo addbond $idx1 $idx2 

    # Bond between H of PVA and O in GA which was formed during the single cross-linking procedure is now deleted
    topo delbond  $OH1 $O   

    # Bond between O connected to C1 of GA and C1(single cross-link) is removed
    topo delbond $idx1 $O  

    # Changing atom type from 19 to 20 signifies double cross-linked C1 of GA
    $atom14 set type 20  

    # Changing atom type of atom11 to 18 which represents O of PVA connected to C1 in GA    
    $atom11 set type 18  

    # connects OH2 ( H connected to O of PVA in double cross-link procedure) and O which was connected to GA, forming water
    #topo addbond $OH2 $O  

    # Changing type of O of water to later remove them
    set atom17 [atomselect top "index $O"]
    $atom17 set type 23   

    # Changing type of H of water to later remove them
    set atom17 [atomselect top "index $OH2"]
    $atom17 set type 24     
    set atom17 [atomselect top "index $OH1"]
    $atom17 set type 24

    # Double cross-link procedure over, break isn't really necessary
    break                      

                }
            }
        }
    }

# Now chaning names of the atom types as when we form lammps data file it reorders types lexicographically

for {set i 1} {$i <= 24} {incr i} {
set sel [atomselect top "type $i"]
set newtypename -1
if {$i<=9} {
$sel set type "T1$i" 
set newtypename "T1$i"  
} elseif { $i<=19 } {
set newtype [expr $i%10]
$sel set type "T2$newtype"  
set newtypename "T2$newtype"
} else {
set newtype [expr $i%10]
$sel set type "T3$newtype"  
set newtypename "T3$newtype"
}
puts "Type $i changed to $newtypename"
}

# Guessing angles, dihedrals and then retyping their types

topo guessangles
topo guessdihedrals
topo retypebonds
topo retypeangles
topo writelammpsdata PVA_GA_NC_OXY_LARGE_DIHEDRALS_INPUT_WITH_CHARGE_AFTER_LINKING_WITH_WATER_LATEST_temp.data
set sel [atomselect top "type T29"]
set type19 [$sel num]
set sel [atomselect top "type T30"]
set type20 [$sel num]
puts "Number of single cross-linked C- $type19\nNumber of double cross-linked C- $type20\n"