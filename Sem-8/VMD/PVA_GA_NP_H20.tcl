#Loading equilibrated lammps data file to form cross-links
topo readlammpsdata PVA_GA_NP_H2O_INPUT_WITH_CHARGE.data  

#system is loaded in unwrapped condition, thus wrapping it a cube
pbc wrap

# C1 carbon of GA                             
set sel15 [atomselect top "type 15"]    
# O of PVA 
set sel12 [atomselect top "type 12"]        
set cutoff 5.0


# Array of indexes of type 15 atoms
set idx15 [$sel15 get index]   
# Array of indexes of type 12 atoms    
set idx12 [$sel12 get index]  


# Size (length) of idx15 array      
set num15 [$sel15 num]          
# Size (length) of idx12 array   
set num12 [$sel12 num]             


for {set i 0} {$i < $num15} {incr i} {

    # Getting index of atom present at the ith position in the idx15 array
    set idx1 [lindex $idx15 $i]   
    # Getting atom corresponding to idx1                
    set atom15 [atomselect top "index $idx1"]       
    set OH1 -1

    # Getting coordinates of atom15 (C1 of GA)
    set pos1 [$atom15 get {x y z}]  

    # Getting index list, A list of atoms connected to this atom                
    set bonds [$atom15 getbonds]                   
    set b [lindex $bonds 0]

    # Carbon connected to three atoms, 1 double bonded O, 1 C2, 1 H
    set bonded_atom1 [lindex $b 0]                  
    set bonded_atom2 [lindex $b 1]
    set bonded_atom3 [lindex $b 2]

    # Getting atoms from the indexes
    set BA1 [atomselect top "index $bonded_atom1"] 
    set BA2 [atomselect top "index $bonded_atom2"]
    set BA3 [atomselect top "index $bonded_atom3"]

    # Getting type of the atoms connected to atom15
    set BA1_type [$BA1 get type]                   
    set BA2_type [$BA2 get type]
    set BA3_type [$BA3 get type]

    # Finding O connected to atom15 (C1 of GA)
    if {$BA1_type == "18"} {                       
        set O $bonded_atom1
    } elseif {$BA2_type == "18"} {
        set O $bonded_atom2
    } else {
        set O $bonded_atom3
    }

        # Searching for O of PVA (O1/ type 11)               
        for {set j 0} {$j < $num12} {incr j} {   

            # Getting index of atom present at the ith position in the idx12 array               
            set idx2 [lindex $idx12 $j]     
            # Getting atom corresponding to idx2                   
            set atom12 [atomselect top "index $idx2"]  

            # Getting types of atom12 (O) and atom15 to see their cross-linking state        
            set type12 [$atom12 get type]                      
            set type15 [$atom15 get type]

            # type12==12 implies fresh O, type15==15 implies fresh C
            if {$type12 == "12" && $type15=="15"} {  

                # Getting coordinates of atom12 (O of PVA)          
                set pos2 [$atom12 get {x y z}]  

                # Getting distance between atom12 and atom15     
                set dist [measure bond [list $idx1 $idx2]]     
                # if the distance is less than the cutoff, we begin bond-formation procedure

                if {[expr $dist <$cutoff]} {                    
                    puts "first crosslinking dist- $dist"
                    set bonds [$atom12 getbonds]
                    set b [lindex $bonds 0]
                    set bonded_atom1 [lindex $b 0]
                    set bonded_atom2 [lindex $b 1]
                    set BA1 [atomselect top "index $bonded_atom1"]
                    set BA2 [atomselect top "index $bonded_atom2"]
                    set BA1_type [$BA1 get type]
                    set BA2_type [$BA2 get type] 
                    # H connected to this atom12 is found
                    if {$BA1_type == "9"} {                    
                        set OH1 $bonded_atom1
                    } else {
                          set OH1 $bonded_atom2
                    }
    puts " single- H- $OH1 O- $O"

    # H connected to atom12(O of PVA) is connected to O of GA
    topo addbond $OH1 $O   

    # Bond between atom12 and this H is deleted         
    topo delbond $OH1 $idx2  

    # Bond between atom12 and atom15  
    topo addbond $idx1 $idx2 

    # Changing atom type of atom15 to 20 which represents single cross-linked C1 in GA      
    $atom15 set type 20   

    # Changing atom type of atom12 to 19 which represents O of PVA connected to C1 in GA          
    $atom12 set type 19    

    set atom9 [atomselect top "index $OH1"]
    set atom18 [atomselect top "index $O"]

    # Changing type of H connected to atom12 to 22 as it is now connected to O of GA
    $atom9 set type 22  

    # Changing type of O connected to atom15 to 23 representing single bonded O        
    $atom18 set type 23            
                }
            } elseif {$type12 == "12" && $type15=="20"} { 

                # single-crosslinked C1 of GA and a fresh O of PVA       
                # Similiar steps formed as done in single cross-linking procedure

                set pos2 [$atom12 get {x y z}]                    
                set dist [measure bond [list $idx1 $idx2]]
                if {[expr $dist <$cutoff]} {
                    puts "Second crosslinking dist- $dist"
                    puts "prev OH1- $OH1 O- $O "
                    set bonds [$atom12 getbonds]
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
    # Bond is created between atom15 and atom12, forming Second cross-link
    topo addbond $idx1 $idx2 

    # Bond between H of PVA and O in GA which was formed during the single cross-linking procedure is now deleted
    topo delbond  $OH1 $O   

    # Bond between O connected to C1 of GA and C1(single cross-link) is removed
    topo delbond $idx1 $O  

    # Changing atom type from 20 to 21 signifies double cross-linked C1 of GA
    $atom15 set type 21  

    # Changing atom type of atom12 to 19 which represents O of PVA connected to C1 in GA    
    $atom12 set type 19  

    # connects OH2 ( H connected to O of PVA in double cross-link procedure) and O which was connected to GA, forming water
    #topo addbond $OH2 $O  

    # Changing type of O of water to later remove them
    set atom18 [atomselect top "index $O"]
    $atom18 set type 24   

    # Changing type of H of water to later remove them
    set atom18 [atomselect top "index $OH2"]
    $atom18 set type 25     
    set atom18 [atomselect top "index $OH1"]
    $atom18 set type 25

    # Double cross-link procedure over, break isn't really necessary
    break                      

                }
            }
        }
    }

# Now chaning names of the atom types as when we form lammps data file it reorders types lexicographically

for {set i 1} {$i <= 25} {incr i} {
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
topo writelammpsdata temp.data
set sel [atomselect top "type T30"]
set type20 [$sel num]
set sel [atomselect top "type T31"]
set type21 [$sel num]
puts "Number of single cross-linked C- $type20\nNumber of double cross-linked C- $type21\n"