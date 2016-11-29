#Assigns atom type and internal of a SMILE sequence, in order to build the topology file and 3D model.

package provide topology_acyl 1.0
package require smiles_parser 1.0

namespace eval ::topology_acyl:: {
	#IC file ATOM section
	set icname 0
	set icelement 1
	set iclinker 2
	set icconnectivity 3
	set icring 4
	set icbranch 5
	set icbond 6
	set ictype 7
	set iccharge 8
	set ichydrogen 9
	
	#key in ICatom
	set ielement 0
	set ilinker 1
	set iconnectivity 2
	set iring 3
	set ibranch 4
	set ibond 5
	
	#value in ICatom
	set iname 0
	set itype 1
	set icharge 2
	set ihydrogen 3
	
	#hydrogens
	set htype 0
	set hcharge 1
	
	#variable for topology writer
    variable atoms [list]
    variable bonds [list]
    variable ICs [list]
    
    #private variable
    variable heavy_atom_key [list]
    variable heavy_atom_type [list]
    variable ICparameters [list]
    variable ICatom [list]    
    variable ICbond [list]    
    variable ICangle [list]
    variable ICdihedral [list]
    variable ICimproper [list]

}

#initialize all the variables
proc ::topology_acyl::init {} {
	variable atoms
    variable bonds
    variable ICs
    #private variable
    variable heavy_atom_key
    variable heavy_atom_type
    variable hydrogens
    variable ICparameters
    variable ICatom
    variable ICbond
    variable ICangle
    variable ICdihedral
    variable ICimproper
    #
	set atoms [list]
    set bonds [list]
    set ICs [list]
    #private variable
	set heavy_atom_key [list]
	set heavy_atom_type [list]
    set ICparameters [list]
    set ICatom [list]
    set ICbond [list]
    set ICangle [list]
    set ICdihedral [list]
    set ICimproper [list]
}



proc ::topology_acyl::noblanks {mylist} {
  set newlist [list]
  foreach elem $mylist {
    if {$elem != ""} {
      lappend newlist $elem
    }
  }

  return $newlist
}

#Reads IC parameters
#	ATOM [defines atomtype, charge and hydrogens]
#		NAME ELEMENT LINKER CONNECTIVITY RING BRANCH BOND TYPE CHARGE HYDROGENS
#			NAME: name of atom key used to identify BOND, ANGLE, DIHEDRAL, IMPROPER
#			ELEMENT: chemical element [C, N, O, ...]
#			LINKER: first atom in chain [0; 1]
#			CONNECTIVITY: number of heavy atoms bound to atom
#			RING: atom in a ring [0; 1]
#			BOND: bond type (single, double, triple bond): [1; 2; 3] (e.g C-C-C = 1:1)
#			TYPE: atom type in topology [CTL1, CTL2, ...]
#			CHARGE: atom partial charge
#			HYDROGENS: hydrogen type and charge bond to atom [HAL1:0.09, HAL2:0.09,...]
#	BOND [defines the theoretical distance between two atoms]
#		ATOM1 ATOM2 BOND (e.g. CH3 CH2 1.5)
#	ANGLE [defines the theoretical angle between three atoms]
#		ATOM1 ATOM2 ATOM3 ANGLE (e.g. CH2 CH1 CH1 120)
#	DIHEDRAL [defines the theoretical angle between four atoms]
#		ATOM1 ATOM2 ATOM3 ATOM4 ANGLE (e.g. CH2 CH2 CH2 CH3 180)
#	IMPROPER [defines the theoretical angle between four atoms]
#		ATOM1 ATOM2 ATOM3 ANGLE (e.g. CH2 CH2 CH2 CH3 120)

proc ::topology_acyl::read_ICparameters {file} {
	variable ICparameters
	variable ICatom
    variable ICbond 
	variable ICangle
	variable ICdihedral
	variable ICimproper
	
	set type "none"
	set infile [open $file "r"]
        while {[gets $infile line] >= 0} {
                switch -regexp  $line {
			^ATOM {
			set type "atom"
			continue}
			^BOND {
			set type "bond"
			continue}
			^ANGLES {
			set type "angle"
			continue}
			^DIHEDRALS {
			set type "dihedral"
			continue}
			^IMPROPER {
			set type "improper"
			continue}
		}
		if {[string $type "none"] != 0} {
	            	set linearray [split [lindex [split $line !] 0 ]]
        	        set linearray [noblanks $linearray]
                	if {[llength $linearray]} {
				switch $type {
					"atom" {
						if {[llength $linearray] == 8} {
							set key [join [lrange $linearray $::topology_acyl::icelement $::topology_acyl::icbond] -]
							set ICatom_array($key) [concat [lindex $linearray $::topology_acyl::icname] [lrange $linearray $::topology_acyl::ictype $::topology_acyl::ichydrogen]
						} else {
							puts "WARNING invalid atom: $linearray"
						}
					}	
					"bond" {
						if {[llength $linearray] == 3} {
							lassign $linearray a1 a2 b
							set ICbond_array("$a1-$a2") $b
						} else {
							puts "WARNING invalid bond: $linearray"
						}
					}
					"angle" {
						 if {[llength $linearray] == 4} {
                                		        lassign $linearray a1 a2 a3 b
		                                        set ICangle_array("$a1-$a2-$a3") $b
                		                } else {
                                		        puts "WARNING invalid angle: $linearray"
		                                }
					}
					"dihedral" {
						if {[llength $linearray] == 5} {
		                                        lassign $linearray a1 a2 a3 a4 b
                		                        set ICdihedral_array("$a1-$a2-$a3-$a4") $b
                                		} else {
			                                     puts "WARNING invalid dihedral: $linearray"
                        			}
                    }
                    "improper" {
						if {[llength $linearray] == 5} {
		                                        lassign $linearray a1 a2 a3 a4 b
                		                        set ICimproper_array("$a1-$a2-$a3-$a4") $b
                                		} else {
			                                     puts "WARNING invalid dihedral: $linearray"
                        			}

					}
        		}
  			}
		}
	}
	close $infile
	#Convert arrays to lists
	set ICatom [array get $ICatom_array]
	set ICbond [array get $ICbond_array]
	set ICangle [array get $ICangle_array]
	set ICdihedral [array get $ICdihedral_array]
	set ICimproper [array get $ICimproper_array]
	set ICparameters [list $ICatom $ICbond $ICangle $ICdihedral $ICimproper]
}	


#########################################################################
#					assign connectivity of atoms						#
#########################################################################



#parse the SMILE flatten structure and assigns the atom key
proc ::topology_acyl::assign_atom_keys {structure} {
	variable heavy_atom_key
	
	lassign $structure atom bond chains
	
	foreach a $atom {
		set aelement [lindex $a $::smiles_parser::atom_symbol];
		set alinker 0;
		if {[lindex $a $::smiles_parser::atom_id] == 0} {
			set alinker 1
		}
		set bondlist [lindex $a $::smiles_parser::atom_linkedby]
		set aconnectivity [llength $bondlist]
		set abonds [list]
		set aring 0
		set abranch 0
		if {[llength [lindex $a $::smiles_parser::atom_branches]]> 0} {
			set abranch 1
		}
		
		foreach b $bondlist {
			lappend abonds [lindex [lindex $bond $b] $::smiles_parser::bond_count]
		}
		
		set key [list $aelement $alinker $aconnectivity $aring $abranch $abonds]
		lappend heavy_atom_key $key
	}
	correct_ringatoms $structure
	correct_chirality $structure
}		


proc ::topology_acyl::detect_ringatoms {structure} {
	set index [list]
	lassign $structure atom bond chain
	
	#find all ring bond
	foreach c $chain {
		set ringatoms [lsearch -all [lindex $c $::smiles_parser::chain_list] *ringbond*]
		foreach r $ringatoms {
			#start from first atom of ring
			set sourceatom [lindex [lindex $c $::smiles_parser::chain_list] $r]
			set ringbond [lindex $sourceatom $::smiles_parser::atom_ringbonds]
			if {[llength $ringbond]>0} { 
				set bond_in_chain [list]
				#save all bonds in chain
				foreach e [lindex $c $::smiles_parser::chain_list] {
					if {[lindex $e 0] == "bond" } {
						lappend bond_in_chain [lindex $e $::smiles_parser::bond_id]
					}
				}
				#extract ringbond
				set ringbond [lindex [lindex [join $ringbond] $::smiles_parser::ringbond_bond] $::smiles_parser::bond_linking]
				set sourceatom [lindex $ringbond 0]
				lappend index $sourceatom
				set targetatom [lindex $ringbond 1]
				set currentatom $sourceatom
				
				#find path connecting source to target
				set count 0
				while { $currentatom != $targetatom } {			
					set bondlist [lindex [lindex $atom $currentatom] $::smiles_parser::atom_linkedby]
					#check if bond belongs to chain
					foreach currentbond $bondlist {
						if {[lsearch $bond_in_chain $currentbond] !=-1 } {
							set currentatom [lindex [lindex [lindex $bond $currentbond] $::smiles_parser::bond_linking] 0]
							lappend index $currentatom
							break
						}	
					}
					incr count
					if {$count>5} {
						puts "Oooops, I got lost trying to connect atom $sourceatom to atom $targetatom"
						break
					}
				}
			}
		}
	}
	
	return $index
}

proc ::topology_acyl::correct_ringatoms {structure} {
	variable heavy_atom_key
	set index [detect_ringatoms $structure]
	foreach i $index {
		set key [lindex $heavy_atom_key $i]
		set key [lreplace $key $::topology_acyl::iring $::topology_acyl::iring 1]
		set heavy_atom_key [lreplace $heavy_atom_key $i $i $key]
	}
}

proc ::topology_acyl::correct_chirality {structure} {
	variable heavy_atom_key
	lassign $structure atom bond chain
	foreach c $chain {
		set chirality [lsort -integer [concat [lsearch -all [lindex $c $::smiles_parser::chain_list] */*] [lsearch -all [lindex $c $::smiles_parser::chain_list] *\\\\*]]]
		if { [llength $chirality]>1 } {
			set c1 [lindex $chirality 0]
			set chirality [lreplace $chirality 0 0]
			foreach c2 $chirality {
				set b1 [lindex [lindex $c $::smiles_parser::chain_list] $c1]
				set b2 [lindex [lindex $c $::smiles_parser::chain_list] $c2]
				set a1 [lindex [lindex $b1 $::smiles_parser::bond_linking] 1]
				set a2 [lindex [lindex $b2 $::smiles_parser::bond_linking] 0]
				set abond [list $a1 $a2]
				if {[lsearch $bond "*$abond*"] !=-1} {
				#cis bond
					if {[lindex $b1 $::smiles_parser::bond_direction] != [lindex $b2 $::smiles_parser::bond_direction]} {
						set chirality1 [lindex [lindex [lindex $heavy_atom_key [lindex [lindex $b1 $::smiles_parser::bond_linking] 1]] $::topology_acyl::ibond] 1]
						set i [lindex [lindex $b2 $::smiles_parser::bond_linking] 0]
						set key2 [lindex $heavy_atom_key $i]
						set chirality2 [lindex $key2 $::topology_acyl::ibond]
						set chirality2 [lreplace $chirality2 1 1 [expr -1*($chirality1/abs($chirality1))*[lindex $chirality2 1]]] 
						set key2 [lreplace $key2 $::topology_acyl::ibond $::topology_acyl::ibond $chirality2]
						set heavy_atom_key [lreplace $heavy_atom_key $i $i $key2]
					}
				}
				set c1 $c2
			}
		}
	}
}


#########################################################################
#				proc for assigning topology parameters					#
#########################################################################

proc ::topology_acyl::guess_quadruplet {structure} {
	lassign $structure atom bond chain
    foreach bond $bonds {
    	lassign [lindex $bond $::smiles_parser::bond_linking] a1 a2
    	set bond1 [lindex $atom {$a1 $::smiles_parser::atom_linkedby}]
    	set bond2 [lindex $atom {$a2 $::smiles_parser::atom_linkedby}]
    	foreach b1 $bond1 {
    		set o1 [lsearch -not -inline [lindex $bond {$b1 $::smiles_parser::bond_linking}] $a1]
    		foreach b2 $bond2 {
	    		set o2 [lsearch -not -inline [lindex $bond {$b2 $::smiles_parser::bond_linking}] $a2]
				lappend quadruplet [list $o1 $a1 $a2 $o2]    		
    		}
    	}
    }
   return $quadruplet
}


proc ::topology_acyl::assign_atomtype {} {
	variable heavy_atom_key
	variable heavy_atom_type
	variable ICatom
	


}

proc ::topology_acyl::build_ICquadruplet {structure} {
	
}


proc ::topology_acyl::assign_name {prefix_C suffix_H} {
	variable heavy_atom_key
	variable atoms
    variable ICatom
    #initialize list of atom
    set atoms [list]
    array set ICatom_array $ICatom
    set index 2
	foreach key $heavy_atom_key {
		if {[info exists ICatom_array($key)]} {
			set value $ICatom_array($key)
			set nameC "${prefix}${index}"
			set typeC [lindex $value $::topology_acyl::itype]
			set chargeC [lindex $value $::topology_acyl::icharge]
			set hydrogens [lindex $value $::topology_acyl::ihydrogen]
			set hydrogenC [list]
			foreach h $hydrogens sh $suffix_H{
				lappend hydrogenC [list "H${index}$sh" [lindex $h $htype] [lindex $h $hcharge]]
			}			
			lappend atoms [list [list $nameC $typeC $chargeC] $hydrogensC]
		} else {
			puts "ERROR, unknown atom $key"
			return -1
		}
	}
	return 1
}

proc ::topology_acyl::assign_bonds {structure} {
	lassign $structure atom bond chain
	variable atoms
	varaible bonds
	foreach b $bond {
		set atoms [lindex $b $::smiles_parser::bond_linking]
		set nC [list]
		set nH [list]
		foreach a $atoms {
			set nameC [lindex [lindex [lindex atoms $a] 0] 0]
			lappend nC $nameC 
			set H [lindex [lindex atoms $a] 1]
			foreach h $H {
				lappend [list $nameC [lindex $h 0]]
			}
			lappend bonds [list $nameC $nH]
		}
	}
	set bonds [lsort -unique $bonds]
}

proc ::topology_acyl::get_atoms {} {
	
}

proc ::topology_acyl::set_ICs {} {

}