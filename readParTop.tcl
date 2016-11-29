package provide readParTop 0.11

namespace eval ::readParTop:: {
	variable parameters [list]
    variable pairparameters [list]
    variable ellipsoidparameters [list]
    variable dipoleparameters [list]
    variable bondparameters [list]
    variable angleparameters [list]
    variable dihedralparameters [list]
    variable improperparameters [list]
    variable topology [list]
    variable impropertopology [list]
    variable atomtopology [list]
    variable chargetopology [list]
    variable masses [list]
    variable patches [list]
}

################# Read Topology file  #############################
proc ::readParTop::read_topology {file} {
	variable topology
	variable impropertopology
	variable atomtopology
	variable chargetopology
	variable masses
	variable patches
	set improperatoms [list]
	set improperlist [list]
	set infile [open $file "r"]
	set resname "none"
	array set improper_array $impropertopology
	while {[gets $infile line] >= 0} {
  		if {[regexp {^RESI} $line]} {
  			if { $resname != "none" } {
  				if {[info exists improper_array($resname)]} {
  				set improper_array($resname) [concat $improper_array($resname) $improperlist]
  				} else {
  					set improper_array($resname) $improperlist
  				}
	  			set improperlist [list]
  			}
  			set resname [lindex [noblanks [split $line ]] 1]
	  	}
  		if {[regexp {^IMPR} $line]} {
				set improperatoms [lreplace [noblanks [split $line ]] 0 0]
				foreach {a1 a2 a3 a4} $improperatoms {
					lappend improperlist [list $a1 $a2 $a3 $a4]
				}
		}
		if {[regexp {^ATOM} $line]} {
			set atoms [lreplace [noblanks [split $line ]] 0 0]
			set atom_array($resname-[lindex $atoms 0]) [lindex $atoms 1]
			set charge_array($resname-[lindex $atoms 0]) [lindex $atoms 2]
		}
		if {[regexp {^MASS} $line]} {
			set mass [lreplace [noblanks [split $line ]] 0 0]
			set masses_array([lindex $mass 1]) [lindex $mass 2]
		}
		if {[regexp {^PRES} $line]} {
			set resname "none"
		}
	}
	close $infile
	set improper_array($resname) [concat improper_array($resname) $improperlist]
	#PATCHES
	set infile [open $file "r"]
	set atoms [list]
	set patch "none"
	set residue 0
	while {[gets $infile line] >= 0} {
		if {[regexp {^PRES} $line]} {
			set residue 0
			if { $patch != "none" } {
  				set patches_array($patch) $atoms
				set atoms [list]
  			}
  			
			set patch [lreplace [noblanks [split $line ]] 0 0]
			set cg [lindex $patch 1]
			if { $cg == "CG" } {
				set patch [lindex $patch 0]
			}
		}
		if {[regexp {^ATOM} $line] && $residue==0} {
			lappend atoms [lreplace [noblanks [split $line ]] 0 0]
		}
		if {[regexp {^RESI} $line]} {
			if { $patch != "none" } {
  				set patches_array($patch) $atoms
				set atoms [list]
  			}
			set residue 1
			set patch "none"
		}
	}
	close $infile
	set impropertopology [array get improper_array]
	set atomtopology [concat $atomtopology [array get atom_array]]
	set chargetopology [concat $chargetopology [array get charge_array]]
	set masses [concat $masses [array get masses_array]]
	set patches [concat $patches [array get patches_array]]
	set topology [list $atomtopology $chargetopology $masses $patches $impropertopology]
}


################# Read Parameter file #############################
proc ::readParTop::read_parameters {file} {
	variable parameters
	variable pairparameters
	variable dipoleparameters
	variable ellipsoidparameters
	variable bondparameters
	variable angleparameters
	variable improperparameters
	variable dihedralparameters
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	if {[regexp {^NONBONDED} $line]} {
			read_pairs $infile
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	if {[regexp {^ELLIPSOID} $line]} {
			read_ellipsoid $infile
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	if {[regexp {^DIPOLE} $line]} {
			read_dipole $infile
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  		if {[regexp {^BONDS} $line] && ![regexp {^CMAP} $line]} {
  			read_bonds $infile
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	 if {[regexp {^ANGLES} $line]} {
  	 	read_angles $infile
  	}
  }
 	close $infile
  set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	 if {[regexp {^DIHEDRALS} $line]} {
  	 		read_dihedrals $infile
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
		if {[regexp {^IMPROPER} $line]} {
			read_impropers $infile
		}
	}
	close $infile
	set parameters [list $pairparameters $dipoleparameters $ellipsoidparameters $bondparameters $angleparameters $improperparameters $dihedralparameters]
}

proc ::readParTop::read_pairs {infile} {
	variable pairparameters
	while {[gets $infile line] >= 0} {
		if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
			break
		} {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]) [lreplace $linearray 0 0]
    		}
  	}
  }
  set pairparameters [concat $pairparameters [array get parameters]]
}

proc ::readParTop::read_ellipsoid {infile} {
	variable ellipsoidparameters
	while {[gets $infile line] >= 0} {
		if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
			break
		} {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]) [lreplace $linearray 0 0] 
    		}
  	}
  }
  set ellipsoidparameters [concat $ellipsoidparameters [array get parameters]]
}

proc ::readParTop::read_dipole {infile} {
	variable dipoleparameters
	while {[gets $infile line] >= 0} {
		if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
			break
		} {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]) [lreplace $linearray 0 0] 
    		}
  	}
  }
  set dipoleparameters [concat $dipoleparameters [array get parameters]]
}

proc ::readParTop::read_bonds {infile} {
	variable bondparameters
	while {[gets $infile line] >= 0} {
		if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
			break
		} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]-[lindex $linearray 1]) [lreplace $linearray 0 1]
    		}
  	}
  }
 set bondparameters [concat $bondparameters [array get parameters]]
}

proc ::readParTop::read_angles {infile} {
	variable angleparameters
	while {[gets $infile line] >= 0} {
		if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
			break
		} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
	    		set parameters([lindex $linearray 0]-[lindex $linearray 1]-[lindex $linearray 2]) [lreplace $linearray 0 2]
  			}
  	}
  }
  set angleparameters [concat $angleparameters [array get parameters]]
}

proc ::readParTop::read_dihedrals {infile} {
	variable dihedralparameters 
	while {[gets $infile line] >= 0} {
			if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
				break
			} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]-[lindex $linearray 1]-[lindex $linearray 2]-[lindex $linearray 3]) [lreplace $linearray 0 3]
  			}
  	}
  }
  set dihedralparameters [concat $dihedralparameters [array get parameters]]
}

proc ::readParTop::read_impropers {infile} {
	variable improperparameters 
	while {[gets $infile line] >= 0} {
			if {[regexp {^DIPOLE} $line] || [regexp {^ELLIPSOID} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line] || [regexp {^END} $line]} {
				break
			} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]-[lindex $linearray 1]-[lindex $linearray 2]-[lindex $linearray 3]) [lreplace $linearray 0 3]
  			}
  	}
  }
  set improperparameters [concat $improperparameters [array get parameters]]
}

################# Read Parameter file #############################


proc ::readParTop::format_parameters {file out} {
	read_parameters $file
	set outfile [open $out w]
	puts $outfile "!Reformated $file parameter file"
	writebonds $outfile
	puts $outfile "\n"
	writeangles $outfile
	puts $outfile "\n"
	writedihedrals $outfile
	puts $outfile "\n"
	writeimpropers $outfile
	puts $outfile "\n"	
	writepairs $outfile
	puts $outfile "\n"	
	close $outfile
}


proc ::readParTop::writebonds {outfile} {
	variable bondparameters
	puts $outfile "BONDS"
	puts $outfile "!"
	puts $outfile "!V(bond) = Kb(b - b0)**2"
	puts $outfile "!"
	puts $outfile "!atom type Kb\t\tb0"
	foreach {type par} $bondparameters {
		set type [split $type "-"]
		puts $outfile "$type\t\t[join $par "\t\t"]"
	}
}


proc ::readParTop::writeangles {outfile } {
	variable angleparameters
	puts $outfile "ANGLES"
	puts $outfile "!"
	puts $outfile "!V(angle) = Ktheta(Theta - Theta0)**2+UB"
	puts $outfile "!"
	puts $outfile "!atom types\t\tTheta0\t\tKtheta1\t\tKtheta2\t\tKtheta3"
	foreach {type par} $angleparameters {
		set type [split $type "-"]
		lassign $par k1 t  u b
		if {$u=="" || $b==""} {
            		set u 0.0
            		set b 0.0
        }
        set par [list $t $k1 0.0 0.0 $u $b]
		puts $outfile "$type\t\t[join $par "\t\t"]"
	}
}

proc ::readParTop::writedihedrals {outfile } {
	variable dihedralparameters
	puts $outfile "DIHEDRALS"
	puts $outfile "!"
	puts $outfile "!V(dihedral) = sum \[kn cos(phin -dn)\]"
	puts $outfile "!"
	puts $outfile "!atom types\t\tk1\t\tk2\t\tk3\t\ttheta1\t\ttheta2\t\ttheta3"
	foreach {type par} $dihedralparameters {
		set type [split $type "-"]
		lassign $par k n phi
		switch $n {
			1
				{set par [list $k 0.0 0.0 $phi 0.0 0.0]}
			2
				{set par [list 0.0 $k 0.0 0.0 $phi 0.0]}
			default
				{set par [list 0.0 0.0 $k 0.0 0.0 $phi]}
		}
		puts $outfile "$type\t\t[join $par "\t\t"]"
	}
}

proc ::readParTop::writeimpropers {outfile } {
	variable improperparameters
	puts $outfile "IMPROPER"
	puts $outfile "!"
	puts $outfile "!V(improper) = kn(phio-phin)**n"
	puts $outfile "!"
	puts $outfile "!atom types\t\tk1\t\tk2\t\tk3\t\ttheta1\t\ttheta2\t\ttheta3"
	foreach {type par} $improperparameters {
		set type [split $type "-"]
		lassign $par k n phi
		switch $n {
			1
				{set par [list $k 0.0 0.0 $phi 0.0 0.0]}
			2
				{set par [list 0.0 $k 0.0 0.0 $phi 0.0]}
			default
				{set par [list 0.0 0.0 $k 0.0 0.0 $phi]}
		}
		puts $outfile "$type\t\t[join $par "\t\t"]"
	}
}

proc ::readParTop::writepairs {outfile } {
	variable pairparameters
	puts $outfile "NONBONDED"
	puts $outfile "!"
	puts $outfile "!V(nonbonded) = LJ"
	puts $outfile "!"
	puts $outfile "!!atom\t\tignored\t\tepsilon\t\tsigma"
	foreach {type par} $pairparameters {
		set type [split $type "-"]
		lassign $par i e r2 i2 e2 r
		set par [list 0.0 [expr -1*$e] [expr $r/pow(2, 1/6)]]
		puts $outfile "$type\t\t[join $par "\t\t"]"
	}
}
