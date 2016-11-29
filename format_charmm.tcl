proc format_charmm {file output} {
	lassign [read_parameters $file] pairparameters bondparameters angleparameters dihedralparameters improperparameters
	set fl [open $output w]
	puts $fl "*>>>>>> CHARMM 27  Force field for LAMMPS <<<<<<<<<"
	puts $fl ""
	write_bonds $fl $bondparameters
	puts $fl ""
	write_angles $fl $angleparameters
	puts $fl ""
	write_dihedrals $fl $dihedralparameters
	puts $fl ""
	write_impropers $fl $improperparameters
	puts $fl ""
	write_pair $fl $pairparameters
	
	close $fl
}

proc write_bonds {fl bondparameters} {
	puts $fl "BONDS"
	puts $fl "!" 
	puts $fl "!V(bond) = Kb(b - b0)**2"
	puts $fl "!"
	puts $fl "!Kb: kcal/mole/A**2"
	puts $fl "!b0: A"
	puts $fl "!"
	puts $fl "!atom type Kb          b0"
	array set parameter_array $bondparameters
	foreach n [array names parameter_array] {
		set nn [split $n "-"]
		lassign $nn a1 a2
		lassign $parameter_array($n) k b
		puts $fl [format "%s\t%s\t %6.2f %6.2f" $a1 $a2 $k $b]
	}
}

proc write_angles {fl angleparameters} {
	puts $fl "ANGLES"
	puts $fl "!" 
	puts $fl "!V(angle) = Ktheta(Theta - Theta0)**2 + Ktheta(Theta - Theta0)**3 +Ktheta(Theta - Theta0)**4"
	puts $fl "!"
	puts $fl "!V(Urey-Bradley) = Kub(S - S0)**2"
	puts $fl "!atom types     Theta0    Ktheta1   Ktheta2   Ktheta3 Kub     S0"
	array set parameter_array $angleparameters
	foreach n [array names parameter_array] {
		set nn [split $n "-"]
		lassign $nn a1 a2 a3
		lassign $parameter_array($n) k t u b
		if {$u=="" || $b==""} {
			set u 0.0
			set b 0.0
		}
		puts $fl [format "%5s %5s %5s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" $a1 $a2 $a3 $t $k 0 0 $u $b]
	}
}

proc write_dihedrals {fl dihedralparameters} {
	puts $fl "DIHEDRALS"
	puts $fl "!" 
	puts $fl "!V(dihedral) = Kchi1(1 + cos(chi - chi1)) + Kchi2(1 + cos(2(chi) - chi2))) + Kchi3(1 + cos(3(chi) - chi3))"
	puts $fl "!"
	puts $fl "!atom types     Kchi1     Kchi2      Kchi3      chi1   chi2   chi3"
	array set parameter_array $dihedralparameters
	foreach n [array names parameter_array] {
		set nn [split $n "-"]
		lassign {0.0 0.0 0.0 0.0 0.0 0.0} Kchi1 Kchi2  Kchi3 chi1 chi2 chi3
		lassign $nn a1 a2 a3 a4
		lassign $parameter_array($n) k n chi
		if {$n==1} {
			set Kchi1 $k
			set chi1 $chi
		} elseif {$n==2} {
			set Kchi2 $k
			set chi2 $chi
		} elseif {$n==3} {
			set Kchi3 $k
			set chi13 $chi
		} else {
			continue
		}
		puts $fl [format "%5s %5s %5s %5s %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f" $a1 $a2 $a3 $a4 $Kchi1 $Kchi2 $Kchi3 $chi1 $chi2 $chi3]
	}
}

proc write_impropers {fl improperparameters} {
	puts $fl "IMPROPER"
	puts $fl "!" 
	puts $fl "!V(improper) = Kpsi1(psi - psi0)**2 + Kpsi2(psi - psi0)**3 + Kpsi3(psi - psi0)**4"
	puts $fl "!"
	puts $fl "!atom types            psi0    Kpsi1   Kpsi2   Kpsi3"
	array set parameter_array $improperparameters
	foreach n [array names parameter_array] {
		set nn [split $n "-"]
		lassign $nn a1 a2 a3 a4
		lassign $parameter_array($n) k n psi
		puts $fl [format "%5s %5s %5s %5s %6.2f %6.2f %6.2f %6.2f" $a1 $a2 $a3 $a4 $psi $k 0.0 0.0]
	}
}

proc write_pair {fl pairparameters} {
	puts $fl "NONBONDED"
	puts $fl "!" 
	puts $fl "V(Lennard-Jones) = Eps,i,j\[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6\]"
	puts $fl "!"
	puts $fl "!atom  ignored    epsilon      Rmin/2"
	array set parameter_array $pairparameters
	foreach n [array names parameter_array] {
		set nn [split $n "-"]
		lassign $nn a1
		lassign $parameter_array($n) i e R
		puts $fl [format "%5s %6.2f %6.2f %6.2f" $a1 $i $e $R]
	}
}

proc read_parameters {file} {
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	if {[regexp {^NONBONDED} $line]} {
			set pairparameters [read_pairs $infile]
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  		if {[regexp {^BONDS} $line] && ![regexp {^CMAP} $line]} {
  			set bondparameters [read_bonds $infile]
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	 if {[regexp {^ANGLES} $line]} {
  	 	set angleparameters [read_angles $infile]
  	}
  }
 	close $infile
  set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
  	 if {[regexp {^DIHEDRALS} $line]} {
  	 		set dihedralparameters [read_dihedrals $infile]
		}
	}
	close $infile
	set infile [open $file "r"]
	while {[gets $infile line] >= 0} {
		if {[regexp {^IMPROPER} $line]} {
			set improperparameters [read_impropers $infile]
		}
	}
	close $infile 
	return [list $pairparameters $bondparameters $angleparameters $dihedralparameters $improperparameters]
}

proc read_pairs {infile} {
	while {[gets $infile line] >= 0} {
		if {[regexp {^HBOND} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line]} {
			break
		} {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]) [lreplace $linearray 0 0]
    		}
  	}
  }
  return [array get parameters]
}

proc read_bonds {infile} {
	while {[gets $infile line] >= 0} {
		if {[regexp {^HBOND} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line]} {
			break
		} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]-[lindex $linearray 1]) [lreplace $linearray 0 1]
    		}
  	}
  }
 return [array get parameters]
}

proc read_angles {infile} {
	while {[gets $infile line] >= 0} {
		if {[regexp {^HBOND} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line]} {
			break
		} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
	    		set parameters([lindex $linearray 0]-[lindex $linearray 1]-[lindex $linearray 2]) [lreplace $linearray 0 2]
  			}
  	}
  }
  return [array get parameters]
}

proc read_dihedrals {infile} {
	while {[gets $infile line] >= 0} {
			if {[regexp {^HBOND} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line]} {
				break
			} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]-[lindex $linearray 1]-[lindex $linearray 2]-[lindex $linearray 3]) [lreplace $linearray 0 3]
  			}
  	}
  }
  return [array get parameters]
}

proc read_impropers {infile} {

	while {[gets $infile line] >= 0} {
			if {[regexp {^HBOND} $line] || [regexp {^NONBONDED} $line] || [regexp {^BONDS} $line] || [regexp {^ANGLES} $line] || [regexp {^DIHEDRALS} $line] || [regexp {^IMPROPER} $line] || [regexp {^CMAP} $line]} {
				break
			} else {
  	    set linearray [split [lindex [split $line !] 0 ]]
    		set linearray [noblanks $linearray]
    		if {[llength $linearray]} {
    			set parameters([lindex $linearray 0]-[lindex $linearray 1]-[lindex $linearray 2]-[lindex $linearray 3]) [lreplace $linearray 0 3]
  			}
  	}
  }
  return [array get parameters]
}