# mutliscaletools is a vmd package that allows conversion to and from coarse grain 
#  representations. See the docs for details on file formats.
# set lstColor [array get color]                ;#(convert array to list)
# array set color $lstColor                     ;#(convert list to array)


package provide multiscaletools 0.11
if {[catch {package require topotools 1.0} ver]} {
   vmdcon -error "$ver. This script requires at least TopoTools v1.0. Exiting..."
   quit
}
if {[catch {package require pbctools 2.3} ver]} {
   vmdcon -error "$ver. This script requires at least pbctools v2.3. Exiting..."
   quit
}


namespace eval ::multiscaletools:: {
  # List of currently recognized conversions
  variable convbeads [list]
  variable exceptions [list]
  variable idxarray
  variable startposarray
  variable aagroup [list]
  variable cggroup [list]
  variable hybgroup [list]
  variable sc_dipole [list]
  variable standard [list]
  variable parameters [list]
  variable DEBUG 0
  variable dipoles [list]
  variable ellipsoids [list]
  namespace export read_db apply_reversal apply_database
}

proc ::multiscaletools::resetconvs {} {
  #Proc to reset the list of conversions currently being held
  variable convbeads
  set convbeads [list]
}

proc ::multiscaletools::noblanks {mylist} {
  set newlist [list]
  foreach elem $mylist {
    if {$elem != ""} {
      lappend newlist $elem
    }
  }

  return $newlist
}

proc ::multiscaletools::read_db {file} {
  #Read all cg conversions in file and add them to the convbeads list
  variable convbeads
  variable exceptions
  variable DEBUG
  set infile [open $file "r"]
  
  while {[gets $infile line] >= 0} {
  	if {[regexp {^EXBEGIN} $line]} {
  		set exceptions [join [lappend exceptions [join [read_exception $infile]]]]
  	} elseif {[regexp {^CGBEGIN} $line]} {
      set newbead [read_bead $infile]
      lappend convbeads $newbead
      if {$DEBUG==1} {
      	puts $newbead
      }
    }
  }

  close $infile
}

proc ::multiscaletools::read_exception {fstream} {
	set	myexception [list]
	while {[gets $fstream line] && ![regexp {^EXEND} $line]} {
    	if {[regexp {^\#} $line]} { 
     		continue 
    	}
    	set linearray [split $line]
    	set linearray [noblanks $linearray]
    	lappend myexception [lindex $linearray 0]
   }
   return $myexception
}

proc ::multiscaletools::read_bead {fstream} {
  #Given a file stream currently starting a new bead, read the bead's components
  # and return the new atom list
  # The bead "object" is simply a list of beads and atoms, as defined in make_atom
  # The first entry in this list is the "bead"; the others are the target atoms

  set mybead [list]
  
  while {[gets $fstream line] && ![regexp {^CGEND} $line]} {
    if {[regexp {^\!} $line]} { 
     continue 
    }

    #split the line up into fields and make a new atom
    set linearray [split $line]
    set linearray [noblanks $linearray]
    set newatom [make_atom [lindex $linearray 0] [lindex $linearray 1] [lindex $linearray 2] [lindex $linearray 3]]
    lappend mybead $newatom
  }

  return $mybead
}

proc ::multiscaletools::make_atom {resname atomname resoffset dipole} {
  #Create a new cg bead/atom with atomname in element 0,resname in element 1, and a resid
  # offset in element 2

  set newatom [list]

#  set newatom(Resname) $resname
#  set newatom(Atomname) $atomname
#  set newatom(Offset) $resoffset
  
  lappend newatom $atomname
  lappend newatom $resname
  lappend newatom $resoffset
   lappend newatom $dipole
  return $newatom
}

proc ::multiscaletools::apply_bead {cgbead molid revcgfile DIPOLE} {
  #Applies the conversion specified in CGBEAD to the molecule MOLID
  # This means going though the molecule, matching everything that
  # corresponds to the first element of cgbead, and then building
  # each of those beads in turn

  set beadname [lindex [lindex $cgbead 0] 0]
  set beadres [lindex [lindex $cgbead 0] 1]
  set beadoff [lindex [lindex $cgbead 0] 2]
  set cgbead [lreplace $cgbead 0 0]
  set headname [lindex [lindex $cgbead 0] 0]
  set headres [lindex $[lindex $cgbead 0] 1]
  set headoff [lindex $[lindex $cgbead 0] 2]
  set headdipole [lindex $[lindex $cgbead 0] 3]
  set cgbead [lreplace $cgbead 0 0]
  variable exceptions


#puts "DEBUG: Looking to make cgbead $beadname $beadres $beadoff with head atom $headname $headres $headoff"
  #Find all of the atoms matching the head definition
  if {$headres == "*"} {
    set headbeads [atomselect $molid "name $headname and not resname $exceptions"]
  } else {
    set headbeads [atomselect $molid "name $headname and resname $headres"]
  }
  $headbeads set occupancy 1
#  $headbeads set name $beadname
  if {$beadres != "*"} {
     $headbeads set resname $beadres
  }

  #Make three arrays of the qualifying characteristics of subordinate beads
  set names [list]
  set resnames [list]
  set offsets [list]
  set weight_list [list]
  lappend weight_list $headname
  lappend weight_list $headoff
  
  foreach subatom $cgbead {
    lappend names [lindex $subatom 0]
    lappend resnames [lindex $subatom 1]
    lappend offsets [lindex $subatom 2]
    lappend weight_list [lindex $subatom 0]
    lappend weight_list [lindex $subatom 2]
  }

  array set weight_array $weight_list
  
  #Now, for each head atom we've found, adjust its position according to 
  # where its children are
  foreach index [$headbeads get index] resid [$headbeads get resid] segid [$headbeads get segid] {
    set headatom [atomselect $molid "index $index"]
    if {[join [$headatom get occupancy] ] < 0} {
      $headatom delete
      continue
    }
    #puts "Applying bead with head $index"
    $headatom set resid $resid
    $headatom set beta [expr 0.01*[$headatom get index]]
    set resid [expr $resid]
    set fullbeadsel "occupancy >= 0 and segid $segid and ( index $index "	
    set selmove "occupancy >= 0 and segid $segid and ( index $index "
    foreach name $names resname $resnames offset $offsets {
   	set selchild "\{ name $name and resid $resid and resname "
      
      if {$resname == "*"} {
         set selchild "$selchild [$headatom get resname]\} "
      } else {
         set selchild "$selchild $resname \}"
      }
      set fullbeadsel "$fullbeadsel or $selchild"
      if {$offset == 1} {
     		set selmove "$selmove or $selchild"
      }
    }

    set fullbeadsel " $fullbeadsel ) "
    set selmove " $selmove )"
    set mybeadsel [atomselect $molid "$fullbeadsel"]
    set sel_mass [[atomselect $molid "$selmove"] get mass]
    set sel_names [[atomselect $molid "$selmove"] get name]
    set masses [list]
    set weight [list]
    foreach m $sel_mass n $sel_names {
    	lappend masses [expr $m * $weight_array($n)]
    	lappend weight $weight_array($n)
    }
        
    $headatom moveto [measure center [atomselect $molid "$selmove"] weight $masses]
	$headatom set name $beadname
    puts $revcgfile "[$headatom get index] [$headatom get resname] [$headatom get name] [$headatom get resid] [$headatom get segid] [$mybeadsel get index] :$weight"
    $mybeadsel set occupancy -1
    $headatom set occupancy 1
    $headatom delete
    $mybeadsel delete
  }

  $headbeads delete
    
}

# Each CG bead in CGId molecule has a respective domain
# in the all-atom reference file (refId). Atoms from each
# domain of refId can be corresponded to atoms of AAId
# (for two atoms to be declared identical, segname, resid, 
# and name should match);
# a bead that represents the center of mass of a given refId domain
# should be moved to the center of mass of the corresponding AAId domain.
proc ::multiscaletools::mapCGMolecule { statusProc AAId CGId refId outPDB} {

   set AAsel [atomselect $AAId all]
   set CGsel [atomselect $CGId all]
   set refsel [atomselect $refId all]

   set Naa [$AAsel num]
   set Ncg [$CGsel num]
   set Nref [$refsel num]

# Fill out the array of coordinates for CG beads.
   set cg_pos_0 [$CGsel get {x y z}]
   for {set i 1} {$i <= $Ncg} {incr i} {
      set cg_pos($i) {0.0 0.0 0.0}
      set cg_norm($i) 0.0
   }

# Lists of segnames, resnames, and names for AAId molecule.
   set AAsegname [$AAsel get segname]
   set AAresid [$AAsel get resid]
   set AAname [$AAsel get name]
# Lists of coordinates and masses for AAId molecule.
   set AApos [$AAsel get {x y z}]
   set AAmass [$AAsel get mass]
   
# For AAId, set all betas to 0 initially.
   $AAsel set beta 0
   for {set i 0} {$i < $Naa} {incr i} {
#      puts "i = $i"

      # Beta field for each atom in refId contains the serial number
      # of the CG bead corresponding to the domain containing the atom.
      # Find the atom in refId which is the same as the one currently chosen from AAId.
      set A [atomselect $refId "name [lindex $AAname $i] and resid [lindex $AAresid $i] and segname [lindex $AAsegname $i]"]
      # If refId does not contain this atom, ignore it.
      # If refId contains one or more atoms satisfying the name, resid, and segname
      # of the chosen atom from AAId, use coordinates and mass of the first
      # of such atoms from refId (in principle, there should be only one such atom).
      if {[$A num] > 0} {
         # Get the serial of the bead to whose domain this refId atom belongs.
         set i_CG [expr int([lindex [$A get beta] 0])]
         # The bead will be moved to the center of mass of the corresponding 
         # AAId domain; add the coordinates of the current AAId atom to the
         # variable that will be used to obtain the coordinates of the bead.
         set cg_pos($i_CG) [vecadd $cg_pos($i_CG) [vecscale [lindex $AAmass $i] [lindex $AApos $i]]]
         set cg_norm($i_CG) [expr $cg_norm($i_CG) + [expr [lindex $AAmass $i]]]
      }
      $A delete
   }

# Update coordinates for CG beads.
   for {set i 1} {$i <= $Ncg} {incr i} {
      set cg_pos($i) [vecscale [expr 1.0/$cg_norm($i)] $cg_pos($i)]
      set A [atomselect $CGId "serial $i"]
      $A set x [lindex $cg_pos($i) 0]
      $A set y [lindex $cg_pos($i) 1]
      $A set z [lindex $cg_pos($i) 2]
      $A delete
   }
   set A [atomselect $CGId all]
#Write the output PDB.
   $A writepdb $outPDB
   $A delete

# Move beads back to old coordinates.
   for {set i 0} {$i < $Ncg} {incr i} {
      set A [atomselect $CGId "index $i"]
      $A set x [lindex [lindex $cg_pos_0 $i] 0]
      $A set y [lindex [lindex $cg_pos_0 $i] 1]
      $A set z [lindex [lindex $cg_pos_0 $i] 2]
      $A delete
   }

# clean up the memory that we grabbed at the first
   $AAsel delete
   $CGsel delete
   $refsel delete
#   mol delete $CGId

# Load the created CG structure into VMD.
   mol load pdb $outPDB
   mol modstyle 0 top VDW 5.000000 10.000000
}

proc ::multiscaletools::apply_database {molid outputfile revcgfile DIPOLE} {
  #Applies the contents of the current convbeads database to the
  # selected molecule, and writes the result to OUTPUTFILE
	set tmp_rev "rev.tmp"
  variable convbeads
 
  
  if {[llength $convbeads] == 0} {
     puts "Multiscale Tool Error) No bead definitions loaded.  Can't apply."
     return
  }
  #Open file for reverse coarse graining information
  # format of file is resname beadname segid index1 index2 index3...
  # where the first 3 fields come from the CG bead, and the indices are corresponding
  # all atom indices
  set rcgout [open $tmp_rev "w"]
  #Beads which should be kept and written are tagged with occupancy 1
  # All other atoms should have occupancy 0
  set allsel [atomselect $molid all]
  set oldocc [$allsel get occupancy]
  set oldxyz [$allsel get {x y z}]
  $allsel set occupancy 0

  #Loop through the conversion database and do each bead type
  foreach cgbead $convbeads {
    apply_bead $cgbead $molid $rcgout $DIPOLE
  }

  set writesel [atomselect $molid "occupancy 1"]
  $writesel writepdb $outputfile
  $writesel delete
	
  $allsel set occupancy $oldocc
  $allsel set {x y z} $oldxyz
  $allsel delete
  close $rcgout
  update_revcg $tmp_rev $outputfile $revcgfile
  file delete $tmp_rev
}

proc ::multiscaletools::update_revcg {tmp_rev outputfile revcgfile} {
	set CGpdb [mol new $outputfile]
	set CG [atomselect $CGpdb all]
	set CG_indices [$CG get index]
	set CG_beta [$CG get beta]
	foreach beta $CG_beta index $CG_indices {
		set a([expr round($beta*100.00)]) $index
	}
	set outfile [open  $revcgfile "w"]
	set infile [open $tmp_rev "r"]
	while {[gets $infile line] >= 0} {
		set linelist [split $line]
    set linelist [noblanks $linelist]
    set index [lindex $linelist 0]
    set resname [lindex $linelist 1]
    set beadname [lindex $linelist 2]
    set resid [lindex $linelist 3]
    set segid [lindex $linelist 4]
    set indices [lreplace $linelist 0 4]
    puts $outfile "$a($index) $resname $beadname $resid $segid $indices"
    }
   close $infile
   close $outfile
   mol delete $CGpdb
}

proc ::multiscaletools::read_revcg {revcgfile} {
	set infile [open $revcgfile "r"]
	while {[gets $infile line] >= 0} {
    set linelist [split $line]
    set linelist [noblanks $linelist]
    set index [lindex $linelist 0]
    set resname [lindex $linelist 1]
    set beadname [lindex $linelist 2]
    set resid [lindex $linelist 3]
    set segid [lindex $linelist 4]
    set indices [lindex [split [lreplace $linelist 0 4] ":" ] 0]
    foreach i $indices {
    	set tmp($i) $index 
    }
 	}
 	set mapping [array get tmp]
 	close $infile
 	return $mapping
}

proc ::multiscaletools::find_atom {a lookuplist resid segid} {
	array set lookup $lookuplist
	if { [string compare + [string index $a 0]] == 0 } {
			set key "$segid,[expr $resid+1],[string replace $a 0 0]"
			if {[info exists lookup($key)]} {
				return [list [string replace $a 0 0] $lookup($key)]
			}
	} elseif { [string compare - [string index $a 0]] == 0 } {
			set key "$segid,[expr $resid-1],[string replace $a 0 0]"
			#puts "- $key"
			if {[info exists lookup($key)]} {
				return [list [string replace $a 0 0] $lookup($key)]
			}
	} else {
			set key "$segid,$resid,$a"
			if {[info exists lookup($key)]} {
				return [list $a $lookup($key)]
			}
	}
	return [list -1 -1]
}


proc ::multiscaletools::select_multiscale {revcgfile center_AA cutoff_AA cutoff_CG idAA idCG revhybrid outHybrid topology} {

		set mapping_list [read_revcg $revcgfile]
		[atomselect $idAA "all"] set occupancy 5
		[atomselect $idCG "all"] set occupancy 5
		set AA [atomselect $idAA "all"]
		save_indices $AA
		set CG [atomselect $idCG "all"]
		save_indices $CG
		if { [catch { set pure_AA [atomselect $idAA "same residue as within $cutoff_AA of ($center_AA)"]}] } {
			return -1
		}
		if { [catch { set pure_CG [atomselect $idAA "same residue as (not within $cutoff_CG of ($center_AA))"]}] } 			{
			return -1
		}
		set hybrid [atomselect $idAA "not (index [$pure_CG get index] or index [$pure_AA get index])"]
		if {[$hybrid num] == 0} {
			return -2
		}
		$pure_AA set occupancy 0
		$hybrid set occupancy 1
		
		#find corresponding CG beads
		set mapping_list [read_revcg $revcgfile]
		
		array set mapping $mapping_list
		set indexlist_CG [$pure_CG get index]
		set index_pure_CG [list]
		foreach index $indexlist_CG {
			lappend index_pure_CG $mapping($index)
		}
		#find corresponding hyrbid in CG
		set indexlist_hybrid [$hybrid get index]
		set index_hybrid_CG [list]
		foreach index $indexlist_hybrid {
			lappend index_hybrid_CG $mapping($index)
		}
		
		[atomselect $idCG "index $index_pure_CG"] set occupancy 3
		[atomselect $idCG "index $index_hybrid_CG"] set occupancy 2
		
		######################################################################################
		set AA [atomselect $idAA "occupancy 0 1"]
		set CG [atomselect $idCG "occupancy 2 3"]
		#change segname and chain
		set seg [$CG get segname]
		set chain [$CG get chain]
		set nseg [list]
		set nchain [list]
		foreach s $seg c $chain {
			lappend nseg "C$s"
			lappend nchain "C$c"
		}
		$CG set segname $nseg
		$CG set chain $nchain
		set AAfile "all_atom"
		set CGfile "CG_atom"
		$AA writepdb $AAfile.pdb
		$AA writepsf $AAfile.psf
		$CG writepdb $CGfile.pdb
		$CG writepsf $CGfile.psf
		#####################################################################################
		resetpsf
		readpsf  $AAfile.psf
		coordpdb $AAfile.pdb
		readpsf  $CGfile.psf
		coordpdb $CGfile.pdb
		writepsf $outHybrid.psf
		writepdb $outHybrid.pdb
		save_all $AAfile.pdb $CGfile.pdb tmp.pdb indices.pdb
		write_hybrid_revcg indices.pdb $outHybrid.pdb $mapping_list $revhybrid
		#clean up
		file delete tmp.pdb indices.pdb $AAfile.pdb $AAfile.psf $CGfile.pdb $CGfile.psf
		resetpsf
		# need to set mass and charge correctly. problem with readpsf....
		lassign $topology atomtopology chargetopology masses patches impropertopology
		set i [mol new $outHybrid.psf]
		mol addfile $outHybrid.pdb
		mol addfile $outHybrid.psf
		set sel [atomselect $i all]
		set_charge $sel $chargetopology
		set_mass $sel $masses
		$sel writepsf $outHybrid.psf
		$sel writepdb $outHybrid.pdb
		#mol delete $i
		return 0
}

proc ::multiscaletools::save_indices {sel} {
		set indices [list]
		foreach r [$sel get index] {
			lappend indices [expr 0.01*$r]
		}
		$sel set beta $indices
}

proc ::multiscaletools::fcat { CGfile AAfile CG_AAfile } {
		set fCG [open $CGfile "r"]
		set fAA [open $AAfile "r"]
		set fout [open $CG_AAfile "a"] 
		fcopy $fAA $fout
		fcopy $fCG $fout
		close $fAA
		close $fCG
		close $fout
}

proc ::multiscaletools::save_all { AAfile CGfile infile outfile } {
		fcat $CGfile $AAfile $infile
		set fin [open $infile "r"]
		set fout [open $outfile "w"]
		while {[gets $fin line] >= 0} {
    	if {[regexp {^ATOM} $line]} {
					puts $fout $line
			}
		}
		close $fin
		close $fout
}

proc ::multiscaletools::write_hybrid_revcg { indicesfile pdbfile mapping_list revhybrid } {
		array set mapping $mapping_list
		set currentpdb [mol new $indicesfile]
		set pdbindices [[atomselect $currentpdb all] get beta]
		set pdboccupancy [[atomselect $currentpdb all] get occupancy]
		mol delete $currentpdb
		set indices [list]
		foreach i $pdbindices {
			lappend indices [expr round(100*$i)]
		}
		
		set currentpdb [mol new $pdbfile]
		set hybridpdb [atomselect $currentpdb "all"]
		$hybridpdb set occupancy $pdboccupancy
		set indiceshybrid [$hybridpdb get index]
		set a [list]
		foreach ih $indiceshybrid i $indices {
			lappend a $ih
			lappend a $i
		}
		array set hybridmapping $a 
		set AA [atomselect $currentpdb "occupancy 0 1"]
		set AAindices [$AA get index]
		set CG [atomselect $currentpdb "occupancy 2 3"]
		set CGindices [$CG get index]
		set hybridAA [atomselect $currentpdb "occupancy 1"]
		set hybridCG [atomselect $currentpdb "occupancy 2"]
		set hAAbeta [list]
		set hCGbeta [list]
		set i 1
		foreach a [$hybridCG get index] {
				set mappingCG($hybridmapping($a)) $i
				lappend hCGbeta $i
				incr i
		}
		foreach a [$hybridAA get index] {
				lappend hAAbeta $mappingCG($mapping($hybridmapping($a)))
		}
		$hybridCG set beta $hCGbeta
		$hybridAA set beta $hAAbeta
		$hybridpdb writepdb $pdbfile
		
		#reverse mapping
		set fmapping [open $revhybrid "w"]
		foreach	i $AAindices {
			puts $fmapping "$i $hybridmapping($i)"
		}
		array set rev_mapping [reverse_array $mapping_list]
		foreach i $CGindices {
			puts $fmapping "$i $rev_mapping($hybridmapping($i))"
		}
		close $fmapping
}

proc ::multiscaletools::reverse_array { mapping_list } {
	set value_unique [list] 
	foreach {key value} $mapping_list {
		set new($value) [list]
	}
	set tmp [list]
	foreach {key value} $mapping_list {
		set tmp $new($value)
		lappend tmp $key
		set new($value) $tmp
	}
	return [array get new]
}

proc ::multiscaletools::set_group { id } {
	variable aagroup
	variable cggroup
	variable hybgroup
	set aagroup [[atomselect $id "occupancy 0 1"] get serial]
	set cggroup [[atomselect $id "occupancy 3"] get serial]
	set hybgroup [[atomselect $id "occupancy 1"] get serial]
}

##################################### Topology #########################################
##################################### Topology #########################################
##################################### Topology #########################################

proc ::multiscaletools::psfgen {pdbCG AA psfCG revCG topology} {
		lassign $topology atomtopology chargetopology masses patches impropertopology
		if {[llength $topology]!=0} {
			set bondlist_AA [topo -molid $AA -sel protein getbondlist]
			set CG [mol new $pdbCG]
			set sel [atomselect $CG "all"]
			set bondlist_CG [guessbonds $revCG $bondlist_AA]
			topo -molid $CG setbondlist $bondlist_CG
			
			set bondlist_CG_nonunique [$sel getbonds]
			set anglelist_CG [guessangles  $sel $bondlist_CG_nonunique]
			topo -molid $CG setanglelist $anglelist_CG
			set dihedrallist_CG [guessdihedrals $sel $bondlist_CG]
			topo -molid $CG setdihedrallist $dihedrallist_CG
			set improperlist_CG [set_improperlist $sel $impropertopology]
			topo -molid $CG setimproperlist $improperlist_CG
			set_type $sel $atomtopology
			
			set_charge $sel $chargetopology
			set_mass $sel $masses
			#not very clean... but needed to set correct charge
			add_cap $CG $patches $atomtopology $chargetopology
			
			$sel writepsf $psfCG
			mol delete $CG
		} else {
			puts "In psfgen: Please provide topology file"
			return -1
		}
		return 1
		
}
proc ::multiscaletools::set_type {sel atomtopology} {
    
    if {[llength $atomtopology]==0} {
		puts "No atom types have been defined. Please provide topology file"
		return -1
	}
	
    set name [$sel get name]
	set resname [$sel get resname]
	set types [list]
	array set atom_array $atomtopology
	foreach n $name r $resname {
		lappend types $atom_array($r-$n)
	}
	$sel set type $types
}
proc ::multiscaletools::set_charge {sel chargetopology} {
	
    if {[llength $chargetopology]==0} {
		puts "No charges have been defined. Please provide topology file"
		return -1
	}
	set name [$sel get name]
	set resname [$sel get resname]
	set charges [list]
	array set charge_array $chargetopology
	foreach n $name r $resname {
		lappend charges $charge_array($r-$n)
	}
	$sel set charge $charges
}

proc ::multiscaletools::set_mass {sel masses} {
    if {[llength $masses]==0} {
		puts "No masses have been defined. Please provide topology file"
		return -1
	}

	  set type [$sel get type]
	  set mass [list]
	  array set mass_array $masses
	  foreach t $type {
	  	lappend mass $mass_array($t)
	  }
	  $sel set mass $mass
}

proc ::multiscaletools::add_patch {CG bondlist_CG patches atomtopology chargetopology} {
	add_cap $CG $patches $atomtopology $chargetopology
	add_disulfide $CG $bondlist_CG $patches
}

proc ::multiscaletools::add_cap {CG patches atomtopology chargetopology} {
	array set patches_array $patches
	if {[info exists patches_array(NTER)] && [info exists patches_array(CTER)] } {
		array set atom_array $atomtopology
		array set charge_array $chargetopology
		set sel [atomselect $CG all]
		set segname [lsort -unique [$sel get segname]]
		foreach s $segname {
			set sel [atomselect $CG "segname $s and name $patches_array(NTER) $patches_array(CTER)"]
			set type [$sel get type]
			set charge [$sel get charge]
			set resname [$sel get resname]
			set type [lreplace $type 0 0 $atom_array([lindex $resname 0]N-$patches_array(NTER))]
			set type [lreplace $type end end $atom_array([lindex $resname end]C-$patches_array(CTER))]
			set charge [lreplace $charge 0 0 $charge_array([lindex $resname 0]N-$patches_array(NTER))]
			set charge [lreplace $charge end end $charge_array([lindex $resname end]C-$patches_array(CTER))]
			$sel set type $type
			$sel set charge $charge
		}
	}
}

proc ::multiscaletools::add_disulfide {CG bondlist_CG patches} {
	array set patches_array $patches
	if {[info exists patches_array(DISU)]} {
		set sel [atomselect $CG "name [lindex [lindex $patches_array(DISU) 0] 0]"]
		set type [$sel get type]
		set indices [$sel get index]
		foreach	i $indices {
			foreach ii $indices {
				if {[lsearch $bondlist_CG [list $i $ii]]==""} {
					set type [lreplace $type [lsearch $indices $i] [lsearch $indices $i] [lindex [lindex $patches_array(DISU) 0] 1] ]
					set type [lreplace $type [lsearch $indices $ii] [lsearch $indices $ii] [lindex [lindex $patches_array(DISU) 0] 1] ]
				}
			}
		}	
		$sel set type $type
	}
}

proc ::multiscaletools::guessbonds {revcgfile bondlist_AA} {
	set mapping_list [read_revcg $revcgfile]
	array set mapping $mapping_list

	set bondlist_CG [list]
	foreach bond $bondlist_AA {
		lassign $bond atom1 atom2 
		if { [info exists mapping($atom1)] && [info exists mapping($atom2)] } {
			if { $atom1 != $atom2 && $mapping($atom1) != $mapping($atom2)} {
					lappend bondlist_CG [list $mapping($atom1) $mapping($atom2)]
    		}
    	}
	}
	#puts "DEBUG: guessbonds finished"
	set bondlist_CG [lsort -unique $bondlist_CG]
	return $bondlist_CG
} 

proc ::multiscaletools::guessangles {sel bonddata} {
    set mol [$sel molid]
    set atomtypes [$sel get type]
    set atomindex [$sel list] 
    set newanglelist {}
    foreach bonds $bonddata aidx $atomindex atyp $atomtypes {
        set nbnd [llength $bonds]
        for {set i 0} {$i < $nbnd-1} {incr i} {
            for {set j [expr {$i+1}]} {$j < $nbnd} {incr j} {
                set b1idx [lindex $bonds $i]
                set idx [lsearch -sorted -integer $atomindex $b1idx]
                set b1typ [lindex $atomtypes $idx]
                set b2idx [lindex $bonds $j]
                set idx [lsearch -sorted -integer $atomindex $b2idx]
                set b2typ [lindex $atomtypes $idx]
                if { ([string compare $b1typ $b2typ] > 0) } {
                    set t1 $b1typ; set b1typ $b2typ; set b2typ $t1
                    set t2 $b1idx; set b1idx $b2idx; set b2idx $t2 
                }
                
                #puts "DEBUG: ($i, $j), ($bonds), $aidx, $atyp"
                # append only angles that are full contained in $sel
                if {([lsearch -sorted -integer $atomindex $b1idx] >= 0)          \
                        && ([lsearch -sorted -integer $atomindex $aidx] >= 0)   \
                        && ([lsearch -sorted -integer $atomindex $b2idx] >= 0) } {
                        if {$b1idx < $b2idx} {
	                        set type [join [list $b1typ $atyp $b2typ] "-"]
    	                	lappend newanglelist [list $type $b1idx $aidx $b2idx]
    	                } else {
    	                	set type [join [list $b2typ $atyp $b1typ] "-"]
    	                	lappend newanglelist [list $type $b2idx $aidx $b1idx]
    	                }
                }
            }
        }
    }
    
   return $newanglelist
}

proc ::multiscaletools::uniqueangle {anglelist} {   
	#remove all duplicate angle, e.g. {1 2 3} {2 3 4}
	set newanglelist [list]
	foreach a $anglelist {
		lassign $a a1 a2 a3
		if {[lsearch -all $newanglelist [list $a1 $a2 *]]=="" && [lsearch -all $newanglelist [list * $a2 $a3]]=="" &&\
			[lsearch -all $newanglelist [list $a3 $a2 *]]=="" && [lsearch -all $newanglelist [list * $a2 $a1]]==""} {
			lappend newanglelist $a
		}
	}
	return $newanglelist
}

proc ::multiscaletools::guessdihedrals {sel bondlist} {

    set mol [$sel molid]
    set atomtypes [$sel get type]
    set atomindex [$sel list]
    set newdihedrallist {}
    set bonddata [$sel getbonds]
   

    # a topological dihedral is defined by a bond and atoms
    # bound to it that are not the bond itself
    foreach bond $bondlist {
        lassign $bond b1 b2 
        set b1idx [lsearch -sorted -integer $atomindex $b1]
        set b1typ [lindex $atomtypes $b1idx]
        set b2idx [lsearch -sorted -integer $atomindex $b2]
        set b2typ [lindex $atomtypes $b2idx]
        foreach o1 [lindex $bonddata $b1idx] {
            foreach o2 [lindex $bonddata $b2idx] {
                if {($o1 == $b1) || ($o2 == $b1) || ($o1 == $b2) || ($o2 == $b2)} {
                    continue
                }
                set o1idx [lsearch -sorted -integer $atomindex $o1]
                set o1typ [lindex $atomtypes $o1idx]
                set o2idx [lsearch -sorted -integer $atomindex $o2]
                set o2typ [lindex $atomtypes $o2idx]
                if { ([string compare $b1typ $b2typ] > 0) \
                 || ( [string equal $b1typ $b2typ] 
                      && [string compare $o1typ $o2typ] > 0 ) } {
                    set type [join [list $o2typ $b2typ $b1typ $o1typ] "-"]
                    lappend newdihedrallist [list $type $o2 $b2 $b1 $o1]
                } else {
                    set type [join [list $o1typ $b1typ $b2typ $o2typ] "-"]
                    lappend newdihedrallist [list $type $o1 $b1 $b2 $o2]
                }
            }
        }
    }
   return $newdihedrallist
}

proc ::multiscaletools::uniquedihedral {dihedrallist} {   
	#remove all duplicate dihedral, e.g. {1 2 3 4} {1 2 3 5}
	set newdihedrallist [list]
	foreach d $dihedrallist {
		lassign $d d1 d2 d3 d4
		if {[lsearch -all $newanglelist [list $d1 $d2 $d3 *]]=="" && [lsearch -all $newanglelist [list * $d2 $d3 $d4]]=="" &&\
			[lsearch -all $newanglelist [list $d4 $d3 $d2 *]]=="" && [lsearch -all $newanglelist [list * $d3 $d2 $d1]]==""} {
			lappend newdihedrallist $d
		}
	}
	return $newdihedrallist
}

proc ::multiscaletools::set_improperlist {sel impropertopology} {
	if {[llength $impropertopology]==0} {
		puts "No impropers angles have been defined. Please provide topology file"
		return -1
	}
	array set improper_array $impropertopology
	set residue [$sel get residue]
	set resname [$sel get resname]
	set segid [$sel get segid]
	set resid [$sel get resid]
	set atomtypes [$sel get name]
 	set atomindex [$sel list]
 	foreach s $segid r $residue ri $resid rn $resname at $atomtypes ai $atomindex {
 		set lookup_array($s,$ri,$at) $ai
 		set residue_array($r) [list $s $ri $rn]
 	}
	set  improperlist [list]
	foreach r [lsort -unique $residue] {
		lassign $residue_array($r) segid resid resname
		if { [info exists improper_array($resname)] } {
  			set impropers $improper_array($resname)	
  		} else {
  			set impropers [list]
  		}
  	foreach i $impropers {
  		set a1 [lindex $i 0]
			lassign [find_atom $a1 [array get lookup_array] $resid $segid] a1 a1i
  		if {$a1i ==-1} {
  			continue
  		}
  		set a2 [lindex $i 1]
  		lassign [find_atom $a2 [array get lookup_array] $resid $segid] a2 a2i
  		if {$a2i ==-1} {
  			continue
  		}
  		set a3 [lindex $i 2]
  		lassign [find_atom $a3 [array get lookup_array] $resid $segid] a3 a3i
  		if {$a3i ==-1} {
  			continue
  		}
  		set a4 [lindex $i 3]
  		lassign [find_atom $a4 [array get lookup_array] $resid $segid] a4 a4i
  		if {$a4i ==-1} {
  			continue
  		}
  		lappend improperlist [list "$a1-$a2-$a3-$a4" $a1i $a2i $a3i $a4i]
  	}
  }
	return $improperlist
}




######################Set [bond|angle|dihedral|improper] type###########################
######################Set [bond|angle|dihedral|improper] type###########################
######################Set [bond|angle|dihedral|improper] type###########################

proc ::multiscaletools::set_connectivitytype {id parameters triplet} {
	lassign $parameters pairparameters dipoleparameters ellipsoidparameters bondparameters angleparameters improperparameters dihedralparameters
	variable dipoles [list]
	variable ellipsoids [list]
	set atomtype [[atomselect $id all] get type]
	set atomindex [[atomselect $id all] get index]
	set atomlist [list]
	foreach an $atomtype as $atomindex {
		lappend atomlist $as
		lappend atomlist $an
	} 
	set bondlist [topo -molid $id getbondlist type]
	set bondlist [set_bondtype $bondlist $atomlist $bondparameters]
	topo -molid $id setbondlist type $bondlist
	
	set anglelist [topo -molid $id getanglelist]
	set anglelist [set_angletype $anglelist $atomlist $angleparameters]
	topo -molid $id setanglelist $anglelist
	
	set dihedrallist [topo -molid $id getdihedrallist]
	set dihedrallist [set_dihedraltype $dihedrallist $atomlist $dihedralparameters]
	topo -molid $id setdihedrallist $dihedrallist
	
	set improperlist [topo -molid $id getimproperlist]
	set improperlist [set_impropertype $improperlist $atomlist $improperparameters]
	topo -molid $id setimproperlist $improperlist
	
	set dipoles [set_dipole $atomlist $dipoleparameters $triplet]
	set ellipsoids [set_ellipsoid $atomlist $ellipsoidparameters]
}

proc ::multiscaletools::set_bondtype {bondlist atomlist bondparameters} {
	variable DEBUG
	array set atom_array $atomlist
	if {![llength $bondparameters]} {
			puts "Bond parameters is empty\n"
			return -1
	} else {
		array set bonds_array $bondparameters
		set newbond [list]
		foreach b $bondlist {
			set key1 [list "$atom_array([lindex $b 0])-$atom_array([lindex $b 1])"]
			set key2 [list "$atom_array([lindex $b 1])-$atom_array([lindex $b 0])"]
			if {[info exists bonds_array($key1)]} {
				lappend newbond [list [lindex $b 0] [lindex $b 1] $key1]
			} elseif {[info exists bonds_array($key2)]} {
				lappend newbond [list [lindex $b 1] [lindex $b 0] $key2]
			} else {
				if {$DEBUG==1} {
					puts "WARNING. No paramters were found for bond $key1. The bond has been deleted."
				}
			}
		}
		set newbond [lsort -unique $newbond]
		return $newbond
	}
}

proc ::multiscaletools::set_angletype {anglelist atomlist angleparameters} {
	variable DEBUG
	set cleared 0
	array set atom_array $atomlist
	if {![llength $angleparameters]} {
			puts "Angle parameters is empty\n"
			return -1
	} else {
		array set angle_array $angleparameters
		set newangle [list]
		foreach b $anglelist {
			set key1 [list "$atom_array([lindex $b 1])-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])"]
			set key2 [list "X-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])"]			
			set key3 [list "$atom_array([lindex $b 1])-$atom_array([lindex $b 2])-X"]
			set key4 [list "$atom_array([lindex $b 3])-$atom_array([lindex $b 2])-$atom_array([lindex $b 1])"]
			set key5 [list "X-$atom_array([lindex $b 2])-$atom_array([lindex $b 1])"]			
			set key6 [list "$atom_array([lindex $b 3])-$atom_array([lindex $b 2])-X"]
			if {[info exists angle_array($key1)]} {
				lappend newangle [list $key1 [lindex $b 1] [lindex $b 2] [lindex $b 3]]
			} elseif {[info exists angle_array($key2)]} {
				lappend newangle [list $key2 [lindex $b 1] [lindex $b 2] [lindex $b 3]]
			} elseif {[info exists angle_array($key3)]} {
				lappend newangle [list $key3 [lindex $b 1] [lindex $b 2] [lindex $b 3]]
			} elseif {[info exists angle_array($key4)]} {
				lappend newangle [list $key4 [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} elseif {[info exists angle_array($key5)]} {
				lappend newangle [list $key5 [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} elseif {[info exists angle_array($key6)]} {
				lappend newangle [list $key6 [lindex $b 3] [lindex $b 2] [lindex $b 1]]	
			} else {
				set cleared 1
				if {$DEBUG==1} {
					puts "WARNING. No paramters were found for angle $key1. The angle has been deleted."
				}
			}
		}	
	}
	set newangle [lsort -unique $newangle]
	return $newangle
}

proc ::multiscaletools::set_dihedraltype {dihedrallist atomlist dihedralparameters} {
	variable DEBUG
	array set atom_array $atomlist
	if {![llength $dihedralparameters]} {
			puts "Dihedral parameters is empty\n"
			return -1
	} else {
		array set dihedral_array $dihedralparameters
		set newdihedral [list]
		foreach b $dihedrallist {
			set key1 [list "$atom_array([lindex $b 1])-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])-$atom_array([lindex $b 4])"]
			set key2 [list "X-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])-X"]
			set key3 [list "$atom_array([lindex $b 1])-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])-X"]
			set key4 [list "X-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])-$atom_array([lindex $b 4])"]
			set key5 [list "$atom_array([lindex $b 4])-$atom_array([lindex $b 3])-$atom_array([lindex $b 2])-$atom_array([lindex $b 1])"]
			set key6 [list "X-$atom_array([lindex $b 3])-$atom_array([lindex $b 2])-X"]
			set key7 [list "$atom_array([lindex $b 4])-$atom_array([lindex $b 3])-$atom_array([lindex $b 2])-X"]
			set key8 [list "X-$atom_array([lindex $b 3])-$atom_array([lindex $b 2])-$atom_array([lindex $b 1])"]
			if {[info exists dihedral_array($key1)]} {
				lappend newdihedral [list $key1 [lindex $b 1] [lindex $b 2] [lindex $b 3] [lindex $b 4]]
			} elseif {[info exists dihedral_array($key2)]} {
				lappend newdihedral [list $key2 [lindex $b 1] [lindex $b 2] [lindex $b 3] [lindex $b 4]]
			} elseif {[info exists dihedral_array($key3)]} {
				lappend newdihedral [list $key3 [lindex $b 1] [lindex $b 2] [lindex $b 3] [lindex $b 4]]
			} elseif {[info exists dihedral_array($key4)]} {
				lappend newdihedral [list $key4 [lindex $b 1] [lindex $b 2] [lindex $b 3] [lindex $b 4]]
			} elseif {[info exists dihedral_array($key5)]} {
				lappend newdihedral [list $key5 [lindex $b 4] [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} elseif {[info exists dihedral_array($key6)]} {
				lappend newdihedral [list $key6 [lindex $b 4] [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} elseif {[info exists dihedral_array($key7)]} {
				lappend newdihedral [list $key7 [lindex $b 4] [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} elseif {[info exists dihedral_array($key8)]} {
				lappend newdihedral [list $key8 [lindex $b 4] [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} else {
				if {$DEBUG==1} {
					puts "WARNING. No paramters were found for dihedral $key1. The dihedral has been deleted."
				}
			}
		}
	}
	set newdihedral [lsort -unique $newdihedral]
	return $newdihedral
}

proc ::multiscaletools::set_impropertype {improperlist atomlist improperparameters} {
	variable DEBUG
	array set atom_array $atomlist
	if {![llength $improperparameters]} {
			puts "Improper parameters is empty\n"
			return -1
	} else {
		array set improper_array $improperparameters
		set newimproper [list]
		foreach b $improperlist {
			set key1 [list "$atom_array([lindex $b 1])-$atom_array([lindex $b 2])-$atom_array([lindex $b 3])-$atom_array([lindex $b 4])"]
			set key2 [list "$atom_array([lindex $b 1])-X-X-$atom_array([lindex $b 4])"]
			set key3 [list "X-$atom_array([lindex $b 2])-X-$atom_array([lindex $b 4])"]
			if {[info exists improper_array($key1)]} {
				lappend newimproper [list $key1 [lindex $b 1] [lindex $b 2] [lindex $b 3] [lindex $b 4]]
			} elseif {[info exists improper_array($key2)]} {
				lappend newimproper [list $key2 [lindex $b 4] [lindex $b 3] [lindex $b 2] [lindex $b 1]]
			} elseif {[info exists improper_array($key3)]} {
				lappend newimproper [list $key3 [lindex $b 1] [lindex $b 2] [lindex $b 3] [lindex $b 4]]
			} else {
				if {$DEBUG==1} {
					puts "WARNING. No paramters were found for improper $key1. The improper has been deleted."
				}
			}
		}
	}
	set newimproper [lsort -unique $newimproper]
	return $newimproper
}

proc ::multiscaletools::set_dipole {atomlist dipoleparameters triplet} {
	variable sc_dipole
	variable standard
	#Extract central column triplet
	set BB_dipole []
	foreach row $triplet { lappend BB_dipole [lindex $row 1] }
	set sc_dipole [list]
	set standard [list]
	if {![llength $dipoleparameters]} {
			puts "dipole parameters is empty. Initializing to zeros\n"
			foreach {n a} $atomlist {
				lappend dipoles {0.0 0.0 0.0}
				lappend standard $a
			}
	} else {
		array set dipole_array $dipoleparameters
		set dipoles [list]
		foreach {n a} $atomlist {
			if {[info exists dipole_array($a)]} {
				lappend dipoles $dipole_array($a)
				lappend sc_dipole $a
			} elseif {[lsearch $BB_dipole [expr $n+1]] >= 0} {
				lappend dipoles {0.75 0.0 0.0}
				lappend standard $a
			} else {
				lappend dipoles {0.0 0.0 0.0}
				lappend standard $a
			}
		}
	}
	set standard [lsort -unique $standard]
	set sc_dipole [lsort -unique $sc_dipole]
	return $dipoles
}

proc ::multiscaletools::set_ellipsoid {atomlist ellipsoidparameters} {
	#array set atom_array $atomlist
	set ellipsoids [list]
	#set atoms [array names atom_array]
	if {![llength $ellipsoidparameters]} {
			puts "ellipsoid parameters is empty. Initializing to NULL\n"
	} else {
		array set ellipsoid_array $ellipsoidparameters
		set ellispoids [list]
		set atomid 0
		foreach {n a} $atomlist {
			incr atomid
			if {[info exists ellipsoid_array($a)]} {
				lappend ellipsoids [concat $atomid $ellipsoid_array($a)]
			} 
		}
	}
	return $ellipsoids
}

############################ backbone dipole #################################
############################ backbone dipole #################################
############################ backbone dipole #################################

proc ::multiscaletools::write_dipole {CG anglelist_CG patches} {
	set triplet [list]
	array set patches_array $patches
	if {[info exists patches_array(DIPL)]} {
		set sel [atomselect $CG "name [lindex $patches_array(DIPL) 0]"]
		set indices [$sel get index]
		foreach a $anglelist_CG {
			lassign $a name a1 a2 a3 
			if { [lsearch $indices $a1]!=-1 && [lsearch $indices $a2]!=-1 && [lsearch $indices $a3]!=-1 } {
				if {$a1<$a3} {
					lappend triplet "[expr $a1+1] [expr $a2+1] [expr $a3+1]"
				} else {
					lappend triplet "[expr $a3+1] [expr $a2+1] [expr $a1+1]"
				}
			}
		}
	}
	set fl [open triplet.dat w]
	puts $fl "\# triplet list of the backbone dipoles"
	puts $fl ""
	puts $fl [join [lsort -integer -index 0 $triplet] "\n"]
	close $fl
	return $triplet
}

