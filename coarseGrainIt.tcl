package require multiscaletools
package require readParTop
package require writelammps
package require psfgen

proc coarseGrainIt {cgfile psfR pdbR psfL pdbL } {
	merge_structure $psfR $pdbR $psfL $pdbL R L
	::readParTop::read_topology $::env(MULTISCALE)/top_lbm.rtf
  	::readParTop::read_parameters $::env(MULTISCALE)/par_lbm.prm
  	::multiscaletools::read_db $cgfile
  	set ::multiscaletools::DEBUG 0
  	set AA [mol new complex.psf]
	mol addfile complex.pdb
	mol addfile complex.psf
	::multiscaletools::apply_database $AA CG_complex.pdb out.log 0
	::multiscaletools::psfgen CG_complex.pdb $AA CG_complex.psf out.log $::readParTop::topology
	mol delete $AA
	set CG [mol new CG_complex.psf]
	mol addfile CG_complex.pdb
	mol addfile CG_complex.psf
	set idCG [[atomselect top "all"] molid]
	set allCG [atomselect top "all"]
	::multiscaletools::set_connectivitytype $idCG $::readParTop::parameters
	::multiscaletools::write_dipole $idCG [topo -molid $idCG getanglelist] $::readParTop::patches
	::writelammps::writelammpsdata $idCG CG_lammps.data CG $allCG $::readParTop::parameters $::multiscaletools::dipoles $::multiscaletools::quaternions
}

proc merge_structure {psfR pdbR psfL pdbL chainR chainL} {
	changeChainSeg $psfR $pdbR $chainR
	changeChainSeg $psfL $pdbL $chainL	
	resetpsf
	readpsf $psfR
	coordpdb $pdbR
	readpsf $psfL
	coordpdb $pdbL
	writepsf complex.psf
	writepdb complex.pdb
}

proc changeChainSeg {psf pdb chain} {
	
	set R [mol new $psf]
	mol addfile $pdb 
	set r [atomselect $R all]
	$r set chain $chain
	set ss [list]
	set s [$r get segname]
	set ps [lindex $s 0]
	set i 1
	foreach S $s {
		if {$S!=$ps} {
			incr i
		}
		lappend ss $chain$i
		set ps $S
	}
	$r set segname $ss
	$r writepsf $psf
	$r writepdb $pdb
	mol delete $R
}