package require multiscaletools

proc multi_coarse_grained { wd cgfile name n } {
	::multiscaletools::read_db $cgfile 
	cd $wd
	for {set i 0} {$i < $n} {incr i} {
		set AA [mol new "$name\_$i.psf"]
		mol addfile "$name\_$i.pdb"
		::multiscaletools::apply_database $AA "CG_$name\_$i.pdb" "rev_$i"
		::multiscaletools::psfgen "CG\_$name\_$i.pdb" $AA "CG_$name\_$i.psf" "rev_$i"
		set CG [mol new "CG_$name\_$i.psf"]
    	mol addfile "CG_$name\_$i.pdb"
		set idCG [[atomselect $CG "all"] molid]
		set allCG [atomselect $CG "all"]
		::multiscaletools::writelammpsdata $idCG "CG_$name\_$i.data" CG $allCG
		mol delete $AA
    }
}