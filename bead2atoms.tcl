proc bead2atoms {file} {
	set infile [open $file r]
	while {[gets $infile line] >= 0} {
		set linelist [split $line]
    	set linelist [noblanks $linelist]
    	set index [lindex $linelist 0]
		set beads_array($index) [lreplace $linelist 0 4]
	}
	close $infile
	set beads [array get beads_array]
	writeIdFile $beads
}

proc writeIdFile {beads} {
	array set beads_array $beads
	set fl [open mapping.dat w]	
	set length [llength [array names beads_array]]
	
	for {set i 0} {$i<$length} {incr i} {
		set atoms [split $beads_array($i) ":"]
		set latoms [llength [lindex $atoms 0]]
		puts $fl [format "%6d" $latoms]
		foreach w [lindex $atoms 1] {
			puts -nonewline $fl [format "%6.1f" $w]
		}
		puts $fl "\r"
		foreach w [lindex $atoms 0] {
			puts -nonewline $fl [format "%6d" [expr $w + 1]]
		}
		puts $fl "\r"
	}
	close $fl

}