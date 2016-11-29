##
## Gui for CG/Hybrid Builder
##
## Authors: Lemmin Thomas && Christophe Bovigny
##
## 
##

## Tell Tcl that we're a package and any dependencies we may have

package provide multiscalegui 0.11

package require multiscaletools 0.11
package require readParTop 0.11
package require writelammps 0.11

namespace eval ::multiscalegui:: {

  # window handle
  variable w                                          
  variable toCGResidueFrame
  variable fromCGFrame                                          
  variable toCGShapeFrame                                          
                                 
  
 

  variable toResCGMenu
  variable toShapeCGMenu
  variable fromMolMenu
  variable fromAAMolMenu
  variable cgpath
  variable toCGoutpdbfile
  variable toCGoutallpdbfile
  variable toCGouttopfile
  variable toCGoutparmfile
  variable fromCGoutpdbfile
  variable revcgfile

 

  #LBM STUFF
  variable pdbpath
  variable psfpath
  variable CG_name "CG_model"
  variable CG_rev "cgrev"
  variable Hybrid_name "Hybrid_model"
  variable Hybrid_rev "hybridrev"
  variable toCGResidueHybrid 0
  variable toCGshapeMass 
  variable drawn 0
  variable drawnlbm 
  variable hybridChoice "0"
  variable CG_model 0 
  #  kAtomsPerDXPoint chosen to make hook.dx generate 15 beads.
  #  taken from the number of points in the file (listed as the last
  # integer in the object 3 line).  24955 / 550 / 3 => 15, which is what
  # we want
  

  variable menuChoice 0
  

}

# -------------------------------------------------------------------------
# should the user be able to pick 'next'
proc ::multiscalegui::toggleNextButton { a b c } {
   variable menuChoice
   variable w
   if { $menuChoice > 0 } {   
     $w.chooseFrame.next configure -state normal
   } else {
     $w.chooseFrame.next configure -state disabled
   }
}

# -------------------------------------------------------------------------

# -----------------------------------------------------------------------------
# LBM Hybrid Box
proc ::multiscalegui::toggleHybridBox { a b c } {
variable toCGResidueHybrid
variable toCGResidueFrame
#variable drawnlbm
	#if { $drawnlbm == 0 } {
  # 	return
  # 		}
	if { $toCGResidueHybrid == 1 } {
	$toCGResidueFrame.hybrid.hybcutoffAA configure -state normal
	$toCGResidueFrame.hybrid.cutoffAA configure -state normal
        $toCGResidueFrame.hybrid.hybcutoffCG configure -state normal
        $toCGResidueFrame.hybrid.cutoffCG configure -state normal
        $toCGResidueFrame.hybrid.hybcenter configure -state normal
        $toCGResidueFrame.hybrid.center configure -state normal
        $toCGResidueFrame.hybrid.preview configure -state normal
	$toCGResidueFrame.hybrid.explaintexthybrid configure -state normal
	} else {
	$toCGResidueFrame.hybrid.hybcutoffAA configure -state disabled
	$toCGResidueFrame.hybrid.cutoffAA configure -state disabled
        $toCGResidueFrame.hybrid.hybcutoffCG configure -state disabled
        $toCGResidueFrame.hybrid.cutoffCG configure -state disabled
        $toCGResidueFrame.hybrid.hybcenter configure -state disabled
        $toCGResidueFrame.hybrid.center configure -state disabled
        $toCGResidueFrame.hybrid.preview configure -state disabled
        $toCGResidueFrame.hybrid.explaintexthybrid configure -state disabled
	}
}

# -------------------------------------------------------------------------

### LBM STUFF
proc ::multiscalegui::calcCutoffValues { } {
#if {[string trim $::multiscalegui::tohybridcutoffAA] == "" || \
#       ! [string is double $::multiscalegui::tohybridcutoffAA] || \
#       $::multiscalegui::tohybridcutoffAA < 1 || [string trim $::multiscalegui::tohybridcutoffCG] == "" || \
#       ! [string is double $::multiscalegui::tohybridcutoffCG] || \
#       $::multiscalegui::tohybridcutoffCG < 1 } {
#          return 0
#  }
  
}
# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

# --LBM-----------------------------------------------------------------------
proc ::multiscalegui::hybridSourceChoice { } {
   variable toCGResidueFrame
   variable hybridChoice
   variable toCGResidueHybrid

      

   if { $hybridChoice == "0" } {
      # molecule
      $toCGResidueFrame.hybrid.hybcutoffAA configure -state disabled
      $toCGResidueFrame.hybrid.cutoffAA configure -state disabled
      $toCGResidueFrame.hybrid.hybcutoffCG configure -state disabled
      $toCGResidueFrame.hybrid.cutoffCG configure -state disabled
      $toCGResidueFrame.hybrid.hybcenter configure -state disabled
      $toCGResidueFrame.hybrid.center configure -state disabled
      $toCGResidueFrame.hybrid.preview configure -state disabled
      $toCGResidueFrame.hybrid.explaintexthybrid configure -state disabled
   
#      $toCGShapeFrame.mol.menu configure -state normal

     set toCGResidueHybrid 0
     
   } elseif { $hybridChoice == "1" } {
  #    # EDM
      $toCGResidueFrame.hybrid.hybcutoffAA configure -state normal
      $toCGResidueFrame.hybrid.cutoffAA configure -state normal
      $toCGResidueFrame.hybrid.hybcutoffCG configure -state normal
      $toCGResidueFrame.hybrid.cutoffCG configure -state normal
      $toCGResidueFrame.hybrid.hybcenter configure -state normal
      $toCGResidueFrame.hybrid.center configure -state normal
      $toCGResidueFrame.hybrid.preview configure -state normal
      $toCGResidueFrame.hybrid.explaintexthybrid configure -state normal
      set toCGResidueHybrid 0
      
   }

}


# -------------------------------------------------------------------------
#
# Create the window and initialize data structures
#
proc ::multiscalegui::multiscalegui {} {
#   puts "very beginning of ::multiscalegui::cggui";
  variable toResCGMenu
  variable toShapeCGMenu
  
  variable w
  variable chooseFrame
  variable chooseCGFrame
  variable toCGResidueFrame
  variable toCGShapeFrame
  variable fromCGFrame



#  ::multiscalegui::init_vars

#  puts "before everything starts\n";

  #trace add variable ::multiscalegui::currentToCGMol write ::multiscalegui::tomolmenuaux
  #trace add variable ::multiscalegui::currentFromCGMol write ::multiscalegui::frommolmenuaux
  #trace add variable ::multiscalegui::currentFromCGAAMol write ::multiscalegui::fromaamolmenuaux

  #trace add variable ::multiscalegui::currentMapCGMol write ::multiscalegui::mapcgmenuaux
  #trace add variable ::multiscalegui::currentMapAAMol write ::multiscalegui::mapaamenuaux
  #trace add variable ::multiscalegui::currentMapRefMol write ::multiscalegui::maprefmenuaux

  #trace add variable ::multiscalegui::ljparamCurrentMol write ::multiscalegui::ljparammenuaux

  trace remove variable ::multiscalegui::menuChoice write ::multiscalegui::toggleNextButton
  trace add variable ::multiscalegui::menuChoice write ::multiscalegui::toggleNextButton

  #trace add variable ::multiscalegui::toCGShapeUseMass write ::multiscalegui::toggleMassBox
  #LBM STUFF
  trace add variable ::multiscalegui::toCGResidueHybrid write ::multiscalegui::toggleHybridBox


  # If already initialized, just turn on
  if { [winfo exists .multiscalegui] } {
    wm deiconify $w
    return
  }

  set w [toplevel ".multiscalegui"]
  wm title $w "Multiscale Builder - Main Menu"

  #Add a menubar
  frame $w.menubar -relief raised -bd 2
  #grid  $w.menubar -padx 1 -column 0 -columnspan 5 -row 0 -sticky ew
  pack $w.menubar -padx 1 -fill x

  menubutton $w.menubar.help -text "About Authors" -underline 0 \
    -menu $w.menubar.help.menu
  $w.menubar.help config -width 12
  pack $w.menubar.help -side right




  ## help menu
  menu $w.menubar.help.menu -tearoff no
  $w.menubar.help.menu add command -label "LBM" \
    -command {tk_messageBox -type ok -title "LBM" \
              -message "This tools was implemented by Lemmin/Bovigny Modeling (LBM) Compagny"}
  $w.menubar.help.menu add command -label "Laboratory of Biomolecular Modeling" \
    -command "vmd_open_url http://lbm.epfl.ch"
  
 


# now, we define a few frames


## ---------------------------------------------------------------------------
## ---------------------- START choose order FRAME ---------------------------
  set chooseFrame [frame $w.chooseFrame]
  set row 0

  #puts "before intro text"
  # intro text
  grid [label $chooseFrame.introText -text "Coarse Graining"] \
     -row $row -column 0 -columnspan 2 -sticky ew
  incr row

  # -----------------------------------------------
  grid [labelframe $chooseFrame.resCg -bd 2 -relief ridge \
            -text "Coarse Grain and Hybrid Plugin : " \
            -padx 1m -pady 1m] -row $row -column 0 -columnspan 2 -sticky nsew
  incr row

  grid [radiobutton $chooseFrame.resCg.r1 -text \
                         "Create a Coarse Grain or Hybrid Model for LAMMPS " \
                   -variable ::multiscalegui::menuChoice -value 1 ] \
        -row 0 -column 0 -sticky ew

 grid [radiobutton $chooseFrame.resCg.r2 -text \
                  "Write Lammps data from pdb/psf" \
                  -variable ::multiscalegui::menuChoice -value 3 ] \
       -row 1 -column 0 -sticky w
  # ---

  #grid [labelframe $chooseFrame.sCg -bd 2 -relief ridge \
  #          -text "Hybrid CG-AA Tools" \
  #          -padx 1m -pady 1m] -row $row -column 0 -columnspan 2 -sticky nsew
  #incr row

  #grid [radiobutton $chooseFrame.sCg.r1 -text \
  #                       "Create CG-AA Model" \
  #                 -variable ::multiscalegui::menuChoice -value 2 ] \
  #      -row 0 -column 0 -sticky w


  grid [button $chooseFrame.next -text "Next ->" -state disabled \
        -command {
           if { $::multiscalegui::menuChoice > 0 } {   
              pack forget $::multiscalegui::chooseFrame   
              if { $::multiscalegui::menuChoice == 1 } {   
                 pack $::multiscalegui::toCGResidueFrame   
                 wm title $::multiscalegui::w "CG Builder - LBM Based"
             # } elseif { $::multiscalegui::menuChoice == 2 } {
                 #pack $::multiscalegui::toCGResidueFrame   
             #    pack $::multiscalegui::toCGShapeFrame   
             #    wm title $::multiscalegui::w "CG Builder - Shape-Based CG"
              } elseif { $::multiscalegui::menuChoice == 3} {
                 pack $::multiscalegui::fromCGFrame   
                 wm title $::multiscalegui::w "Reverse LBM-BASED CG" 
             #  } elseif { $::multiscalegui::menuChoice == 4 } {
             #    pack $::multiscalegui::mapSBFrame   
             #    wm title $::multiscalegui::w "Map Shaped-Based CG To All Atom"
              #} elseif { $::multiscalegui::menuChoice == 5 } {
              #   pack $::multiscalegui::ljParamFrame   
              #   wm title $::multiscalegui::w "Get Lennard-Jones Params for CG model"
              #} elseif { $::multiscalegui::menuChoice == 6 } {
              #   pack $::multiscalegui::baFromSimFrame 
              #   wm title $::multiscalegui::w "Extract Bond/Angle Params from AA Sim"
              #} elseif { $::multiscalegui::menuChoice == 7 } {
              #   pack $::multiscalegui::scaleBAFrame 
              #   wm title $::multiscalegui::w "Scale Bond/Angle Spring Constants"
              } else {
                 tk_messageBox -type ok -message \
                    "Unknown Menu Choice.  Contact Developers" -title "Error!"
              }
           }   
        }] \
        -row $row -column 0 -columnspan 2 -sticky ew
  incr row

## ---------------------- END   choose order FRAME ---------------------------
## ---------------------------------------------------------------------------

#   ::multiscalegui::createToCGResidueFrame {}
#puts "before doing residue TO frame"
## ---------------------------------------------------------------------------
## ---------------------- START to CG via RESIDUE FRAME ---------------------
  set toCGResidueFrame [frame $w.toCGResidueFrame]
  set row 0

  #puts "before intro text"
  # intro text
  grid [label $toCGResidueFrame.introText -text "Coarse Grain / Hybrid Plugin"] \
     -row $row -column 0 -columnspan 2 -sticky w
  incr row
  grid [label $toCGResidueFrame.introTextDesc1 -text \
          "Convert All Atom pdb/psf files to CG/Hybrid pdb/psf and mapping files"] \
     -row $row -column 0 -columnspan 2 -sticky w
  incr row
  grid [label $toCGResidueFrame.introTextDesc2 -text \
          "If you're vegetarian : close this plugin, it could be dangerous for you !!!"] \
     -row $row -column 0 -columnspan 2 -sticky w
  incr row

  
  grid [labelframe $toCGResidueFrame.inp -bd 2 -relief ridge \
            -text "Input" \
            -padx 1m -pady 1m] -row $row -column 0 -columnspan 2 -sticky nsew

  incr row
 
  grid [label $toCGResidueFrame.inp.inpdblabel -text "input PDB"] \
      -row $row -column 0 -sticky w

  grid [entry $toCGResidueFrame.inp.inpdbpath -width 56 \
                                      -textvariable ::multiscalegui::pdbpath ] \
      -row $row -column 1 -sticky ew
  grid [button $toCGResidueFrame.inp.inpdbbutton -text "Browse" \
        -command {
           set tempfile [tk_getOpenFile]
           if {![string equal $tempfile ""]} { set ::multiscalegui::pdbpath $tempfile }
        }] -row $row -column 2 -sticky w

  incr row
	grid [label $toCGResidueFrame.inp.inpsflabel -text "input PSF"] \
      -row $row -column 0 -sticky w

  grid [entry $toCGResidueFrame.inp.inpsfpath -width 56 \
                                      -textvariable ::multiscalegui::psfpath ] \
      -row $row -column 1 -sticky ew
  grid [button $toCGResidueFrame.inp.inpsfbutton -text "Browse" \
        -command {
           set tempfile [tk_getOpenFile]
           if {![string equal $tempfile ""]} { set ::multiscalegui::psfpath $tempfile }
        }] -row $row -column 2 -sticky w
	incr row

  #puts "before cg database"
  # -----------------------------------------------
  # deal with the CG database
  grid [labelframe $toCGResidueFrame.database -bd 2 -relief ridge \
            -text "CG Database" \
            -padx 1m -pady 1m] -row $row -column 0 -columnspan 2 -sticky nsew
  incr row


  set rowDB 0

  grid [label $toCGResidueFrame.database.pLabel -text "Proteins"] \
      -row $rowDB -column 0 -sticky w
  grid [label $toCGResidueFrame.database.pPath \
                   -text "([file join $::env(MULTISCALE) lbm.cgc])"] \
      -row $rowDB -column 1 -columnspan 2 -sticky w
  grid [button $toCGResidueFrame.database.paddbutton -text "Add" \
        -command {
           ::multiscaletools::read_db [file join $::env(MULTISCALE) lbm.cgc]
           $::multiscalegui::numBeadsLabel configure \
                                  -text "[llength $::multiscaletools::convbeads]"
           $::multiscalegui::toCGResidueFrame.database.paddbutton \
                                                    configure -state disabled
           $::multiscalegui::toCGResidueFrame.database.paddbutton \
                                                    configure -text "Added!"
        }] -row $rowDB -column 3 -sticky e
  incr rowDB

  grid [label $toCGResidueFrame.database.wLabel -text "Water"] \
      -row $rowDB -column 0 -sticky w
  grid [label $toCGResidueFrame.database.wPath \
                       -text "([file join $::env(CGTOOLSDIR) water.cgc])"] \
      -row $rowDB -column 1 -columnspan 2 -sticky w
  grid [button $toCGResidueFrame.database.waddbutton -text "Add" \
        -command {
           ::multiscaletools::read_db [file join $::env(CGTOOLSDIR) water.cgc]
           $::multiscalegui::numBeadsLabel configure \
                                  -text "[llength $::multiscaletools::convbeads]"
           $::multiscalegui::toCGResidueFrame.database.waddbutton \
                                                    configure -state disabled
           $::multiscalegui::toCGResidueFrame.database.waddbutton \
                                                    configure -text "Added!"
        }] -row $rowDB -column 3 -sticky e
  incr rowDB

  grid [label $toCGResidueFrame.database.udLabel -text "User Defined"] \
      -row $rowDB -column 0 -sticky w

  grid [entry $toCGResidueFrame.database.udpath -width 46 \
                                      -textvariable ::multiscalegui::cgpath ] \
      -row $rowDB -column 1 -sticky ew
  grid [button $toCGResidueFrame.database.udbutton -text "Browse" \
        -command {
           set tempfile [tk_getOpenFile]
           if {![string equal $tempfile ""]} { set ::multiscalegui::cgpath $tempfile }
        }] -row $rowDB -column 2 -sticky w
  grid [button $toCGResidueFrame.database.udaddbutton -text "Add" \
        -command {
           ::multiscaletools::read_db $::multiscalegui::cgpath
           $::multiscalegui::numBeadsLabel configure \
                                  -text "[llength $::multiscaletools::convbeads]"
        }] -row $rowDB -column 3 -sticky e

  incr rowDB


  grid [label $toCGResidueFrame.database.beadsText \
                              -text "Bead Definitions Currently Loaded:"] \
     -row $rowDB -column 1 -sticky w
  set ::multiscalegui::numBeadsLabel [label $toCGResidueFrame.database.beadsNum \
                                  -text "[llength $::multiscaletools::convbeads]"]
  grid $::multiscalegui::numBeadsLabel -row $rowDB -column 2 -sticky e
  incr rowDB

 

#### LBM LAB HYBRID PART

# deal with the Hybrid System
   grid [labelframe $toCGResidueFrame.hybrid -bd 2 -relief ridge \
            -text "Hybrid Part" \
            -padx 1m -pady 1m] -row $row -column 0 -columnspan 5 -sticky nsew
  incr row

set rowHB 0
  grid [checkbutton $toCGResidueFrame.hybrid.choiceButton -variable \
            [namespace current]::toCGResidueHybrid] -row 0 -column 0 -sticky w 
	

incr rowHB

   grid [label $toCGResidueFrame.hybrid.hybcutoffAA -state disabled -text \
               "Cutoff All Atom :"] -row $rowHB -column 0 -sticky w
   grid [entry $toCGResidueFrame.hybrid.cutoffAA -state disabled -width 7 \
                                      -textvariable ::multiscalegui::tohybridcutoffAA ] \
        -row $rowHB -column 1 -sticky e 
                                      
   

  incr rowHB

 grid [label $toCGResidueFrame.hybrid.hybcutoffCG -state disabled -text \
               "Cutoff Coarse Grain :"] -row $rowHB -column 0 -sticky w
  grid [entry $toCGResidueFrame.hybrid.cutoffCG  -state disabled -width 7 \
               -textvariable ::multiscalegui::tohybridcutoffCG] \
    -row $rowHB -column 1 -sticky e
incr rowHB
 grid [label $toCGResidueFrame.hybrid.explaintexthybrid -state disabled -text "Give a selection to center your hybrid system (ex : resid 5 and backbone)"] \
     -row $rowHB -column 1 -sticky e
incr rowHB
 grid [label $toCGResidueFrame.hybrid.hybcenter -state disabled -text \
               "Center :" ] -row $rowHB -column 0 -sticky w
 grid [entry $toCGResidueFrame.hybrid.center -state disabled -width 46 \
               -textvariable ::multiscalegui::tohybridcenter] \
    -row $rowHB -column 1 -sticky e
    incr rowHB
 #grid [label $toCGResidueFrame.hybrid.void -text "                                      "] -row $rowHB -column 2
 grid [button $toCGResidueFrame.hybrid.preview -state disabled -text "Preview For Hybrid" \
        -command ::multiscalegui::buildHybridExecute]  \
        -row $rowHB -column 1 -sticky e
incr row
 
	
incr rowHB


# Output filenames
  grid [label $toCGResidueFrame.outpdblabel -text "Output PDB: "] \
    -row $row -column 0 -sticky w
  grid [entry $toCGResidueFrame.toCGoutpdbfile -width 30 -textvariable ::multiscalegui::toCGoutpdbfile] \
    -row $row -column 1 -sticky ew
  incr row
  grid [label $toCGResidueFrame.revcglabel -text "Rev CG File: "] -row $row \
               -column 0 -sticky w 
  grid [entry $toCGResidueFrame.revcgfile -width 30 \
           -textvariable ::multiscalegui::revcgfile] -row $row -column 1 -sticky ew 
  incr row
  
  grid [button $toCGResidueFrame.back -text "Back To Previous Screen" \
        -command { \
            pack forget $::multiscalegui::toCGResidueFrame   
            set ::multiscalegui::menuChoice 0
            wm title $::multiscalegui::w "CG Builder - Method Selection"
            pack $::multiscalegui::chooseFrame }] \
        -row $row -column 0 -sticky w
 grid [button $toCGResidueFrame.applyDB -text "Build Model" \
        -command ::multiscalegui::build_model] \
        -row $row -column 1 -sticky e
  incr row


## ---------------------- END    to residue CG rep FRAME -------------------
## ---------------------------------------------------------------------------


## ---------------------------------------------------------------------------
## ---------------------- START fromCGFrame ---------------------

set fromCGFrame [frame $w.fromCGFrame]
  set row 0

  #puts "before intro text"
  # intro text
  grid [label $fromCGFrame.introText -text "Convert pDB/PSF to Lammps data"] \
     -row $row -column 0 -columnspan 2 -sticky w
  incr row
  grid [label $fromCGFrame.introTextDesc1 -text \
          "Creation of data files for LAMMPS, from CoarseGrained files pdb, psf"] \
     -row $row -column 0 -columnspan 2 -sticky w
 incr row
	grid [labelframe $fromCGFrame.inp -bd 2 -relief ridge \
            -text "Input" \
            -padx 1m -pady 1m] -row $row -column 0 -columnspan 2 -sticky nsew

  incr row
 
  grid [label $fromCGFrame.inp.inpdblabel -text "input PDB"] \
      -row $row -column 0 -sticky w

  grid [entry $fromCGFrame.inp.inpdbpath -width 56 \
                                      -textvariable ::multiscalegui::pdbpath ] \
      -row $row -column 1 -sticky ew
  grid [button $fromCGFrame.inp.inpdbbutton -text "Browse" \
        -command {
           set tempfile [tk_getOpenFile]
           if {![string equal $tempfile ""]} { set ::multiscalegui::pdbpath $tempfile }
        }] -row $row -column 2 -sticky w

  incr row
	grid [label $fromCGFrame.inp.inpsflabel -text "input PSF"] \
      -row $row -column 0 -sticky w

  grid [entry $fromCGFrame.inp.inpsfpath -width 56 \
                                      -textvariable ::multiscalegui::psfpath ] \
      -row $row -column 1 -sticky ew
  grid [button $fromCGFrame.inp.inpsfbutton -text "Browse" \
        -command {
           set tempfile [tk_getOpenFile]
           if {![string equal $tempfile ""]} { set ::multiscalegui::psfpath $tempfile }
        }] -row $row -column 2 -sticky w
	incr row
# Output filenames
  grid [label $fromCGFrame.outpdblabel -text "Output Lammps: "] \
    -row $row -column 0 -sticky w
  grid [entry $fromCGFrame.toCGoutpdbfile -width 30 -textvariable ::multiscalegui::toCGoutpdbfile] \
    -row $row -column 1 -sticky ew
  incr row
 
 
 

  grid [button $fromCGFrame.back -text "Back To Previous Screen" \
        -command { \
            pack forget $::multiscalegui::fromCGFrame   
            set ::multiscalegui::menuChoice 0
            wm title $::multiscalegui::w "CG Builder - Method Selection"
            pack $::multiscalegui::chooseFrame }] \
        -row $row -column 0 -sticky w
 grid [button $fromCGFrame.applyDB -text "Build DataFile" \
        -command ::multiscalegui::build_data] \
        -row $row -column 1 -sticky e
  incr row 


## ---------------------- END    to residue CG rep FRAME -------------------
## ---------------------------------------------------------------------------


  pack $chooseFrame 

  # this trace lets the plugin determine when you have loaded a molecule
  # in VMD
  trace add variable ::vmd_initialize_structure write \
    ::multiscalegui::vmd_init_struct_trace

}

# -------------------------------------------------------------------------







# -------------------------------------------------------------------------
proc ::multiscalegui::addStatusLine {strText} {
   variable toCGShapeFrame
   ::multiscalegui::addGenericStatusLine $toCGShapeFrame $strText
}

# -------------------------------------------------------------------------
proc ::multiscalegui::addGenericStatusLine {frameName strText} {
   $frameName.statFrame.statBox configure -state normal
   $frameName.statFrame.statBox insert end "$strText\n"
   $frameName.statFrame.statBox see end
   update idletask
   $frameName.statFrame.statBox configure -state disabled
}


# -------------------------------------------------------------------------


# -------------------------------------------------------------------------

# -------------------------------------------------------------------------
proc multiscalegui_tk {} {
#   puts "before invoking the constructor";
  ::multiscalegui::multiscalegui
#   puts "after invoking the constructor";
	# READ PARAMETERS and TOPOLOGY files  
  ::readParTop::read_topology $::env(MULTISCALE)/top_lbm.rtf
  ::readParTop::read_topology $::env(MULTISCALE)/top_K51.rtf
  ::readParTop::read_topology $::env(MULTISCALE)/top_amber2charmm.inp
  ::readParTop::read_parameters $::env(MULTISCALE)/par_lbm.prm
  #AA Amber parameters
  #::readParTop::read_topology $::env(MULTISCALE)/top_amber2charmm.inp
  ::readParTop::read_parameters $::env(MULTISCALE)/par_amber_lammps.inp
  # parameters and topology of ligand K49
  ::readParTop::read_parameters $::env(MULTISCALE)/par_K51.inp
  set ::multiscaletools::DEBUG 0
  return $::multiscalegui::w
}

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------

# -------------------------------------------------------------------------


# -------------------------------------------------------------------------
proc ::multiscalegui::viewAA { molid atom1 cutoff} {
    
    set sel [atomselect $molid "same residue as within $cutoff of ($atom1)"]
    #mol delrep 0 top
    mol selection "[$sel text]"
    mol color Type 
    mol representation CPK
    mol material Opaque
    mol addrep top
    #mol selupdate 0 top 0
    #mol colupdate 0 top 0
    #mol showrep top 0 0

    
}

proc ::multiscalegui::viewCG { molid atom1 cutoff } {
    set sel [atomselect $molid "same residue as (not within $cutoff of ($atom1))"]
    #mol delrep 0 top
    mol selection "[$sel text]"
    mol color ColorID 9
    mol material Transparent
    mol representation Beads 0.9 8
    mol addrep top
    #mol selupdate 1 top 0
    #mol colupdate 1 top 0
    #mol showrep top 1 0
}

proc ::multiscalegui::viewHybrid {molid atom1 cutoffAA cutoffCG } {
set pure_AA [atomselect $molid "same residue as within $cutoffAA of ($atom1)"]
set pure_CG [atomselect $molid "same residue as (not within $cutoffCG of ($atom1))"]
set hybrid [atomselect $molid "not (index [$pure_CG get index] [$pure_AA get index])"]
    #mol selection "[$pure_CG text]"
    
    #mol representation licorice
    #mol addrep top
    #mol selection "[$pure_AA text]"
    #mol color type
    #mol representation CPK
    #mol addrep top
    mol selection "[$hybrid text]"
    #mol color ColorID 7
    mol color Type    
    mol material Steel
    mol representation licorice
    mol addrep top
}

proc ::multiscalegui::view { molidAA atomAA cutoffAA cutoffCG } {

viewAA $molidAA $atomAA $cutoffAA


viewCG $molidAA $atomAA $cutoffCG

viewHybrid $molidAA $atomAA $cutoffAA $cutoffCG

}
proc ::multiscalegui::test { } {
mol delrep 0 top
mol delrep 0 top
mol delrep 0 top
set molidAA 0
set atomAA "resid 5 and backbone"
set cutoffAA 5
set cutoffCG 10
::multiscalegui::view $molidAA $atomAA $cutoffAA $cutoffCG

}

proc ::multiscalegui::build_model { } {
         # Test if output or not
         set outpdb $::multiscalegui::toCGoutpdbfile
         if { [string trim $::multiscalegui::toCGoutpdbfile] == "" } {
			 tk_messageBox -message "Put an output \n Gerce ignoble" \
                   -type ok -title "Error!"
         return
         }
	
	 # Test if rcg or not
         set outrevcg $::multiscalegui::revcgfile
         if { [string trim $::multiscalegui::revcgfile] == "" } {
			 tk_messageBox -message "Put a REVCG => " \
                   -type ok -title "Error!"
         return
         }

	if { $::multiscalegui::toCGResidueHybrid == 1 } {
		file delete $::multiscalegui::Hybrid_name.pdb $::multiscalegui::Hybrid_name.psf $::multiscalegui::Hybrid_rev
		set ::multiscalegui::Hybrid_rev $::multiscalegui::revcgfile
		set ::multiscalegui::Hybrid_name $::multiscalegui::toCGoutpdbfile
		::multiscalegui::build_hybrid
                
		set hyb [mol new $::multiscalegui::Hybrid_name.psf]
		mol addfile $::multiscalegui::Hybrid_name.pdb
		mol addfile $::multiscalegui::Hybrid_name.psf
		set id [[atomselect top "all"] molid]
		puts $id
		set all [atomselect top "all"]
		set triplet [::multiscaletools::write_dipole $id [topo -molid $id getanglelist] $::readParTop::patches]
		::multiscaletools::set_connectivitytype $id $::readParTop::parameters $triplet
		::writelammps::writelammpsdata $id $::multiscalegui::Hybrid_name.data multiscale $all $::readParTop::parameters $::multiscaletools::dipoles $::multiscaletools::ellipsoids
		::multiscaletools::set_group $id
		::writelammps::writelammpsinput $::multiscalegui::Hybrid_name.data $::multiscaletools::sc_dipole $::multiscaletools::standard

	} else {
		file delete $::multiscalegui::CG_name.pdb $::multiscalegui::CG_name.psf $::multiscalegui::CG_rev
		set ::multiscalegui::CG_rev $::multiscalegui::revcgfile
		set ::multiscalegui::CG_name $::multiscalegui::toCGoutpdbfile
		::multiscalegui::build_CG
        set CG [mol new $::multiscalegui::CG_name.psf]
		mol addfile $::multiscalegui::CG_name.pdb
		mol addfile $::multiscalegui::CG_name.psf
		set idCG [[atomselect top "all"] molid]
		set allCG [atomselect top "all"]
		[atomselect top "all"] set occupancy 3
		set triplet [::multiscaletools::write_dipole $idCG [topo -molid $idCG getanglelist] $::readParTop::patches]
		::multiscaletools::set_connectivitytype $idCG $::readParTop::parameters $triplet		
		::writelammps::writelammpsdata $idCG $::multiscalegui::CG_name.data multiscale $allCG $::readParTop::parameters $::multiscaletools::dipoles $::multiscaletools::ellipsoids
		::multiscaletools::set_group $idCG
        ::writelammps::writelammpsinput $::multiscalegui::CG_name.data $::multiscaletools::sc_dipole $::multiscaletools::standard
	}
}

###################################################################################
#													Coarse grained procs																		#
#																																									#
###################################################################################

proc ::multiscalegui::build_CG { } {
		
		set AA [mol new $::multiscalegui::psfpath]
		mol addfile $::multiscalegui::pdbpath
		#::multiscaletools::read_db $CG_topology
		::multiscaletools::apply_database $AA $::multiscalegui::CG_name.pdb $::multiscalegui::CG_rev 1
		
		::multiscaletools::psfgen $::multiscalegui::CG_name.pdb $AA $::multiscalegui::CG_name.psf $::multiscalegui::CG_rev $::readParTop::topology
		mol delete $AA
}

proc ::multiscalegui::load_AA_CG {psfAA pdbAA pdbCG psfCG} {
		set AA [mol new $psfAA]
		mol addfile $pdbAA
		set CG [mol new $psfCG]
		mol addfile $pdbCG
		return [list $AA $CG]
}


proc ::multiscalegui::display_hybrid {psfhybrid pdbhybrid} {
 mol delete all
 set current [mol new $psfhybrid]
 mol addfile $pdbhybrid
 mol delrep 0 $current
 #pure AA
 mol selection "occupancy 0"
 mol color ColorID 1
 mol material Opaque
 mol representation CPK 1.0 0.3 8 6
 mol addrep $current
 #hybrid AA
 mol selection "occupancy 1"
 mol color ColorID 8
 mol material Opaque
 mol representation Licorice
 mol addrep $current
 #hybrid CG
 mol selection "occupancy 2"
 mol color ColorID 8
 mol material Transparent
 mol representation CPK 3.0 1.5 8.0 8.0
 mol addrep $current 
 #pure CG
 mol selection "occupancy 3"
 mol color ColorID 19
 mol material Opaque
 mol representation CPK 3.0 1.5 8.0 8.0
 mol addrep $current 
}


###################################################################################
#													Hybrid procs																						#
#																																									#
###################################################################################

proc ::multiscalegui::buildHybridExecute { } {

if { $::multiscalegui::CG_model == 0 } {
	::multiscalegui::build_CG
	set ::multiscalegui::CG_model 1
}
set cutAA -1
   
      # need to have a positive value for the CA cutoff
      if {[string trim $::multiscalegui::tohybridcutoffAA] == "" || \
          ! [string is double $::multiscalegui::tohybridcutoffAA] || \
          $::multiscalegui::tohybridcutoffAA < 0 } {
         tk_messageBox -message "Cutoff AA must be an double or a pokemon." \
                   -type ok -title "Error!"
         return
      }
      set cutAA $::multiscalegui::tohybridcutoffAA
   

set cutCG -1
   
      # need to have a positive value for the CG cutoff
      if {[string trim $::multiscalegui::tohybridcutoffCG] == "" || \
          ! [string is double $::multiscalegui::tohybridcutoffCG] || \
          $::multiscalegui::tohybridcutoffCG < 0 } {
         tk_messageBox -message "Cutoff CG must be an double or a pokemon." \
                   -type ok -title "Error!"
         return
      }
     
      set cutCG $::multiscalegui::tohybridcutoffCG
   	if { $cutCG < $cutAA } {
   			tk_messageBox -message "Cut off CG has to be greater than cut off all atom!" \
                   -type ok -title "Error!"
         return
    }
		set center -1 
		if { [string trim $::multiscalegui::tohybridcenter] == "" } {
			 tk_messageBox -message "Center selection is not defined. Encule!" \
                   -type ok -title "Error!"
         return
    }
    set center $::multiscalegui::tohybridcenter

	if { [::multiscalegui::build_hybrid]==0 } {
		::multiscalegui::display_hybrid $::multiscalegui::Hybrid_name.psf $::multiscalegui::Hybrid_name.pdb
	}
      
}


proc ::multiscalegui::build_hybrid { } {
	set l [load_AA_CG $::multiscalegui::psfpath $::multiscalegui::pdbpath $::multiscalegui::CG_name.pdb \
              $::multiscalegui::CG_name.psf]
	lassign $l idAA idCG
	set out [::multiscaletools::select_multiscale $::multiscalegui::CG_rev $::multiscalegui::tohybridcenter \
                $::multiscalegui::tohybridcutoffAA $::multiscalegui::tohybridcutoffCG \
		$idAA $idCG $::multiscalegui::Hybrid_rev $::multiscalegui::Hybrid_name $::readParTop::topology]
      if { $out == -1 } {
	tk_messageBox -message "Unvalid selection Do not bullshit with me, man!" \
        -type ok -title "Error!"
       return -1
       } 
       if { $out == -2 } {
       tk_messageBox -message "No atoms in hybrid selection. Change cutoffs." \
       -type ok -title "Error!"
       return -1
       }
       
  return 0
}

proc ::multiscalegui::build_data { } {
		 
		# Test if pdb or not
         set outrevcg $::multiscalegui::pdbpath
         if { [string trim $::multiscalegui::pdbpath] == "" } {
			 tk_messageBox -message "YOU HAVE TO PUT A PDB FILE IN, OTHERWISE THOMAS LEMMIN COMES AND BEAT YOU => http://www.facebook.com/#!/profile.php?id=573598702  " \
                   -type ok -title "!!!!!!!INPUT PDB ERROR!!!!!!!"
         return
         }
	  set outpsf $::multiscalegui::psfpath
         if { [string trim $::multiscalegui::psfpath] == "" } {
			 tk_messageBox -message "YOU HAVE TO PUT A PSF FILE IN, OTHERWISE THOMAS LEMMIN COMES AND BEAT YOU => http://www.facebook.com/#!/profile.php?id=573598702 " \
                   -type ok -title "!!!!!!!INPUT PSF ERROR!!!!!!"
         return
         }

		  # Test if output or not
       		  set outpdb $::multiscalegui::toCGoutpdbfile
        	  if { [string trim $::multiscalegui::toCGoutpdbfile] == "" } {
		 tk_messageBox -message "PUT AN OUTPUT FILENAME, PLEASE. WE ARE SURE THAT YOU ARE  \n Une Gerce Ignoble" \
                  -type ok -title "!!!!!!!OUTPUT ERROR!!!!!!"
         		return
         	}

		
		set CG [mol new $::multiscalegui::psfpath]
		mol addfile $::multiscalegui::pdbpath
		mol addfile $::multiscalegui::psfpath
		set id [[atomselect top "all"] molid]
		puts $id
		set all [atomselect top "all"]
		set triplet [::multiscaletools::write_dipole $id [topo -molid $id getanglelist] $::readParTop::patches]
		::multiscaletools::set_connectivitytype $id $::readParTop::parameters $triplet
		::writelammps::writelammpsdata $id $::multiscalegui::toCGoutpdbfile.data multiscale $all $::readParTop::parameters $::multiscaletools::dipoles $::multiscaletools::ellipsoids

		::multiscaletools::set_group $id
		::writelammps::writelammpsinput $::multiscalegui::toCGoutpdbfile.data $::multiscaletools::sc_dipole $::multiscaletools::standard


}




