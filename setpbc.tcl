#!/usr/bin/env tclsh

set beta 90.0
set pi 3.14159265359
set fp [open [lindex $argv 0]]
set sp [open pbcs.tcl "w"]
set count 0

foreach line [split [read -nonewline $fp] "\n"] {
  if {[lindex $line 1] != "Step"} {
    puts $sp "pbc set {[lindex $line 2] [lindex $line 6] [lindex $line 10] 90.0 $beta 90.0} -first $count -last $count"
    set count [expr $count+1]
  }
}
puts $sp "pbc box"
close $fp
close $sp
