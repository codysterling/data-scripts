set type "mapi"
set x "I"
set x2 "i"
set delta 0.01
set pbc 1
set rmax 10
set first 4000
set last -1

proc writerdf {rdf name} {
  set rdfout [open "rdf-40ps-$name.rdf" w]
  foreach x [lindex $rdf 0] y [lindex $rdf 1] {puts $rdfout "$x $y"}
  close $rdfout
  set rdfout [open "rdf-40ps-$name.int" w]
  foreach x [lindex $rdf 0] y [lindex $rdf 2] {puts $rdfout "$x $y"}
  close $rdfout
}

# all atoms
set vals [measure gofr [atomselect top "all"] [atomselect top "all"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals all-$type

# N-H
set vals [measure gofr [atomselect top "name N"] [atomselect top "name H and within 1.5 of name N"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals nh-$type

# H-X
set vals [measure gofr [atomselect top "name H and within 1.5 of name N"] [atomselect top "name $x"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals h$x2-$type

# N-X
set vals [measure gofr [atomselect top "name N"] [atomselect top "name $x"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals n$x2-$type

# N-C
set vals [measure gofr [atomselect top "name N"] [atomselect top "name C"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals nc-$type

# N-Pb
set vals [measure gofr [atomselect top "name N"] [atomselect top "name Pb"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals npb-$type

# H-Pb
set vals [measure gofr [atomselect top "name H and within 1.5 of name N"] [atomselect top "name Pb"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals hpb-$type

# Pb-X
set vals [measure gofr [atomselect top "name Pb"] [atomselect top "name $x"] delta $delta rmax $rmax usepbc $pbc selupdate 1 first $first last $last]
writerdf $vals pb$x2-$type

