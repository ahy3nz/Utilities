set frames [molinfo top get numframes]
set start 1
set every 2
set end 301

for {set i $start} {$i < $end} {incr i $every} {
	animate goto $i
    set j [expr ($i - $start) / $every]
    puts $i
    puts $j
	render TachyonInternal foo-$j.tga
}
