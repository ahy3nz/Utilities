# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# Log file 'C:/Users/ayang/Documents/Research/ZeroTiltDSPC-100/tclcode', created by user ayang
mol delrep 0 0
set rep 0
mol addrep 0
mol modselect 0 0 water
mol modstyle 0 0 Points 1.000000
mol modmaterial 0 0 Transparent
mol color Name
mol representation Points 1.000000
mol selection water
#mol material BrushedMetal
mol addrep 0
mol modselect 1 0 not water and not hydrogen
mol modmaterial 1 0 Opaque
mol modstyle 1 0 Licorice 0.500000 15.000000 12.000000
mol color Name
mol representation Licorice 0.500000 15.000000 12.000000
mol selection not water
mol material Opaque
mol modrep 1 0
mol smoothrep 0 1 2
mol smoothrep 0 0 2
display resetview
rotate x by -90
axes location off
display projection Orthographic
pbc box
scale to 0.03
animate goto 1
# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# end of log file.
