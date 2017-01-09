# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# Log file 'C:/Users/ayang/Documents/Research/ZeroTiltDSPC-100/tclcode', created by user ayang
mol addrep 0
mol modselect 0 0 water
mol modstyle 0 0 Points 1.000000
mol modmaterial 0 0 Transparent
mol color Name
mol representation Points 1.000000
mol selection water
#mol material BrushedMetal
mol addrep 0
mol modselect 1 0 not water
mol modmaterial 1 0 Opaque
mol modstyle 1 0 Licorice 0.500000 15.000000 12.000000
mol color Name
mol representation Licorice 0.500000 15.000000 12.000000
mol selection not water
mol material Opaque
mol modrep 1 0
mol smoothrep 0 1 5
mol smoothrep 0 0 5
display resetview
rotate x to 90
axes location off
display projection Orthographic
pbc box
scale to 0.03
animate goto 5
# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# end of log file.
