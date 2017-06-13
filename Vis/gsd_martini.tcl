# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# Log file 'C:/Users/ayang/Documents/Research/ZeroTiltDSPC-100/tclcode', created by user ayang
#mol addrep 0
mol modselect 0 0 name P4 or name BP4
mol modstyle 0 0 CPK 1.00 0.300 10.00 10.00
mol modmaterial 0 0 Transparent
mol color Name
mol representation Points 1.000000
mol selection name P4 or name BP4
#mol material BrushedMetal
mol addrep 0
mol modselect 1 0 not name P4 and not name BP4
mol modmaterial 1 0 Opaque
mol modstyle 1 0 Licorice 0.300000 10.000000 10.000000
mol color Name
mol representation Licorice 0.300000 10.000000 10.000000
mol selection not name P4 and not name BP4
mol material Opaque
mol modrep 1 0
mol smoothrep 0 1 5
mol smoothrep 0 0 5
display resetview
rotate x to 90
axes location off
display projection Orthographic
pbc box -center origin
pbc unwrap
scale to 0.2
animate goto 5
# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# end of log file.
