# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# Log file 'C:/Users/ayang/Documents/Research/Stage2_Strong0/log', created by user ayang
#mol addrep 0
set Trace1 2380
set Trace2 1749
set Trace3 2522
set Trace4 383
set Trace5 1460
set Trace6 1436
set Trace7 2068
set Trace8 2299


mol addrep 0
mol modselect 0 0 water
#mol modstyle 0 0 Points 1.000000
mol modstyle 0 0 CPK 2.500000 0.300000 12.000000 12.000000
mol modmaterial 0 0 AOShiny
mol color Name
#mol representation Points 1.000000
mol representation CPK 1.000000 0.300000 12.000000 12.000000
mol selection water
mol addrep 0
mol modselect 1 0 not water
mol modmaterial 1 0 Transparent
mol modstyle 1 0 Licorice 0.500000 15.000000 12.000000
mol color Name
mol representation Licorice 0.500000 15.000000 12.000000
mol selection not water
mol material Transparent
mol modrep 1 0
mol smoothrep 0 1 5
mol smoothrep 0 0 5

mol modselect 2 0 resid $Trace1
mol modmaterial 2 0 AOShiny
mol modstyle 2 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 2 0 Points 3.000000 
mol modcolor 2 0 ColorID 0
mol smoothrep 0 2 5
mol addrep 0
mol modselect 3 0 resid $Trace2
mol modmaterial 3 0 AOShiny
mol modstyle 3 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 3 0 Points 3.000000 
mol modcolor 3 0 ColorID 1
mol smoothrep 0 3 5
mol addrep 0
mol modselect 4 0 resid $Trace3
mol modmaterial 4 0 AOShiny
mol modstyle 4 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 4 0 Points 3.000000
mol modcolor 4 0 ColorID 2
mol smoothrep 0 4 5
#mol modmaterial Opaque
mol addrep 0
mol modselect 5 0 resid $Trace4
mol modmaterial 5 0 AOShiny
mol modstyle 5 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 5 0 Points 3.000000
mol modcolor 5 0 ColorID 3
mol smoothrep 0 5 5
#mol modmaterial Opaque
mol addrep 0
mol modselect 6 0 resid $Trace5
mol modmaterial 6 0 AOShiny
mol modstyle 6 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 6 0 Points 3.000000 
mol modcolor 6 0 ColorID 4
mol smoothrep 0 6 5
#mol modmaterial Opaque
mol addrep 0
mol modselect 7 0 resid $Trace6
mol modmaterial 7 0 AOShiny
mol modstyle 7 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 7 0 Points 3.000000
mol modcolor 7 0 ColorID 5
mol smoothrep 0 7 5
#mol modmaterial Opaque
mol addrep 0
mol modselect 8 0 resid $Trace7
mol modmaterial 8 0 AOShiny
mol modstyle 8 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 8 0 Points 3.000000
mol modcolor 8 0 ColorID 6
mol smoothrep 0 8 5
#mol modmaterial Opaque
mol addrep 0
mol modselect 9 0 resid $Trace8
mol modmaterial 9 0 AOShiny
mol modstyle 9 0 CPK 6.000000 0.300000 12.000000 12.000000
#mol modstyle 9 0 Points 3.000000
mol modcolor 9 0 ColorID 7
mol smoothrep 0 9 5
#mol modmaterial AOShiny
display projection orthographic
axes location off
scale to 0.03
pbc box
rotate x to 90
animate goto 5
# VMD for WIN32, version 1.9.3beta1 (July 21, 2016)
# end of log file.
