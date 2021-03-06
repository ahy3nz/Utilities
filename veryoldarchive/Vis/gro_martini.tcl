# VMD for MACOSXX86, version 1.9.1 (February 4, 2012)
# Log file '/Users/ahy3nz/Programs/setup/CoarseGraining/CG_Vis', created by user ahy3nz
source ~/Programs/Utilities/Vis/cg_bonds.tcl

mol modselect 0 0 resid 1 to 64
mol modstyle 0 0 Licorice 1.000000 15.000000 12.000000
mol color Name
mol representation VDW 1.000000 12.000000
mol selection resid 1 to 64

mol addrep 0
mol modselect 1 0 resid 65 to 128
mol modstyle 1 0 Licorice 1.000000 15.000000 12.000000
mol color Name
mol selection resid 65 to 128

mol addrep 0
mol modselect 2 0 resid 129 to 10000
mol modcolor 2 0 ColorID 1
mol modmaterial 2 0 Glass2
mol color ColorID 1
mol representation VDW 1.000000 12.000000
#mol selection resid 129 to 768
#mol modrep 2 0
mol smoothrep 0 1 5
mol smoothrep 0 0 5
mol smoothrep 0 2 5

display resetview
rotate x to 90.000000
scale to 0.02
pbc box
display projection Orthographic
axes location off
# VMD for MACOSXX86, version 1.9.1 (February 4, 2012)
# end of log file.
