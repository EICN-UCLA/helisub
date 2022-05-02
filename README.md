# helisub

For compiling:  g++ -I. lib/librelion_lib.a lib/librelion_gui_lib.a -o helisub helisub.C lib/librelion_lib.a lib/librelion_gui_lib.a

(helisub is the program after compiling)

By typing ./helisub, you will see the usage of the program: 

Usage: helisub <in> <out> <apix> <hrise> <hturn> <vectorX> <vectorY> <vectorZ> <offset> <#take> <boxsize> [<csym>]

in: input star file from 3D classification
  
out: root name of the output star file and image stack
  
apix: pixel size
 
hrise: helical rise in angstrom
  
hturn: helical turn angle in degree
  
vectorX, vectorY, vectorZ: vector from center of the helix to the first sub-particle you want to locate. If the sub-particle is on the X axis, the vector is (R, 0,0). R is the radius of the helix. 
  
offset: 0
  
#take: number of sub-particles in each helical segment
  
boxsize: size of the sub-particle
  
csym: symmetry of the sub-particle (optional)
  
