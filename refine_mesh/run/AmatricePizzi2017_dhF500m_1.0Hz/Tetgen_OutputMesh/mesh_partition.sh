\rm -f GEODG3D.meshDATA GEODG3D.meshPART

ulimit -s unlimited

./mesh_partition <<EOF
mesh_AmatricePizzi2017_dhF500m -> file from mesher
2        -> mesher (0=GMSH, 1=LAGRIT, 2=TETGEN) 
128       -> nb of requested partitions
0        -> graph type for partition (0=nodal, 1=dual graph)
0        -> weighted partioning (0=no, 1=yes)
EOF
