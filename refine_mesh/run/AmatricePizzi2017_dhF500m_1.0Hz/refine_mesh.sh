#--------------------------------------------------------------
#
# Build & refine mesh for the Guerrero Subduction Zone 
#
#--------------------------------------------------------------

#!/bin/sh

# ------------------------------------------------------
#               Variables definitions
# ------------------------------------------------------

# -------------- Refinement parameters ---------------------------

ratio_first=3.0        # *** ratio (lambda / cell radius) 
ratio_second=5.0
freq=1.05	       # frequency

nb_iter=5              # number of iterations
prefxtet='mesh.'       # prefix of output tetgen files names

mesh_name='/AmatricePizzi2017_dhF500m'
topo_res='1000'    # Resolution for topography
fault_res='500'    # Resolution for .var file
bin_res='100'      # Resolution of binary files
output_mesh='Tetgen_OutputMesh'
output_properties='Output_properties'

# ------------------ Files for TetGen ----------------------------

filepoly='./Poly_files/'${mesh_name}'.poly'            # .poly file name
filepolyvar='./Poly_files/'${mesh_name}'.var'          #.var file name (must be the same that filepoly)
fileout='refine_mesh.out'                             # refine mesh output file
str_name=${mesh_name}

# ----------------- Files for refine_mesh.f90 ----------------------
toposup_bin='./Bin_files/TopoSup_'${mesh_name}'_'${bin_res}'m.bin'         # file name of inland topography
topobasin_bin='./Bin_files/TopoBasin_'${mesh_name}'_'${bin_res}'m.bin'     # file name of bottom topography (basin)
topotmvb_bin='./Bin_files/TopoTMVB_'${mesh_name}'_'${bin_res}'m.bin'       # file name of TMVB topography (inferior)
output_field='vs'                                   # Output field to visualize in Paraview

velf='ItalyCenterBianchi_1000m.txt'             # 3D velocity file name
nx=76                                   # number of nodes in x direction binary files
ny=76                                   # number of nodes in y direction binary files
dh=100.                                 # distance between nodes binary files
xi=-50000.0                             # coordinate of first x node
xf=50000.0                             # coordinate of last x node
yi=-50000.0                             # coordinate of first y node
yf=50000.0                             # coordinate of last y node
zi=-40000.0                              # coordinate of max vertical node
zf=0.0                               # coordinate of min vertical node
srcx=0.0                                 # x coordinate of the source
srcy=0.0                                 # y coordinate of the source
srcz=-15000.0                            # z coordinate of the source
pml=0                                    # considering PML for refinement (1 = yes,0 = no)   
ncpml=10                               # number of cpml nodes 
dcpml=10000.0                             # cpml thickness
zref=50000.0                            # depth of surface refinement
P1=100.0                                 # P1 limit for Pk adapt
P2=70.0                                  # P2 limit for Pk adapt

echo 'OK reading variables'

rm -f refine_mesh.out
rm -f $fileout
rm -f mesh.* *.mesh

cp $filepoly mesh.poly
cp $filepolyvar mesh.var

echo 'OK Copy .poly file'

#cp $toposup_bin TopoSup.bin
#cp $topobasin_bin TopoBasin.bin
#cp $topotmvb_bin TopoTMVB.bin

#echo 'Ok Copy .bin files'

echo '---------------------------------------'
echo '    Build initial mesh with TETGEN     '
echo '---------------------------------------'
echo ''

#./../../../BIN/tetgen -pq1.4a5e9 mesh.poly >> $fileout
tetgen -pq1.4a5e9 mesh.poly >> $fileout

echo '          First mesh DONE!             '
echo '---------------------------------------'
echo ''
echo '     Calling first time refine_mesh    '
echo '---------------------------------------'

#./../bin/refine_mesh_Acambay >> $fileout <<EOF
./../../bin/refine_mesh_LAquila >> $fileout <<EOF

1                               -> 1 if first iteration
$nb_iter                        -> number of iterations
mesh.1                          -> prefix of the files produced by TETGEN
'TopoSup.bin'                   -> file name of inland topography
'TopoBasin.bin'                 -> file name of bottom topography (basin)
'TopoTMVB.bin'                  -> file name of TMVB topography (inferior)
$nx     $ny     $dh             -> nb points in x and y, and distance between grid points 
$xi  $xf  $yi  $yf  $zi  $zf    -> extreme values of .poly  
0 150. 500.                     -> =0(=1) == no (yes) source, size of elements, delta around. 
$srcx  $srcy  $srcz             -> source coordinates 
$pml  $ncpml  $dcpml               -> =0(=1) == no (yes) cpml, number of elements, thickness of cpml 
$xi $yi $zi                     -> origin of the cartesian grid (xi,yi,zi)
0                               -> direction of axis z in file (0=upwards, 1=downwards)
$freq                           -> max frequency to consider
$ratio_first                    -> requested ratio for elements (ratio = lambda_min / tetra inradius) 
0.                              -> requested ratio for subsurface refinement (ratio = lambda_min / tetra inradius) 
$zref                           -> depth of the subsurface refinement
$P1  $P2                        -> P1 and P2 limits for Pk adapt
$velf                           -> file name of 3D velocity model
$output_field                   -> field to write in vtk file
EOF

iter=1

# ------------------------------------------------------------
#        Loop for refinenment mesh to n-iterations
# ------------------------------------------------------------

until [ ! "$iter" -le $nb_iter ]

do
echo '------------------------------------'
echo '    Iteration ' $iter ' of ' $nb_iter 
echo '------------------------------------'
echo ''

#./tetgen -rq1.2ak $prefxtet$iter >> $fileout
#./../../../BIN/tetgen -rqak $prefxtet$iter >> $fileout 
tetgen -rqak $prefxtet$iter >> $fileout 

iter=$(expr $iter + 1)


#./../bin/refine_mesh_Acambay >> $fileout <<EOF
./../../bin/refine_mesh_LAquila >> $fileout <<EOF

$iter                               -> 1 if first iteration
$nb_iter                        -> number of iterations
mesh.$iter                          -> prefix of the files produced by TETGEN
'TopoSup.bin'                   -> file name of inland topography
'TopoBasin.bin'                 -> file name of bottom topography (basin)
'TopoTMVB.bin'                  -> file name of TMVB topography (inferior)
$nx     $ny     $dh             -> nb points in x and y, and distance between grid points 
$xi  $xf  $yi  $yf  $zi  $zf    -> extreme values of .poly  
0 150. 500.                     -> =0(=1) == no (yes) source, size of elements, delta around. 
$srcx  $srcy  $srcz             -> source coordinates 
$pml  $ncpml  $dcpml               -> =0(=1) == no (yes) cpml, number of elements, thickness of cpml 
$xi $yi $zi                     -> origin of the cartesian grid (xi,yi,zi)
0                               -> direction of axis z in file (0=upwards, 1=downwards)
$freq                           -> max frequency to consider
$ratio_second                    -> requested ratio for elements (ratio = lambda_min / tetra inradius) 
0.                              -> requested ratio for subsurface refinement (ratio = lambda_min / tetra inradius) 
$zref                           -> depth of the subsurface refinement
$P1  $P2                        -> P1 and P2 limits for Pk adapt
$velf                           -> file name of 3D velocity model
$output_field                   -> field to write in vtk file
EOF

if [ $iter -ge $nb_iter ]
then
mv vp.mesh vp.$iter.mesh
mv vs.mesh vs.$iter.mesh
mv rho.mesh rho.$iter.mesh
mv qs.mesh qs.$iter.mesh
mv qp.mesh qp.$iter.mesh
fi

if [ $iter -ge 2 ]
then

# To remove previous iteration outputs
iter=$(expr $iter - 1 )
rm *.$iter.*

iter=$(expr $iter + 1 )
fi

done

echo 'This is the final number of iter: ' $iter

rm -f *.edge
rm -f *.vol

mv vp.$iter.mesh $output_properties/vp_$str_name.mesh
mv vs.$iter.mesh $output_properties/vs_$str_name.mesh
mv rho.$iter.mesh $output_properties/rho_$str_name.mesh
mv qs.$iter.mesh $output_properties/qs_$str_name.mesh
mv qp.$iter.mesh $output_properties/qp_$str_name.mesh

mv mesh.$iter.node $output_mesh/mesh_$str_name.node
mv mesh.$iter.face $output_mesh/mesh_$str_name.face
mv mesh.$iter.ele  $output_mesh/mesh_$str_name.ele
mv mesh.$iter.vtk  $output_mesh/mesh_$str_name.vtk

exit

