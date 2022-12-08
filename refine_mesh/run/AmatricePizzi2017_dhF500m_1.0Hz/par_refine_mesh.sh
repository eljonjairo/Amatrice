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

nb_iter=3              # number of iterations
prefxtet='mesh.'       # prefix of output tetgen files names

mesh_name='AmatricePizzi2017_dhF500m'
output_mesh='Tetgen_OutputMesh'
output_properties='Output_properties'

# ------------------ Files for TetGen ----------------------------
filepoly='./Poly_files/'${mesh_name}'.poly'            # .poly file name
filepolyvar='./Poly_files/'${mesh_name}'.var'          #.var file name (must be the same that filepoly)
fileout='par_refine_mesh.out'                             # refine mesh output file
str_name=${mesh_name}

# ----------------- Files for refine_mesh.f90 ----------------------
output_field='vs'                                   # Output field to visualize in Paraview
velf='ItalyCenterBianchi_1000m.txt'             # 3D velocity file name
xi=300000.0                            # coordinate of first x node
xf=400000.0                            # coordinate of last x node
yi=4680000.0                            # coordinate of first y node
yf=4780000.0                           # coordinate of last y node
zi=-60000.0                              # coordinate of max vertical node
zf=2842.33                           # coordinate of min vertical node
pml=0                                    # considering PML for refinement (1 = yes,0 = no)   
ncpml=10                               # number of cpml nodes 
dcpml=20000.0                             # cpml thickness
P1=100.0                                 # P1 limit for Pk adapt
P2=70.0                                  # P2 limit for Pk adapt

echo 'OK reading variables'

rm -f refine_mesh.out
rm -f $fileout
rm -f mesh.* *.mesh

cp $filepoly mesh.poly
cp $filepolyvar mesh.var

echo 'OK Copy .poly file'

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

mpirun -np 2 ./../../bin/PARTest_RefMesh >> $fileout <<EOF

1                               -> 1 if first iteration
$nb_iter                        -> number of iterations
mesh.1                          -> prefix of the files produced by TETGEN
$xi  $xf  $yi  $yf  $zi  $zf    -> extreme values of .poly  
$pml  $ncpml  $dcpml               -> =0(=1) == no (yes) cpml, number of elements, thickness of cpml 
$xi $yi $zi                     -> origin of the cartesian grid (xi,yi,zi)
0                               -> direction of axis z in file (0=upwards, 1=downwards)
$freq                           -> max frequency to consider
$ratio_first                    -> requested ratio for elements (ratio = lambda_min / tetra inradius) 
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

mpirun -np 2 ./../../bin/PARTest_RefMesh >> $fileout <<EOF

$iter                           -> 1 if first iteration
$nb_iter                        -> number of iterations
mesh.$iter                      -> prefix of the files produced by TETGEN
$xi  $xf  $yi  $yf  $zi  $zf    -> extreme values of .poly
$pml  $ncpml  $dcpml               -> =0(=1) == no (yes) cpml, number of elements, thickness of cpml
$xi $yi $zi                     -> origin of the cartesian grid (xi,yi,zi)
0                               -> direction of axis z in file (0=upwards, 1=downwards)
$freq                           -> max frequency to consider
$ratio_first                    -> requested ratio for elements (ratio = lambda_min / tetra inradius)
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
