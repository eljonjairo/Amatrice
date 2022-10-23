
 PROGRAM PARTEST_REFMESH 

 USE bspline_sub_module        

 IMPLICIT NONE

 INCLUDE 'mpif.h'

 INTEGER :: nprocs, iproc, ierror, tag, status(MPI_STATUS_SIZE)   
 INTEGER :: iter, niter, pml, n_pml, zaxis, len_name, nnode, mtetra, ntetra, ltetra, ntetraori 
 INTEGER :: i, j, k, icont, ii, jj, kk, inode, itetra, jtetra, ypos, zpos, strlen, iprint, nprint
 INTEGER :: ivel, nvel, nx_vel, ny_vel, nz_vel, ineg, nb_P0, nb_P1, nb_P2, cont_ratio, n_bad
 INTEGER :: inbvx, inbvy, inbvz, iloy, iloz, ncell, ncell_sou, ncell_pml, ncell_control, ncell_control2
 INTEGER :: locncell, locncell_pml, locncell_control, locncell_control2, loccont_ratio, locn_bad
 INTEGER :: locineg, locnb_P0, locnb_P1, locnb_P2
 INTEGER, PARAMETER :: kx=2, ky=2, kz=2, iknot=0, idx=0, idy=0, idz=0, master=0

 REAL :: xxmin, xxmax, yymin, yymax, zzmin, zzmax, cpml, xor, yor, zor
 REAL :: fmax, ratio_lim, ratio_surf, l_surf, lim_P1, lim_P2, xmin, xmax, ymin, ymax, zmin, zmax
 REAL :: opt_radius, opt_vol, opt_len, l_pml, vol_pml, cfl, dt, locdt_min, dt_min, locxmin, locxmax, locymin, locymax
 REAL :: locmin_radius, locmax_radius, locmin_ratio, locmax_ratio, locmean_ratio, loctot_vol, locmin_vol, locmax_vol 
 REAL :: vel_ztop, vel_zbottom, x(4), y(4), z(4), min_vol, max_vol, volume, min_surf, max_surf, surface, locmin_surf, locmax_surf
 REAL :: inradius, min_ratio, max_ratio, min_radius, max_radius, min_len, max_len, tot_vol, qp, qs, loczmin, loczmax 
 REAL :: x21, y21, z21, x31, y31, z31, vnx, vny, vnz, ratio, lambda_min, mean_ratio
 REAL*8 :: xgt, ygt, zgt, vp, vs, rho, zgt_tmp      ! Because of bspline module these variables have to be real*8 

 CHARACTER(LEN=80) :: infile, velmodel3D_name, output_opt, infile1, infile2, infile3,infile4, outfile1, outfile2

 INTEGER, DIMENSION(:), POINTER :: TETRA1, TETRA2, TETRA3, TETRA4 ! tetra (tetra index, node number)
 INTEGER, DIMENSION(:), POINTER :: LOCTETRA1, LOCTETRA2, LOCTETRA3, LOCTETRA4 ! tetra (tetra index, node number)

 LOGICAL :: vol_check, arithmetic, ind_pml

 REAL, DIMENSION(:), POINTER :: NODEX, NODEY, NODEZ ! node (node index, cartesian coordinate)
 REAL, DIMENSION(:), POINTER :: VP_F, VS_F, RHO_F, QS_F, QP_F, BADTETRA ! physical properties (cell index)
 REAL, DIMENSION(:), POINTER :: LOCVP_F, LOCVS_F, LOCRHO_F, LOCQS_F, LOCQP_F, LOCBADTETRA ! physical properties (cell index)
 REAL, DIMENSION(:), POINTER :: LOCOUTVOL, OUTVOL ! Array for Output Volume File
 REAL*8, DIMENSION(:), POINTER :: XVEL, YVEL, ZVEL ! coordinate vectors  of the velocity properties (cell index)
 REAL*8, DIMENSION(:), POINTER :: TX, TY, TZ ! Just bspline arrays are real*8, the arrays use it for MPI have to be real
 
 REAL*8, DIMENSION(:,:,:), POINTER :: VP_MODEL3D, VS_MODEL3D, RHO_MODEL3D, BCOEF_VP, BCOEF_VS, BCOEF_RHO
 REAL, DIMENSION(:,:), POINTER :: VELPOINTS ! coordinates and physical properties of the 3D velocity model

 INTEGER, DIMENSION(3) :: int_flag, eval_flag

 nprint = 5                                  ! Number of array elements to be printed for checking

 CALL MPI_INIT(ierror)
 CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
 CALL MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierror)

 IF ( iproc==master ) THEN
   WRITE(6,*) '==============================='
   WRITE(6,*) '      PAR_REFINE_MESH'
   WRITE(6,*) '==============================='
   
   ! read input parameters
   READ(*,*) iter               ! First iteration?
   READ(*,*) niter              ! Number of iterations
   READ(*,*) infile             ! prefix of the files produced by TETGEN
   READ(*,*) xxmin, xxmax, yymin, yymax, zzmin, zzmax !extreme values of .poly
   READ(*,*) pml, n_pml, cpml   ! =0(=1) == no (yes) cpml, number of elements, thick of cpml 
   READ(*,*) xor, yor, zor      ! origin of the cartesian grid
   READ(*,*) zaxis              ! direction of axis z in file (0=upwards, 1=downwards)
   READ(*,*) fmax               ! max frequency to consider
   READ(*,*) ratio_lim          ! requested ratio for elements (ratio = lambda_min / tetra inradius) 
   READ(*,*) lim_P1, lim_P2     ! P1 and P2 limits for Pk adapt
   READ(*,*) velmodel3D_name    ! file name of 3D velocity model
   READ(*,*) output_opt         ! Field to write in vtk file (vp,vs,rho,qp,qs,badtetra)

  !=====================================
  !-------------------------------------
  !     Step 1 Read the mesh file
  !-------------------------------------
  !=====================================

   WRITE(*,*) ''
   WRITE(*,*) ' Read file in TETGEN format'
   WRITE(*,*) '----------------------------'
 ENDIF

 ! Share parameters

 CALL MPI_BCAST(iter,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(niter,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror) 
 strlen = LEN(infile)
 CALL MPI_BCAST(infile,strlen,MPI_CHARACTER,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(xxmin,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(yymin,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(zzmin,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(xxmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(yymax,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(zzmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(pml,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(n_pml,1,MPI_INTEGER,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(cpml,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 strlen = LEN(velmodel3D_name)
 CALL MPI_BCAST(velmodel3D_name,strlen,MPI_CHARACTER,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(ratio_lim,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(fmax,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(lim_P1,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 
 CALL MPI_BCAST(lim_P2,1,MPI_REAL,master,MPI_COMM_WORLD,ierror) 

 len_name = LEN_TRIM(infile)
 infile1  = infile(1:len_name) // '.node'
 infile2  = infile(1:len_name) // '.ele'
 infile3  = infile(1:len_name) // '.var'
 infile4  = infile(1:len_name) // '.vtk'
 outfile1 = infile(1:len_name) // '.vol'
 outfile2 = infile(1:len_name) // '.mtr'

 OPEN(UNIT=1, FILE=infile1)

 !-------
 ! NODES
 !-------
 READ(1,*) nnode

 ! allocate node table
 ALLOCATE( NODEX(nnode), NODEY(nnode), NODEZ(nnode) )
 ! read node coordinates
 DO inode =1, nnode
   READ(1,*) ii, NODEX(inode), NODEY(inode), NODEZ(inode)
 ENDDO
 CLOSE(1)

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
 IF ( iproc==master ) THEN
   WRITE(6,*) ''
   WRITE(6,*) ' Read nodes OK, nb nodes', nnode

   OPEN(UNIT=2, FILE=infile2)
   READ(2,*) ntetraori

   mtetra = MOD(ntetraori,nprocs)
   ntetra = ntetraori+nprocs-mtetra
   ltetra = ntetra/nprocs

   ! Add mterra elements to scatter arrays of the same size, in the last proccesor the last mtetra elements are the same
   ! first mtetra elements, these repeated elements dont have to be included in the output files

   ! allocate tetra table
   ALLOCATE( TETRA1(ntetra), TETRA2(ntetra), TETRA3(ntetra), TETRA4(ntetra), OUTVOL(ntetra) )
   ALLOCATE( VP_F(ntetra), VS_F(ntetra), RHO_F(ntetra), QS_F(ntetra), QP_F(ntetra), BADTETRA(ntetra) );

   ! read elements characteristics
   DO itetra =1, ntetraori
     READ(2,*)  ii, TETRA1(itetra), TETRA2(itetra), TETRA3(itetra), TETRA4(itetra) 
   ENDDO
   CLOSE(2)

   ! Repeate the first melements
   jtetra = 0
   DO itetra = ntetraori+1,ntetra
     jtetra = jtetra+1
     TETRA1(itetra) = TETRA1(jtetra)
     TETRA2(itetra) = TETRA2(jtetra)
     TETRA3(itetra) = TETRA3(jtetra)
     TETRA4(itetra) = TETRA4(jtetra)
   END DO
   WRITE(6,*) ''
   WRITE(6,*) ' Read tetrahedra OK, nb tetrahedra', ntetraori
 ENDIF

 CALL MPI_BCAST(ltetra,1,MPI_INTEGER,master,MPI_COMM_WORLD, ierror) 
 CALL MPI_BCAST(mtetra,1,MPI_INTEGER,master,MPI_COMM_WORLD, ierror) 

 ! allocate tetra table
 ALLOCATE( LOCTETRA1(ltetra), LOCTETRA2(ltetra), LOCTETRA3(ltetra), LOCTETRA4(ltetra), LOCOUTVOL(ltetra) )
 ALLOCATE( LOCVP_F(ltetra), LOCVS_F(ltetra), LOCRHO_F(ltetra), LOCQS_F(ltetra), LOCQP_F(ltetra), LOCBADTETRA(ltetra) );

 ! Scatter TETRAS array 
 CALL MPI_SCATTER(TETRA1,ltetra,MPI_INTEGER,LOCTETRA1,ltetra,MPI_INTEGER,master,MPI_COMM_WORLD, ierror) 
 CALL MPI_SCATTER(TETRA2,ltetra,MPI_INTEGER,LOCTETRA2,ltetra,MPI_INTEGER,master,MPI_COMM_WORLD, ierror) 
 CALL MPI_SCATTER(TETRA3,ltetra,MPI_INTEGER,LOCTETRA3,ltetra,MPI_INTEGER,master,MPI_COMM_WORLD, ierror) 
 CALL MPI_SCATTER(TETRA4,ltetra,MPI_INTEGER,LOCTETRA4,ltetra,MPI_INTEGER,master,MPI_COMM_WORLD, ierror) 

 IF ( iproc==master ) THEN
   WRITE(*,*) '----------------------------'
   WRITE(*,*) ' Reading 3D velocity model  '
   WRITE(*,*) '----------------------------'
 ENDIF

 !print*,velmodel3D_name
 OPEN(UNIT=42, FILE=velmodel3D_name)

 ! Reading number of velocity model points and lengths of x, y and z vectors
 READ(42,*) nvel,nx_vel,ny_vel,nz_vel

 ALLOCATE( VELPOINTS(nvel,6) )

 ! Reading coordinates and velocity properties of the 3D model
 DO ivel = 1,nvel
   READ(42,*) (VELPOINTS(ivel,jj), jj=1,6)
 ENDDO
 CLOSE(42)

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
 IF ( iproc==master ) THEN
   WRITE(*,*) ''
   WRITE(*,*) ' Read velocity structure'
   WRITE(*,*) ' Number or vel points', nvel
   WRITE(*,*) ' Number or x points', nx_vel
   WRITE(*,*) ' Number or y points', ny_vel
   WRITE(*,*) ' Number or z points', nz_vel
   WRITE(*,*) ''
 END IF

 ALLOCATE( TX(nx_vel+kx), TY(ny_vel+ky), TZ(nz_vel+kz) )
 ALLOCATE( VP_MODEL3D(nx_vel,ny_vel,nz_vel), VS_MODEL3D(nx_vel,ny_vel,nz_vel), RHO_MODEL3D(nx_vel,ny_vel,nz_vel) )
 ALLOCATE( BCOEF_VP(nx_vel,ny_vel,nz_vel), BCOEF_VS(nx_vel,ny_vel,nz_vel), BCOEF_RHO(nx_vel,ny_vel,nz_vel) )
 ALLOCATE( XVEL(nx_vel), YVEL(ny_vel), ZVEL(nz_vel) )

 ! Arranging 3D velocity matrixes to compute bspline interpolant
 DO i = 1,nx_vel
   XVEL(i) = VELPOINTS(i,1)
   DO j = 1,ny_vel
     ypos = nx_vel*(j-1) + 1
     YVEL(j) = VELPOINTS(ypos,2)
     DO k = 1,nz_vel
       zpos = nx_vel*ny_vel*(k-1) + 1
       ZVEL(k) = VELPOINTS(zpos,3)
     ENDDO
   ENDDO
 ENDDO

 icont = 0
 DO k = 1,nz_vel
   DO j = 1,ny_vel
     DO i = 1,nx_vel
       icont = icont + 1
       VP_MODEL3D(i,j,k) = VELPOINTS(icont,4)
       VS_MODEL3D(i,j,k) = VELPOINTS(icont,5)
       RHO_MODEL3D(i,j,k) = VELPOINTS(icont,6)
     ENDDO
   ENDDO
 ENDDO

 vel_ztop = MAXVAL(zvel)
 vel_zbottom = MINVAL(zvel)

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)
 IF ( iproc==master ) THEN
   WRITE(6,*) '---------------------------------------'
   WRITE(6,*) ' Read velocity model OK'
   WRITE(6,*) '---------------------------------------'

   WRITE(6,*) ' Read Limits of zVel vector from vel file'
   WRITE(6,*) 'zvel_top',vel_ztop,'zvel_bottom',vel_zbottom
   WRITE(6,*) '-------------------------------------------'
 END IF

 !Computing 3D interpolant for velocity properties
 call db3ink( XVEL,nx_vel,YVEL,ny_vel,ZVEL,nz_vel,VP_MODEL3D,kx,ky,kz,iknot,TX,TY,TZ,BCOEF_VP,int_flag(1) )
 call db3ink( XVEL,nx_vel,YVEL,ny_vel,ZVEL,nz_vel,VS_MODEL3D,kx,ky,kz,iknot,TX,TY,TZ,BCOEF_VS,int_flag(2) )
 call db3ink( XVEL,nx_vel,YVEL,ny_vel,ZVEL,nz_vel,RHO_MODEL3D,kx,ky,kz,iknot,TX,TY,TZ,BCOEF_RHO,int_flag(3) )

 if (any(int_flag/=0)) then
   do i=1,3
     if (int_flag(i)/=0) then
       write(6,*) 'Error initializing 3D spline in :',i,' property'
       write(6,*) 'Check with "get_status_messag(',int_flag(i),')" function'
     end if
   end do
 end if

 IF ( iproc==master ) THEN
   WRITE(6,*) '---------------------------------------------------'
   WRITE(6,*) ' Computing 3D interpolant for velocity properties, OK'
   WRITE(6,*) '---------------------------------------------------'
   WRITE(*,*) ''
 END IF

  !have to set these before the first evaluate call in the interpolation:
 inbvx = 1
 inbvy = 1
 inbvz = 1

 iloy  = 1
 iloz  = 1

 !==============================
 !------------------------------
 !     Step 4 mesh checking
 !------------------------------
 !==============================

 IF ( iproc==master ) THEN
  WRITE(*,*) ''
  WRITE(*,*) ' Mesh checking'
  WRITE(*,*) '---------------'
  WRITE(*,*) 'ITERATION = ',ITER
 END IF

 !-----------------------------------------------------
 ! output file contains the ideal volume for each cell
 !-----------------------------------------------------

 locncell = 0
 locncell_pml = 0
 locncell_control = 0
 locncell_control2 = 0
 locineg = 0

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

  ! initialization
 locmin_vol    = huge(min_vol)
 locmax_vol    = -huge(max_vol)
 min_surf      = huge(min_surf)
 max_surf      = -huge(max_surf)
 locmin_ratio  = huge(min_surf)
 locmax_ratio  = -huge(max_surf)
 locmin_radius = huge(min_surf)
 locmax_radius = -huge(max_surf)
 min_len    = huge(min_len)
 max_len    = -huge(max_len)
 loctot_vol = 0.
 locxmax    = -huge(xmax)
 locymax    = -huge(xmax)
 loczmax    = -huge(xmax)
 locxmin    = huge(xmax)
 locymin    = huge(xmax)
 loczmin    = huge(xmax)
 locdt_min  = huge(xmax)
 locnb_P0   = 0
 locnb_P1   = 0
 locnb_P2   = 0

 locn_bad      = 0
 loccont_ratio = 0

 DO itetra = 1, ltetra

   ! retrieve node
   x(1) = NODEX(LOCTETRA1(itetra))
   y(1) = NODEY(LOCTETRA1(itetra))
   z(1) = NODEZ(LOCTETRA1(itetra))
  
   ! extremum coordinates
   locxmax = MAX(locxmax, x(1))
   locymax = MAX(locymax, y(1))
   loczmax = MAX(loczmax, z(1))

   locxmin = MIN(locxmin, x(1))
   locymin = MIN(locymin, y(1))
   loczmin = MIN(loczmin, z(1))

   ! retrieve node
   x(2) = NODEX(LOCTETRA2(itetra))
   y(2) = NODEY(LOCTETRA2(itetra))
   z(2) = NODEZ(LOCTETRA2(itetra))

   ! extremum coordinates
   locxmax = MAX(locxmax, x(2))
   locymax = MAX(locymax, y(2))
   loczmax = MAX(loczmax, z(2))

   locxmin = MIN(locxmin, x(2))
   locymin = MIN(locymin, y(2))
   loczmin = MIN(loczmin, z(2))

   ! retrieve node
   x(3) = NODEX(LOCTETRA3(itetra))
   y(3) = NODEY(LOCTETRA3(itetra))
   z(3) = NODEZ(LOCTETRA3(itetra))

   ! extremum coordinates
   locxmax = MAX(locxmax, x(3))
   locymax = MAX(locymax, y(3))
   loczmax = MAX(loczmax, z(3))

   locxmin = MIN(locxmin, x(3))
   locymin = MIN(locymin, y(3))
   loczmin = MIN(loczmin, z(3))

   ! retrieve node
   x(4) = NODEX(LOCTETRA4(itetra))
   y(4) = NODEY(LOCTETRA4(itetra))
   z(4) = NODEZ(LOCTETRA4(itetra))

   ! extremum coordinates
   locxmax = MAX(locxmax, x(4))
   locymax = MAX(locymax, y(4))
   loczmax = MAX(loczmax, z(4))

   locxmin = MIN(locxmin, x(4))
   locymin = MIN(locymin, y(4))
   loczmin = MIN(loczmin, z(4))

   ! compute volume
   volume = ABS ( (x(2)-x(1)) * ( (y(3)-y(1)) * (z(4)-z(1)) - (y(4)-y(1)) * (z(3)-z(1)) ) + &
                  (y(2)-y(1)) * ( (x(4)-x(1)) * (z(3)-z(1)) - (z(4)-z(1)) * (x(3)-x(1)) ) + &
                  (z(2)-z(1)) * ( (x(3)-x(1)) * (y(4)-y(1)) - (x(4)-x(1)) * (y(3)-y(1)) ) ) / 6.

  ! total volume (for checking)
   loctot_vol = loctot_vol + volume

   ! extremum volume
   locmin_vol = MIN(locmin_vol, volume)
   locmax_vol = MAX(locmax_vol, volume)

   ! compute surface
   surface = 0.

   ! face 1 (nodes 123)
   ii = 1
   jj = 2
   kk = 3

   x21 = x(ii) - x(jj)
   y21 = y(ii) - y(jj)
   z21 = z(ii) - z(jj)

   x31 = x(ii) - x(kk)
   y31 = y(ii) - y(kk)
   z31 = z(ii) - z(kk)

   vnx = 0.5*(z21*y31 - y21*z31)
   vny = 0.5*(x21*z31 - z21*x31)
   vnz = 0.5*(y21*x31 - x21*y31)
   surface = surface + SQRT(vnx**2 + vny**2 + vnz**2)

   ! face 2 (nodes 134)
   ii = 1
   jj = 3
   kk = 4

   x21 = x(ii) - x(jj)
   y21 = y(ii) - y(jj)
   z21 = z(ii) - z(jj)

   x31 = x(ii) - x(kk)
   y31 = y(ii) - y(kk)
   z31 = z(ii) - z(kk)

   vnx = 0.5*(z21*y31 - y21*z31)
   vny = 0.5*(x21*z31 - z21*x31)
   vnz = 0.5*(y21*x31 - x21*y31)

   surface = surface + SQRT(vnx**2 + vny**2 + vnz**2)

   ! face 3 (nodes 142)
   ii = 1
   jj = 4
   kk = 2

   x21 = x(ii) - x(jj)
   y21 = y(ii) - y(jj)
   z21 = z(ii) - z(jj)

   x31 = x(ii) - x(kk)
   y31 = y(ii) - y(kk)
   z31 = z(ii) - z(kk)

   vnx = 0.5*(z21*y31 - y21*z31)
   vny = 0.5*(x21*z31 - z21*x31)
   vnz = 0.5*(y21*x31 - x21*y31)

   surface = surface + SQRT(vnx**2 + vny**2 + vnz**2)
   ! face 4 (nodes 234)
   ii = 2
   jj = 3
   kk = 4

   x21 = x(ii) - x(jj)
   y21 = y(ii) - y(jj)
   z21 = z(ii) - z(jj)

   x31 = x(ii) - x(kk)
   y31 = y(ii) - y(kk)
   z31 = z(ii) - z(kk)

   vnx = 0.5*(z21*y31 - y21*z31)
   vny = 0.5*(x21*z31 - z21*x31)
   vnz = 0.5*(y21*x31 - x21*y31)

   surface = surface + SQRT(vnx**2 + vny**2 + vnz**2)

   min_surf = MIN(min_surf, surface)
   max_surf = MAX(max_surf, surface)

   !------------------
   ! compute inradius
   !------------------
   inradius = 3. * volume / surface

   locmin_radius = MIN(locmin_radius, inradius)
   locmax_radius = MAX(locmax_radius, inradius)
   !==============================================
   !
   ! Retrieve the physical property for each cell
   !
   !==============================================

   ! Modified by Victor
   ! Determining the centroid of the i- tetrahedrum

   xgt = SUM(x(1:4))/4.
   ygt = SUM(y(1:4))/4.
!     zgt = SUM(z(1:4))/4.
!     zgt_temp = SUM(z(1:4))/4.
   IF (ITER <= 7) THEN
     zgt = MAXVAL(z(1:4))
     ELSE
     zgt = SUM(z(1:4))/4.
   END IF

!! Modified by Carlos
   IF (zgt >= vel_ztop) THEN
      zgt_tmp = vel_ztop !Then, the interpolator will consider the velocity of the nearest upper boundary element 
     ELSEIF (zgt <= vel_zbottom) THEN
       zgt_tmp = vel_zbottom !Then, the interpolator will consider the velocity of the nearest bottom boundary element
     ELSE
       zgt_tmp = zgt
   END IF

   ! Evaluate the interpolant for vp, vs and rho
   call db3val(xgt,ygt,zgt_tmp,idx,idy,idz,tx,ty,tz,nx_vel,ny_vel,nz_vel,kx,ky,kz,&
               bcoef_vp,vp,eval_flag(1),inbvx,inbvy,inbvz,iloy,iloz)

   call db3val(xgt,ygt,zgt_tmp,idx,idy,idz,tx,ty,tz,nx_vel,ny_vel,nz_vel,kx,ky,kz,&
               bcoef_vs,vs,eval_flag(2),inbvx,inbvy,inbvz,iloy,iloz)

   call db3val(xgt,ygt,zgt_tmp,idx,idy,idz,tx,ty,tz,nx_vel,ny_vel,nz_vel,kx,ky,kz,&
               bcoef_rho,rho,eval_flag(3),inbvx,inbvy,inbvz,iloy,iloz)

   IF (any(eval_flag/=0)) THEN
     DO i = 1,3
       IF (eval_flag(i)/=0) THEN
         WRITE(*,*) 'Error evaluating 3D spline in : ',i,' property for the ',itetra,' element'
         WRITE(6,*) 'Check with "get_status_messag(',eval_flag(i),')" function'
       END IF
     END DO
   END IF

   ! check if vs = 0 (water)
   IF (vs /= 0.) THEN
     lambda_min = vs / fmax
   ELSE
     lambda_min = vp / fmax
   ENDIF

   IF (vs < 0.) THEN
     locineg = locineg + 1
     WRITE(*,*) iproc,itetra,locineg,xgt,ygt,zgt,vs
   ENDIF

   ! Central Italy Qp and Qs Ameri 2012 Table 2.
   IF (zgt <= 1000)  THEN
     qp = 200
     qs = 100

    ELSEIF ((zgt <= 27000) .AND. (zgt > 1000)) THEN
      qp = 400
      qs = 200

    ELSEIF ((zgt <= 42000) .AND. (zgt > 27000)) THEN
      qp = 300
      qs = 600

    ELSE
      qp = 800
      qs = 400
   ENDIF

   ! Built the arrays of properties for each element (vs,vp,rho,qs,qp)
   LOCVP_F(itetra)  = vp
   LOCVS_F(itetra)  = vs
   LOCRHO_F(itetra) = rho
   LOCQP_F(itetra)  = qp
   LOCQS_F(itetra)  = qs

   !----------------------------------------------------------
   !              COMPUTING OPTIMAL RADIUS
   !----------------------------------------------------------
   opt_radius = lambda_min / ratio_lim
   ! Control to know the number of elements in the effective domain (no cmpl)
   IF (((xgt >= xxmin+cpml .AND. xgt <= xxmax-cpml) .AND. (ygt >= yymin+cpml .AND.&
     ygt <= yymax-cpml)) .AND. zgt > zzmin+cpml)  THEN
     locncell_control = locncell_control+1
   ENDIF

   ! -------------------------------------------------------------------
   !      Determine  elements to refine in the useful domain. (no cpml)
   ! --------------------------------------------------------------------

   ! check inradius
   IF (inradius > opt_radius) THEN
     locncell = locncell + 1
     opt_len = SQRT(24.) * opt_radius ! Regular tetrahedron (Wikipedia)
     !opt_vol = SQRT(2.)/12. * (opt_len**3.) ! Regular tetrahedron (Wiki)
     opt_vol = volume * 0.9
    ELSE
      opt_vol = -1.
   ENDIF

   ! Determining wether the element belongs to CPML region (just for statitsitcs)
   ind_pml = .FALSE.

   IF (xgt <= xxmin+cpml .OR. xgt >= xxmax-cpml) THEN
     locncell_control2 = locncell_control2+1
     ind_pml = .TRUE.
   ENDIF
   IF (ygt <= yymin+cpml .OR. ygt >= yymax-cpml) THEN
     locncell_control2 = locncell_control2+1
     ind_pml = .TRUE.
   ENDIF
   IF (zgt <= zzmin+cpml) THEN
     locncell_control2 = locncell_control2+1
     ind_pml = .TRUE.
   ENDIF

   ! -------------------------------------------------------------------
   !      Determine  elements to refine in cpml
   ! --------------------------------------------------------------------

   IF (ind_pml) THEN

     l_pml = cpml/n_pml
     !vol_pml = SQRT(3.)/10.*(l_pml**3.) ! (Previous option)
     vol_pml = SQRT(2.)/12. * (l_pml**3.) ! Regular tetrahedron (Wikipedia)

     IF(vol_pml < volume) THEN
       locncell_pml = locncell_pml + 1
       IF (pml == 1) THEN
         opt_vol = volume*0.9
       ENDIF
       IF (pml == 1) THEN
         opt_vol = -1
       ENDIF
     ENDIF
   ENDIF

   IF(opt_vol > 0.) THEN
     LOCBADTETRA(itetra) = 1.
    ELSE
      LOCBADTETRA(itetra) = 0.
   ENDIF

   LOCOUTVOL(itetra) = opt_vol

   ratio = lambda_min / inradius

   locmin_ratio  = MIN(locmin_ratio, ratio)
   locmax_ratio  = MAX(locmax_ratio, ratio)
   locmean_ratio = locmean_ratio + ratio

   IF (locmax_ratio > 20) THEN
     loccont_ratio = loccont_ratio+1
   ENDIF

   IF (ratio < 3.) THEN ! To respect at least the P2 criterion
     locn_bad = locn_bad+1
   ENDIF

   ! -----------------------------------------------------------------------
   !                       Retrieving min dt
   ! ----------------------------------------------------------------------

   ! taken into account adapt interpolation order

   IF (ratio > lim_P1) THEN
     cfl   = 1.
     locnb_P0 = locnb_P0 + 1
    ELSEIF ((ratio <= lim_P1) .AND. (ratio > lim_P2)) THEN
      cfl   = 1/3.
      locnb_P1 = locnb_P1 + 1
    ELSE
     cfl   = 1/5.
     locnb_P2 = locnb_P2 + 1
   ENDIF

   ! compute time step
   dt        = cfl * 2. * inradius / vp
   locdt_min = MIN(locdt_min, dt)

 END DO

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

 ! Sending local results to master
 CALL MPI_GATHER(LOCBADTETRA,ltetra,MPI_REAL,BADTETRA,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)
 CALL MPI_GATHER(LOCOUTVOL,ltetra,MPI_REAL,OUTVOL,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)
 CALL MPI_GATHER(LOCVP_F,ltetra,MPI_REAL,VP_F,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)
 CALL MPI_GATHER(LOCVS_F,ltetra,MPI_REAL,VS_F,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)
 CALL MPI_GATHER(LOCRHO_F,ltetra,MPI_REAL,RHO_F,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)
 CALL MPI_GATHER(LOCQP_F,ltetra,MPI_REAL,QP_F,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)
 CALL MPI_GATHER(LOCQS_F,ltetra,MPI_REAL,QS_F,ltetra,MPI_REAL,master,MPI_COMM_WORLD,ierror)

 CALL MPI_REDUCE(locdt_min,dt_min,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locn_bad,n_bad,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locncell_control,ncell_control,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locncell_control2,ncell_control2,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locncell,ncell,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locncell_pml,ncell_pml,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmin_vol,min_vol,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmax_vol,max_vol,1,MPI_REAL,MPI_MAX,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(loctot_vol,tot_vol,1,MPI_REAL,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmin_radius,min_radius,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmax_radius,max_radius,1,MPI_REAL,MPI_MAX,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locxmin,xmin,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locymin,ymin,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(loczmin,zmin,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locxmax,xmax,1,MPI_REAL,MPI_MAX,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locymax,ymax,1,MPI_REAL,MPI_MAX,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(loczmax,zmax,1,MPI_REAL,MPI_MAX,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmin_ratio,min_ratio,1,MPI_REAL,MPI_MIN,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmax_ratio,max_ratio,1,MPI_REAL,MPI_MAX,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locmean_ratio,mean_ratio,1,MPI_REAL,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(loccont_ratio,cont_ratio,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locnb_P0,nb_P0,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locnb_P1,nb_P1,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locnb_P2,nb_P2,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 CALL MPI_REDUCE(locineg,ineg,1,MPI_INTEGER,MPI_SUM,master,MPI_COMM_WORLD,ierror)
 
! DO ii = 1, nprint
!   print*,iproc,ii,LOCVS_F(ii)
! END DO
! print*
!
! IF(iproc==master) THEN
!   DO ii = 1,nprocs      
!     jj = (ii-1)*ltetra
!     DO kk = 1, nprint
!       iprint = jj+kk
!       print*,ii-1,iprint,VS_F(iprint)
!     END DO
!   END DO  
! END IF         

 IF(iproc==master) THEN

   OPEN(UNIT=13, FILE=outfile1)
   WRITE(13,*) ntetraori
   DO itetra = 1,ntetraori
     WRITE(13,*) itetra,OUTVOL(itetra)
   END DO
   CLOSE(13)

   WRITE(6,*) ''
   WRITE(6,*) ' Frequency to mesh', fmax
   WRITE(6,*) ' Elements per wavelength', ratio_lim

   WRITE(6,*) ''
   WRITE(6,*) ' Nb bad-shaped tetrahedra', n_bad

   WRITE(6,*) ''
   WRITE(6,*) ' Nb tetrahedra', ntetra
   WRITE(6,*) 'ITERACION', infile

   WRITE(6,*) ''
   WRITE(6,100) ncell_control, (REAL(ncell_control)/REAL(ntetra)) * 100
   100  FORMAT('  Nb cells in the useful domain ', i10, ' [', f5.2, ' %]')

   WRITE(6,*) ''
   WRITE(6,300) ncell, (REAL(ncell)/REAL(ntetra)) * 100
   300  FORMAT('  Nb cells of useful domain to refine ', i10, ' [', f5.2, ' %]')

   WRITE(6,*) ''
   WRITE(6,200) ncell_control2, (REAL(ncell_control2)/REAL(ntetra)) * 100
   200  FORMAT('  Nb cells in the PML ', i10, ' [', f5.2, ' %]')

   WRITE(6,*) ''
   WRITE(6,400) ncell_pml, (REAL(ncell_pml)/REAL(ntetra)) * 100
   400  FORMAT('  Nb cells of PML to refine ', i10, ' [', f5.2, ' %]')

   WRITE(6,*) ''
   WRITE(6,*) ' Min vol', min_vol
   WRITE(6,*) ' Max vol', max_vol
   WRITE(6,*) ' Tot vol', tot_vol
   WRITE(6,*) ''
   WRITE(6,*) ' Min radius', min_radius
   WRITE(6,*) ' Max radius', max_radius
   WRITE(6,*) ''
   WRITE(6,*) ' Min VP', MINVAL(VP_F(:))
   WRITE(6,*) ' Max VP', MAXVAL(VP_F(:))
   WRITE(6,*) ''
   WRITE(6,*) ' Min VS', MINVAL(VS_F(:))
   WRITE(6,*) ' Max VS', MAXVAL(VS_F(:))
   WRITE(6,*) ''
   WRITE(6,*) ' Min QP', MINVAL(QP_F(:))
   WRITE(6,*) ' Max QP', MAXVAL(QP_F(:))
   WRITE(6,*) ''
   WRITE(6,*) ' Min QS', MINVAL(QS_F(:))
   WRITE(6,*) ' Max QS', MAXVAL(QS_F(:))
   WRITE(6,*) ''
   WRITE(6,*) ' xmin / xmax', xmin, xmax
   WRITE(6,*) ' ymin / ymax', ymin, ymax
   WRITE(6,*) ' zmin / zmax', zmin, zmax
   WRITE(6,*) ''
   WRITE(6,*) ' Min ratio ', min_ratio
   WRITE(6,*) ' Max ratio ', max_ratio
   WRITE(6,*) ' Mean ratio', mean_ratio/ntetra
   WRITE(6,*) ' Nb Elements with ratio > 20',cont_ratio
   WRITE(6,*) ''
   WRITE(6,*) ' Min dt ', dt_min
   WRITE(6,*) ''
   WRITE(6,*) ' nb P0 ', nb_P0
   WRITE(6,*) ' nb P1 ', nb_P1
   WRITE(6,*) ' nb P2 ', nb_P2
   WRITE(6,*) ''
   WRITE(6,*) ' nb Negative Vs ', ineg

     ! write physical files
   IF (ITER .GE. NITER) THEN

     WRITE(*,*) ''
     WRITE(*,*) ' Write physical files'
     WRITE(*,*) '----------------------'

     OPEN(UNIT=14, FILE='vp.mesh', ACCESS='DIRECT', FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=4)
     OPEN(UNIT=15, FILE='vs.mesh', ACCESS='DIRECT', FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=4)
     OPEN(UNIT=16, FILE='rho.mesh', ACCESS='DIRECT', FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=4)
     OPEN(UNIT=17, FILE='qs.mesh', ACCESS='DIRECT', FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=4)
     OPEN(UNIT=18, FILE='qp.mesh', ACCESS='DIRECT', FORM='UNFORMATTED', STATUS='UNKNOWN', RECL=4)

     DO ii = 1, ntetraori
       WRITE(14, rec=ii) REAL(VP_F(ii),4)
       WRITE(15, rec=ii) REAL(VS_F(ii),4)
       WRITE(16, rec=ii) REAL(RHO_F(ii),4)
       WRITE(17, rec=ii) REAL(QS_F(ii),4)
       WRITE(18, rec=ii) REAL(QP_F(ii),4)
     ENDDO

     CLOSE(14)
     CLOSE(15)
     CLOSE(16)
     CLOSE(17)
     CLOSE(18)

     IF (output_opt .EQ. 'vp') THEN

       WRITE(*,*) '----------------------------------'
       WRITE(*,*) ' Write Vp for the .vtk file'
       WRITE(*,*) '----------------------------------'

       OPEN(UNIT=19,FILE=infile4,STATUS='OLD',POSITION='APPEND')
       WRITE(19,*)'CELL_DATA',ntetraori
       WRITE(19,*)'SCALARS vp float'
       WRITE(19,*)'LOOKUP_TABLE default'

       DO ii = 1, ntetraori
         WRITE(19,*) VP_F(ii)
       ENDDO
       CLOSE(19)

       WRITE(*,*) ''
       WRITE(*,*) ' OK!'
       WRITE(*,*) '----------------------------------'

     ENDIF

     IF (output_opt .EQ. 'vs') THEN

       WRITE(*,*) '----------------------------------'
       WRITE(*,*) ' Write Vs for the .vtk file'
       WRITE(*,*) '----------------------------------'

       OPEN(UNIT=19,FILE=infile4,STATUS='OLD',POSITION='APPEND')
       WRITE(19,*)'CELL_DATA',ntetraori
       WRITE(19,*)'SCALARS vs float'
       WRITE(19,*)'LOOKUP_TABLE default'

       DO ii = 1, ntetraori
         WRITE(19,*) VS_F(ii)
       ENDDO

       WRITE(*,*) ''
       WRITE(*,*) ' OK!'
       WRITE(*,*) '----------------------------------'

       CLOSE(19)

     ENDIF

     IF (output_opt .EQ. 'rho') THEN

       WRITE(*,*) ''
       WRITE(*,*) ' Write Rho in .vtk file OK'
       WRITE(*,*) '----------------------------------'

       OPEN(UNIT=19,FILE=infile4,STATUS='OLD',POSITION='APPEND')
       WRITE(19,*)'CELL_DATA',ntetraori
       WRITE(19,*)'SCALARS rho float'
       WRITE(19,*)'LOOKUP_TABLE default'

       DO ii = 1, ntetraori
         WRITE(19,*) RHO_F(ii)
       ENDDO

       CLOSE(19)

       WRITE(*,*) ''
       WRITE(*,*) ' Ok!'
       WRITE(*,*) '----------------------------------'

     ENDIF


     IF (output_opt .EQ. 'qp') THEN

       WRITE(*,*) ''
       WRITE(*,*) ' Write qp in .vtk file OK'
       WRITE(*,*) '----------------------------------'

       OPEN(UNIT=19,FILE=infile4,STATUS='OLD',POSITION='APPEND')
       WRITE(19,*)'CELL_DATA',ntetraori
       WRITE(19,*)'SCALARS qp float'
       WRITE(19,*)'LOOKUP_TABLE default'

       DO ii = 1, ntetraori
         WRITE(19,*) QP_F(ii)
       ENDDO

       CLOSE(19)

       WRITE(*,*) ''
       WRITE(*,*) ' Ok!'
       WRITE(*,*) '----------------------------------'

     ENDIF

     IF (output_opt .EQ. 'qs') THEN

       WRITE(*,*) ''
       WRITE(*,*) ' Write qs in .vtk file OK'
       WRITE(*,*) '----------------------------------'

       OPEN(UNIT=19,FILE=infile4,STATUS='OLD',POSITION='APPEND')
       WRITE(19,*)'CELL_DATA',ntetraori
       WRITE(19,*)'SCALARS qs float'
       WRITE(19,*)'LOOKUP_TABLE default'

       DO ii = 1, ntetraori
         WRITE(19,*) QS_F(ii)
       ENDDO

       CLOSE(19)

       WRITE(*,*) ''
       WRITE(*,*) ' Ok!'
       WRITE(*,*) '----------------------------------'

     ENDIF

     IF (output_opt .EQ. 'badtetra') THEN

       WRITE(*,*) ''
       WRITE(*,*) ' Write badtetra in .vtk file OK'
       WRITE(*,*) '----------------------------------'

       OPEN(UNIT=19,FILE=infile4,STATUS='OLD',POSITION='APPEND')
       WRITE(19,*)'CELL_DATA',ntetraori
       WRITE(19,*)'SCALARS badtetra float'
       WRITE(19,*)'LOOKUP_TABLE default'
       DO ii = 1, ntetraori
         WRITE(19,*) BADTETRA(ii)
       ENDDO

       CLOSE(19)

       WRITE(*,*) ''
       WRITE(*,*) ' Ok!'
       WRITE(*,*) '----------------------------------'

     ENDIF

   ENDIF

   ! deallocate tables  
   DEALLOCATE( TETRA1,TETRA2,TETRA3,TETRA4,OUTVOL )
   DEALLOCATE( VP_F,VS_F,RHO_F,QS_F,QP_F,BADTETRA );

 END IF         
 
 DEALLOCATE( NODEX, NODEY, NODEZ )
 DEALLOCATE( LOCTETRA1, LOCTETRA2, LOCTETRA3, LOCTETRA4, LOCOUTVOL )
 DEALLOCATE( LOCVP_F, LOCVS_F, LOCRHO_F, LOCQS_F, LOCQP_F, LOCBADTETRA );
 DEALLOCATE( VELPOINTS )
 DEALLOCATE( TX, TY, TZ )
 DEALLOCATE( VP_MODEL3D, VS_MODEL3D, RHO_MODEL3D )
 DEALLOCATE( BCOEF_VP, BCOEF_VS, BCOEF_RHO )
 DEALLOCATE( XVEL, YVEL, ZVEL )

 CALL MPI_BARRIER(MPI_COMM_WORLD, ierror)

 IF(iproc==master) THEN
   WRITE(*,*) '==============================='
   WRITE(*,*) '     REFINE_MESH ENDDED OK'
   WRITE(*,*) '==============================='
 END IF

 CALL MPI_FINALIZE(ierror)

 END PROGRAM PARTEST_REFMESH 
