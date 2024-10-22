program ed_ti_slab
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop
  integer                                       :: Nlat,Nlso,Nso
  integer                                       :: ilat,iorb,ispin
  logical                                       :: converged,PBC
  
  !Bath:
  integer                                       :: Nb
  real(8)                                       :: wmixing
  real(8),allocatable,dimension(:,:)            :: Bath
  
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal_
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal
  
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  complex(8),allocatable,dimension(:,:)         :: tiHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),dimension(:,:),allocatable            :: kpath
  integer                                       :: Nk,Lk,Ly,Nkpath
  real(8)                                       :: dens

  !gamma matrices and fields:
  integer                                       :: BAND_CHECK,iii
  complex(8),dimension(4,4)                     :: emat,soxmat,soymat,sozmat
  real(8)                                       :: e0,mh,lambda,dummymag,dummymag_rescale,SCALE_DISPERSION,MASS_OFFSET,radius
  complex(8),dimension(4,4)                     :: NHmat,GapOpeningMat,magmat
  real(8)                                       :: GapOpeningField,NHfield_up,NHfield_dw,DUMMY
      real(8)                            :: offset

  !misc
  logical                                       :: tridiag,lrsym,do_bands2d,do_eigenbands
  character(len=60)                             :: finput
  character(len=32)                             :: hkfile
  complex(8),dimension(:,:,:),allocatable       :: toconverge

  !mpi
  integer                                       :: comm,rank
  logical                                       :: master,DIAG_SIGMA

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  ! 
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED_ti_SLAB.conf')
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(Nkpath,"NKPATH",finput,default=501)
  call parse_input_variable(tridiag,"TRIDIAG",finput,default=.true.)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(GapOpeningField,"GAPOPENINGFIELD",finput,default=0.0d0)
  call parse_input_variable(NHfield_up,"NHFIELD_UP",finput,default=0.0d0)
  call parse_input_variable(NHfield_dw,"NHFIELD_DW",finput,default=0.0d0)
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(PBC,"PBC",finput,default=.false.)
  call parse_input_variable(lrsym,"LRSYM",finput,default=.true.)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(DIAG_SIGMA,"DIAG_SIGMA",finput,default=.false.)
  call parse_input_variable(MASS_OFFSET,"MASS_OFFSET",finput,default=-3.d0)
  call parse_input_variable(BAND_CHECK,"BAND_CHECK",finput,default=20)
  call parse_input_variable(SCALE_DISPERSION,"SCALE_DISPERSION",finput,default=0.1d0)
  call parse_input_variable(do_bands2d,"DO_BANDS2D",finput,default=.false.)
  call parse_input_variable(do_eigenbands,"DO_EIGENBANDS",finput,default=.false.)
  !
  call ed_read_input(trim(finput))

  !SETUP THE GAMMA MATRICES:
  emat = kron_pauli( pauli_sigma_0, pauli_tau_z)
  !
  soxmat=kron_pauli( pauli_sigma_z, pauli_tau_x)
  soymat=kron_pauli( pauli_sigma_0, pauli_tau_y)
  sozmat=kron_pauli( pauli_sigma_x, pauli_tau_x)
  magmat=kron_pauli( pauli_sigma_z, pauli_tau_0)
  !
  GapOpeningMat=kron_pauli( pauli_sigma_y, pauli_tau_X)
  NHmat =zero
  NHmat(1:2,1:2) = xi * NHFIELD_UP * pauli_tau_0
  NHmat(3:4,3:4) = (2.25 + xi) * NHFIELD_DW * pauli_tau_0
  !!

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  BAND_CHECK=LY*4
  print*,"BAND CHECK =",band_check
  XMU=2.25*NHFIELD_DW/2

  !set the global number of lattice sites equal to the number of layers along the y-axis
  Nlat = Ly
  
  !set the local number of total spin-orbitals (4)
  if(Nspin/=2.OR.Norb/=3)stop "Wrong setup from input file: Nspin=2, Norb=3 -> 6Spin-Orbitals"
  Nso  = Nspin*Norb

  !set the total lattice-spin-orbit dimension:
  Nlso=Nlat*Nspin*Norb
  !
  !Allocate Functions:
  !
  !Buil the Hamiltonian on a grid or on  path
  call build_hkr(trim(hkfile))
  Hloc = lso2nnn(tiHloc,Nlat,Nspin,Norb)
  !

  !
  if(do_eigenbands)then
    call build_eigenbands()

  endif
  
  call get_state_localization()
  
  call finalize_MPI()

contains


  !+---------------------------------------------------------------------------+!
  !PURPOSE: build the ti Hamiltonian H(k_x,kz,R_y) on the STRIPE along Y
  !+---------------------------------------------------------------------------+!
  subroutine build_hkr(file)
    character(len=*),optional          :: file
    integer                            :: i,ik
    !
    Lk=Nk**2
    !
    !SETUP THE H(kx,Ry,kz):
    if(master)then
       write(LOGfile,*)"Build H(kx,y,kz) for ti-stripe:"
       write(*,*)"# of kx and kz points     :",Nk
       write(*,*)"# of y-layers      :",Nlat
    endif
    !
    if(allocated(Hkr))deallocate(Hkr)
    allocate(Hkr(Nlso,Nlso,Lk))
    !
    call TB_set_bk([pi2,0d0,0d0],[0d0,pi2,0d0],[0d0,0d0,pi2])
    call TB_build_model(Hkr,ti_edge_model,Ly,Nso,[Nk,1,Nk],pbc=PBC,wdos=.false.)
    !   
    !
    !SETUP THE LOCAL PART Hloc(Ry)
    if(allocated(tiHloc))deallocate(tiHloc)
    allocate(tiHloc(Nlso,Nlso))
    tiHloc = sum(Hkr(:,:,:),dim=3)/Lk
    where(abs((tiHloc))<1.d-9)tiHloc=0d0
    !call TB_write_Hloc(tiHloc)
    !
  end subroutine build_hkr


  !+---------------------------------------------------------------------------+!
  !PURPOSE: solve H_ti(k_x,R_y,kz) along the 1d -pi:pi path in the BZ.
  !+---------------------------------------------------------------------------+!
  subroutine build_eigenbands(kpath_)
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    type(rgb_color),dimension(:,:),allocatable :: colors
    integer                            :: Npts,iii
    character(len=64)                  :: file
    !
    !PRINT H(kx,Ry) ALONG A -pi:pi PATH
    if(master)write(LOGfile,*)"Solve H(kx,y,kz) along [-Z:Z]:"
    Npts=3
    allocate(Kpath(Npts,3))
    offset=find_x_coordinate(20)/pi
    !
    !
    kpath(1,:)=[offset,0d0,0d0]*pi
    kpath(2,:)=[offset,0d0,1d0]*pi
    kpath(3,:)=[offset,0d0,2d0]*pi
    file="Eigenbands.nint"
    allocate(colors(Ly,Nso))
    colors=grey99
    do iii=1,Ly
      colors(iii,:) = [red1,red1,red1,blue1,blue1,blue1]
    enddo
    call solve_nh_model(ti_edge_model,Ly,Nso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "G","Z","G"],&
         file="Eigenbands.nint",&
         pbc=PBC)
  end subroutine build_eigenbands


  

  !+---------------------------------------------------------------------------+!
  !PURPOSE: bands_2d and other plots
  !+---------------------------------------------------------------------------+!


    subroutine get_state_localization()
    integer                                   :: ipts,icount,ic,unit,iorb,ilat,io,nrot,u1,u2,N,iband,iso
    real(8)                                   :: coeff(Nlat*Nso),coeffreduced(Nlat),siteindex(Nlat)
    complex(8),dimension(Nlat*Nso,Nlat*Nso)   :: h,evec
    complex(8),dimension(Nlat*Nso)            :: Eval
    real(8),dimension(Nlat*Nso)               :: RePartEval
    integer                                   :: labels(Nlat*Nso)
    !
    !
    do io=1,Nlat
      siteindex(io)=io
    enddo
    !
    icount=0
    do iband=BAND_CHECK-1,BAND_CHECK+1
      icount=icount+1
      coeffreduced=zero
      h = ti_edge_model([offset,0.d0,1.d0]*pi,Nlat,Nso,pbc) - XMU*Eye(Nlat*Nso)
      call eig(h,Eval,Evec)
      call sort_eigensystem(Eval,Evec)
      coeff=Evec(:,iband)*conjg(Evec(:,iband))
      coeff=coeff/sqrt(dot_product(Evec(:,iband),conjg(Evec(:,iband))))
      do io=1,Nlat
        do iso=1,Nso
          coeffreduced(io) = coeffreduced(io) + coeff((io-1)*Nso+iso)
        enddo
      enddo
      call splot("state_"//reg(txtfy(icount))//"_localization.dat",siteindex,coeffreduced) 
    enddo
    !
  end subroutine get_state_localization    



  !+---------------------------------------------------------------------------+!
  !PURPOSE: the ti-edge model hamiltonian
  !+---------------------------------------------------------------------------+!
  !ti on a stripe geometry;
  function ti_edge_model(kpoint,Nlat,N,pbc) result(Hrk)
    real(8),dimension(:)                :: kpoint
    real(8)                             :: kx,kz
    integer                             :: Nlat,N
    complex(8),dimension(N,N)           :: Hmat,Tmat,TmatH
    complex(8),dimension(Nlat*N,Nlat*N) :: Hrk
    integer                             :: i,Idmin,Idmax,Itmin,Itmax
    logical                             :: pbc
    kx=kpoint(1)
    kz=kpoint(3)
    Hrk=zero
    Hmat=h0_rk_ti(kx,kz,N)
    Tmat=t0_rk_ti(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nlat
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat 
    enddo
    do i=1,Nlat-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nlat-1)*N
       Itmax=0+Nlat*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
  end function ti_edge_model

  function h0_rk_ti(kx,kz,N) result(H)
    real(8)                    :: kx,kz
    integer                    :: N
    complex(8),dimension(4,4)  :: Hblock
    complex(8),dimension(N,N)  :: H
    !
    Hblock = (Mh - e0*(cos(kx) + cos(kz)))*emat+&
         lambda*(sin(kx)*soxmat + sin(kz)*sozmat)  + GapOpeningField*GapOpeningMat +Nhmat
    !
    H=zero
    !
    H(1,1)=(1+xi)*MASS_OFFSET*e0 - SCALE_DISPERSION*e0*(cos(kx) + cos(kz))
    H(4,4)=(1+xi)*MASS_OFFSET*e0 - SCALE_DISPERSION*e0*(cos(kx) + cos(kz))
    H(2:3,2:3)=Hblock(1:2,1:2)
    H(2:3,5:6)=Hblock(1:2,3:4)
    H(5:6,2:3)=Hblock(3:4,1:2)
    H(5:6,5:6)=Hblock(3:4,3:4)
  end function h0_rk_ti

  function t0_rk_ti(N) result(H)
    integer                    :: N
    complex(8),dimension(4,4)  :: Hblock
    complex(8),dimension(N,N)  :: H
    !
    Hblock = -0.5d0*e0*emat + xi*0.5d0*lambda*soymat
    H=0
    !
    H(1,1)=-SCALE_DISPERSION*0.5d0*e0
    H(2:3,2:3)=Hblock(1:2,1:2)
    H(2:3,5:6)=Hblock(1:2,3:4)
    H(4,4)=-SCALE_DISPERSION*0.5d0*e0
    H(5:6,2:3)=Hblock(3:4,1:2)
    H(5:6,5:6)=Hblock(3:4,3:4)
  end function T0_rk_ti


  function select_block(ip,Matrix) result(Vblock)
    integer                                          :: ip
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: Matrix
    complex(8),dimension(Nspin*Norb,Nspin*Norb)      :: Vblock
    integer                                          :: is,js,ispin,jspin,iorb,jorb
    Vblock=zero
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                is = iorb + (ispin-1)*Norb !spin-orbit stride
                js = jorb + (jspin-1)*Norb !spin-orbit stride
                Vblock(is,js) = Matrix(ip,ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function select_block

  
  function lso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn
  
  function so2j_index(ispin,iorb) result(isporb)
    integer :: ispin,iorb
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ispin-1)*Nspin + iorb
  end function so2j_index
  
  function so2j(fg) result(g)
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(i,j) = fg(ispin,jspin,iorb,jorb)
             enddo
          enddo
       enddo
    enddo
  end function so2j

  function j2so(fg) result(g)
    complex(8),dimension(Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nspin,Nspin,Norb,Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin
    do ispin=1,Nspin
       do jspin=1,Nspin
          do iorb=1,Norb
             do jorb=1,Norb
                i=so2j_index(ispin,iorb)
                j=so2j_index(jspin,jorb)
                g(ispin,jspin,iorb,jorb)  = fg(i,j)
             enddo
          enddo
       enddo
    enddo
  end function j2so
  
  
  subroutine solve_nh_model(hkr_model,Nlat,Nso,kpath,Nkpath,colors_name,points_name,file,pbc,iproject)
    interface 
       function hkr_model(kpoint,Nlat,Nso,pbc)
         real(8),dimension(:)                    :: kpoint
         integer                                 :: Nlat,Nso
         logical                                 :: pbc
         complex(8),dimension(Nlat*Nso,Nlat*Nso) :: hkr_model
       end function hkr_model
    end interface
    integer                                   :: Nlat,Nso,Nlso
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath,Nktot
    type(rgb_color),dimension(Nlat,Nso)       :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: pbc,iproject
    character(len=256)                        :: file_general,file_real,file_imag
    logical                                   :: pbc_,iproject_
    character(len=256)                        :: xtics
    integer                                   :: Npts,Ndim
    integer                                   :: ipts,ik,ic,unit,iorb,ilat,io,nrot,u1,u2
    real(8)                                   :: coeff(Nlat*Nso),klen,ktics(size(Kpath,1))
    type(rgb_color)                           :: corb(Nlat*Nso),c(Nlat*Nso)
    real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff,bk_x,bk_y,bk_z
    complex(8),dimension(Nlat*Nso,Nlat*Nso)   :: h,evec
    complex(8),dimension(Nlat*Nso)            :: Eval
    integer,allocatable                       :: labels(:,:)
    real(8),allocatable                       :: kseg(:)
    complex(8),allocatable                    :: Ekval(:,:)
    integer,allocatable                       :: Ekcol(:,:)
    !
    !
    file_general    = "Eigenbands.tb";if(present(file))file_general=file
    file_real=reg(file_general)//".real"
    file_imag=reg(file_general)//".imag"
    iproject_= .false.        ;if(present(iproject))iproject_=iproject
    pbc_     = .true.         ;if(present(pbc))pbc_=pbc
    !
    Nlso  = Nlat*Nso
    Npts  = size(kpath,1)
    Ndim  = size(kpath,2)
    Nktot = (Npts-1)*Nkpath
    !
    do ilat=1,Nlat
       do io=1,Nso
          corb(io + (ilat-1)*Nso) = colors_name(ilat,io)
       enddo
    enddo
    !
    bk_x=[pi2,0d0,0d0]
    bk_y=[0d0,pi2,0d0]
    bk_z=[0d0,0d0,pi2]
    !
    if(iproject_)then
       select case(Ndim)
       case (1)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x
       case(2)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y
       case (3)
          forall(ipts=1:Npts)kpath(ipts,:) = kpath(ipts,1)*bk_x + kpath(ipts,2)*bk_y + kpath(ipts,3)*bk_z
       end select
    endif
    !
    if(master)then
       write(*,*)"Solving model along the path:"
       do ipts=1,Npts
          write(*,"(A,10(A,1x),A1)")"Point"//str(ipts)//": [",(str(kpath(ipts,ic)),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic=0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nlso))
    allocate(ekcol(Nktot,Nlso))
    klen = 0d0
    call start_timer()
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/Nkpath
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hkr_model(kpoint,Nlat,Nso,pbc) - XMU*eye(Nlat*Nso)
          call eig(h,Eval,Evec)
          call sort_eigensystem(Eval,Evec)
          call eta(ic,Nktot)
          do io=1,Nlso
             coeff(:)=Evec(:,io)*conjg(Evec(:,io))
             c(io) = coeff.dot.corb
             Ekval(ic,io) = Eval(io)
             Ekcol(ic,io) = rgb(c(io))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = Kseg(ic-1)
    call stop_timer()
    !
     open(free_unit(unit),file=str(file_real))
     do io=1,Nlso
        do ic=1,Nktot
           write(unit,*)kseg(ic),real(Ekval(ic,io)),Ekcol(ic,io)
        enddo
        write(unit,*)""
     enddo
     close(unit)
     !
     !
     xtics=""
     xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
     do ipts=2,Npts-1
        xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
     enddo
     xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
     !
     open(unit,file=reg(file_real)//".gp")
     write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
     write(unit,*)"#set out '"//reg(file_real)//".png'"
     write(unit,*)""
     write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
     write(unit,*)"#set out '"//reg(file_real)//".svg'"
     write(unit,*)""
     write(unit,*)"#set term postscript eps enhanced color 'Times'"
     write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_real)//".pdf'"
     write(unit,*)"unset key"
     write(unit,*)"set xtics ("//reg(xtics)//")"
     write(unit,*)"set grid ytics xtics"
     !
     write(unit,*)"plot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb variable"
     write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
     !
     close(unit)
     !
     call system("chmod +x "//reg(file_real)//".gp")
     !
    ic=0
    !allocate(kseg(Nktot))
    !allocate(ekval(Nktot,Nlso))
    !allocate(ekcol(Nktot,Nlso))
    klen = 0d0
    call start_timer()
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/Nkpath
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hkr_model(kpoint,Nlat,Nso,pbc) - XMU*eye(Nlat*Nso)
          call eig(h,Eval,Evec)
          Eval=Imag(Eval)+xi*real(Eval)
          call sort_eigensystem(Eval,Evec)
          call eta(ic,Nktot)
          do io=1,Nlso
             coeff(:)=Evec(:,io)*conjg(Evec(:,io))
             c(io) = coeff.dot.corb
             Ekval(ic,io) = Eval(io)
             Ekcol(ic,io) = rgb(c(io))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = Kseg(ic-1)
    call stop_timer()
    !
     open(free_unit(unit),file=str(file_imag))
     do io=1,Nlso
        do ic=1,Nktot
           write(unit,*)kseg(ic),real(Ekval(ic,io)),Ekcol(ic,io)
        enddo
        write(unit,*)""
     enddo
     close(unit)
     !
     !
     xtics=""
     xtics="'"//reg(points_name(1))//"'"//str(ktics(1))//","
     do ipts=2,Npts-1
        xtics=reg(xtics)//"'"//reg(points_name(ipts))//"'"//str(ktics(ipts))//","
     enddo
     xtics=reg(xtics)//"'"//reg(points_name(Npts))//"'"//str(ktics(Npts))//""
     !
     open(unit,file=reg(file_imag)//".gp")
     write(unit,*)"#set terminal pngcairo size 350,262 enhanced font 'Verdana,10'"
     write(unit,*)"#set out '"//reg(file_imag)//".png'"
     write(unit,*)""
     write(unit,*)"#set terminal svg size 350,262 fname 'Verdana, Helvetica, Arial, sans-serif'"
     write(unit,*)"#set out '"//reg(file_imag)//".svg'"
     write(unit,*)""
     write(unit,*)"#set term postscript eps enhanced color 'Times'"
     write(unit,*)"#set output '|ps2pdf  -dEPSCrop - "//reg(file_imag)//".pdf'"
     write(unit,*)"unset key"
     write(unit,*)"set xtics ("//reg(xtics)//")"
     write(unit,*)"set grid ytics xtics"
     !
     write(unit,*)"plot '"//reg(file_imag)//"' every :::0 u 1:2:3 w l lw 3 lc rgb variable"
     write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
     !
     close(unit)
     !
     call system("chmod +x "//reg(file_imag)//".gp")
  end subroutine solve_nh_model
    


  subroutine sort_eigensystem(eigval,eigvec)
      integer                                   :: ii
      integer,dimension(Nlat*Nso)               :: labels
      real(8),dimension(Nlat*Nso)               :: re_eigval
      complex(8),dimension(Nlat*Nso)            :: eigval,eigval_tmp
      complex(8),dimension(Nlat*Nso,Nlat*Nso)   :: eigvec,eigvec_tmp
      !
      labels=zero
      eigval_tmp=zero
      eigvec_tmp=zero
      !
      re_eigval=REAL(eigval)
      call sort_array(re_eigval,labels)
      !
      do ii=1,Nlat*Nso
        eigval_tmp(ii)=eigval(labels(ii))
        eigvec_tmp(:,ii)=eigvec(:,labels(ii))
      enddo
      !
      eigvec=eigvec_tmp
      eigval=eigval_tmp
      !
  end subroutine sort_eigensystem
  
  function lazy_sort(arr_in) result (arr)
    REAL(8), DIMENSION(:)                :: arr_in
    REAL(8), DIMENSION(:), allocatable   :: arr
    INTEGER                              :: i,j,inc,n
    REAL(8)                              :: v
    n=size(arr_in)
    allocate(arr(n))
    arr=arr_in
    inc=1
    do
      inc=3*inc+1
      if (inc > n) exit
    end do
    do
      inc=inc/3
      do i=inc+1,n
        v=arr(i)
        j=i
        do
          if (arr(j-inc) <= v) exit
          arr(j)=arr(j-inc)
          j=j-inc
          if (j <= inc) exit
        end do
        arr(j)=v
      end do
      if (inc <= 1) exit
    end do
  end function lazy_sort
  
  
    function find_x_coordinate(Nintervals) result(xcoord)
    integer                           :: interval,Nintervals
    real(8)                           :: xcoord,xmin,xmax,xcoord_tmp,yval_tmp,yval
    !
    write(LOGfile,*)"Searching for the kx offset:"
    write(LOGfile,*)"Using ",Nintervals," intervals"
    !
    yval=1000
    do interval = 0,Nintervals-1
      xmin = -pi + 2*pi/Nintervals*interval
      xmax = -pi + 2*pi/Nintervals*(interval+1)
      xcoord_tmp=(xmax+xmin)/2d0
      call  brent(gap_minimum,xcoord_tmp,[xmin,xcoord_tmp,xmax])
      yval_tmp=abs(gap_minimum(xcoord_tmp))
      print*,"from ",xmin," to ",xmax," ymin is ",yval_tmp," at ",xcoord_tmp
      if(yval_tmp.lt.yval)then 
        xcoord=xcoord_tmp
        yval=yval_tmp
      endif
    enddo
    write(LOGfile,*)"The EP stand at kx = ",xcoord
    write(LOGfile,*)"The real energy gap is ",yval
  end function find_x_coordinate

  function gap_minimum(xcoord) result (gap)
    real(8),dimension(3)                                        :: kpoint
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: Hrk,evec
    real(8)                                                     :: gap,xcoord
    integer                                                     :: coeff
    complex(8),dimension(Nlat*Nso)                              :: Eval,eval_sorted
    !
    kpoint=[xcoord,0d0,pi]
    Hrk=ti_edge_model(kpoint,Nlat,Nso,.false.)
    call eig(hrk,Eval,Evec)
    Eval_sorted=lazy_sort(REAL(Eval))
    gap=abs(Eval_sorted(BAND_CHECK + 1)-Eval_sorted(BAND_CHECK))
  end function gap_minimum

  function get_e(xcoord) result (gap)
    real(8),dimension(3)                                        :: kpoint
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb)       :: Hrk,evec
    real(8)                                                     :: gap,xcoord
    complex(8),dimension(Nlat*Nso)                              :: Eval,eval_sorted
    !
    kpoint=[xcoord,0d0,pi]
    Hrk=ti_edge_model(kpoint,Nlat,Nso,.false.)
    call eig(hrk,Eval,Evec)
    Eval_sorted=lazy_sort(REAL(Eval))
    gap=(Eval_sorted(2*Nlat+1)+Eval_sorted(2*Nlat))/2d0
  end function get_e
  

end program ed_ti_slab
