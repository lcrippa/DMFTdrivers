program soc_test_model
  USE SCIFOR
  USE DMFT_TOOLS
  USE EDIPACK2 
  USE MPI
  implicit none

  integer                                       :: iloop,Nb,Lk,Nx,Nso,ik,iorb,irepl,le,iii,jjj,NsymMats
  logical                                       :: converged
  real(8),dimension(:,:),allocatable            :: Dbands
  real(8),dimension(:,:),allocatable            :: Ebands
  real(8)                                       :: wmixing,onsite,lambda,orb_shift,t0,t1,t2,t3,t4
  real(8),dimension(:),allocatable              :: ts,Dband
  real(8),dimension(:),allocatable              :: de,dens
  real(8),dimension(:),allocatable              :: Wband,wfreq,eps_array,wm,wr
  !Bath:
  real(8),allocatable                           :: Bath(:),Bath_(:)
  real(8),dimension(:,:),allocatable            :: lambdasym_vectors
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_t2g
  complex(8),dimension(:,:,:,:),allocatable     :: spinorbit_matrix_d_nn,spinorbit_matrix_t2g_nn
  complex(8),dimension(:,:),allocatable         :: spinorbit_matrix_d,spinorbit_matrix_t2g, u_jbasis
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc_imp(:,:),Hloc_latt(:,:),tmpmat(:,:)
  complex(8),allocatable                        :: Hloc_imp_nn(:,:,:,:),Hloc_latt_nn(:,:,:,:)
  real(8),dimension(:),allocatable              :: H0 !elements on the diagonal of Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:)     :: S0
  real(8),allocatable,dimension(:,:)            :: Zmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkmats
  complex(8),allocatable,dimension(:)           :: Gtest
  complex(8),allocatable,dimension(:,:)         :: tmpmat_so
  !
  character(len=32)                             :: finput,foutput,model,PARAMETRIZATION
  complex(8),allocatable                        :: Hk(:,:,:)
  real(8),allocatable                           :: Wt(:)
  !
  integer                                       :: comm,rank,unit
  logical                                       :: master
  logical                                       :: mixG0,symOrbs,DOSFLAG,ZJBASIS


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  
  t0 = 0.42
  t1 = 0.17
  t2 = 0.30
  t3 = 0.03
  t3 = 0.04
  

  call parse_cmd_variable(finput,"FINPUT",default='input.in')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(DOSFLAG,"DOSFLAG",finput,default=.false.)
  call parse_input_variable(ZJBASIS,"ZJBASIS",finput,default=.false.)
  call parse_input_variable(PARAMETRIZATION,"PARAMETRIZATION",finput,default="diagonal")
  call parse_input_variable(Le,"LE",finput,default=500)
  call parse_input_variable(Nx,"Nx",finput,default=100,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  call parse_input_variable(lambda,"lambda",finput,default=0d0,comment="soc parameter")
  call parse_input_variable(orb_shift,"orb_shift",finput,default=0d0,comment="xy orbital shift")
  !
  !
  call ed_read_input(trim(finput))
  !
  !

  !
  Nso=Nspin*Norb
  allocate(TS(Nso))

  
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)
  

  !
  call parse_input_variable(ts,"TS",finput,default=(/(0.5d0,iii=1,size(TS))/),comment="hopping parameter")

  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(nbath,"nbath")
  call add_ctrl_var(ed_hw_bath,"ed_hw_bath")

  !The spin-orbit matrix, initialized by hand

  if(Nspin/=2.OR.Norb/=3)then
    allocate(spinorbit_matrix_d(10,10));spinorbit_matrix_d=zero
    allocate(spinorbit_matrix_d_nn(Nspin,Nspin,5,5));spinorbit_matrix_d_nn=zero
    allocate(spinorbit_matrix_t2g(2*Norb,2*Norb));spinorbit_matrix_t2g=zero
    allocate(spinorbit_matrix_t2g_nn(Nspin,Nspin,Norb,Norb));spinorbit_matrix_t2g_nn=zero
    ZJBASIS=.false.
    print*,"Wrong setup from input file: Nspin=2 Norb=3 -> 6Spin-Orbitals: You can't use spinorbit here. Everything will be zero, Hlambda will be dummy"
  else
    allocate(spinorbit_matrix_d(10,10));spinorbit_matrix_d=zero
    allocate(spinorbit_matrix_d_nn(Nspin,Nspin,5,5));spinorbit_matrix_d_nn=zero
    allocate(spinorbit_matrix_t2g(6,6));spinorbit_matrix_t2g=zero
    allocate(u_jbasis(6,6));u_jbasis=zero
    allocate(spinorbit_matrix_t2g_nn(Nspin,Nspin,3,3));spinorbit_matrix_t2g_nn=zero
    !
    spinorbit_matrix_t2g(1,:) = lambda*(-0.5)*[  zero,   -xi,   zero,   zero,  zero,   one]
    spinorbit_matrix_t2g(2,:) = lambda*(-0.5)*[    xi,  zero,   zero,   zero,  zero,   -xi]
    spinorbit_matrix_t2g(3,:) = lambda*(-0.5)*[  zero,  zero,   zero,   -one,    xi,  zero]
    spinorbit_matrix_t2g(4,:) = lambda*(-0.5)*[  zero,  zero,   -one,   zero,    xi,  zero]
    spinorbit_matrix_t2g(5,:) = lambda*(-0.5)*[  zero,  zero,    -xi,    -xi,  zero,  zero]  
    spinorbit_matrix_t2g(6,:) = lambda*(-0.5)*[   one,   xi,    zero,   zero,  zero,  zero] 
    spinorbit_matrix_t2g_nn = so2nn(spinorbit_matrix_t2g)
    !
    u_jbasis(1,:) = (1.0/sqrt(6d0))*[ -one*sqrt(3d0),     one,          zero,  zero,           zero, -one*sqrt(2d0)]
    u_jbasis(2,:) = (1.0/sqrt(6d0))*[  -xi*sqrt(3d0),     -xi,          zero,  zero,           zero,   xi*sqrt(2d0)]
    u_jbasis(3,:) = (1.0/sqrt(6d0))*[           zero,    zero,          zero, 2*one, -one*sqrt(2d0),           zero]
    u_jbasis(4,:) = (1.0/sqrt(6d0))*[           zero,    zero, one*sqrt(3d0),  -one, -one*sqrt(2d0),           zero]
    u_jbasis(5,:) = (1.0/sqrt(6d0))*[           zero,    zero, -xi*sqrt(3d0),   -xi,  -xi*sqrt(2d0),           zero]  
    u_jbasis(6,:) = (1.0/sqrt(6d0))*[           zero,   2*one,          zero,  zero,           zero,  one*sqrt(2d0)] 
  endif


  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc_imp(Nso,Nso))
  allocate(Hloc_imp_nn(Nspin,Nspin,Norb,Norb))
  allocate(Hloc_latt(Nso,Nso))
  allocate(Hloc_latt_nn(Nspin,Nspin,Norb,Norb))
  allocate(Zmats(Nso,Nso))

  allocate(S0(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  allocate(Gtest(Lmats))

  Weiss=zero
  Weiss_=zero
  Gmats=zero
  Greal=zero
  Smats=zero
  Sreal=zero
  Hloc_imp=zero
  Hloc_imp_nn=zero
  Hloc_latt=zero
  Hloc_latt_nn=zero
  dens=zero
  Gtest=zero


  !Build Hk and Hloc for impurity and lattice
  call set_hamiltonian()

  !Setup Bath
  call set_bath()

  !Setup solver
  call ed_init_solver(bath)
  !

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
    iloop=iloop+1
    call start_loop(iloop,nloop,"DMFT-loop")
    !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
    call ed_solve(bath) 

    !Set self-energy matrix and the ones used to calculate Gloc
    call set_sigma()
    call get_Zmats()
    !Compute the local gfs on the imaginary axis:
    call get_glocal()

    !Get the Weiss field/Delta function to be fitted
    call get_weissfield()

    !Fit the new bath
    call fit_new_bath()


    !Check convergence (if required change chemical potential)     
    Gtest=zero
    do iorb=1,Norb
      Gtest=Gtest+Weiss(1,1,iorb,iorb,:)/Norb
    enddo
    converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
    
    !converged = check_convergence(Weiss(1,1,1,1,:),dmft_error,nsuccess,nloop,reset=.false.)
    
    call ed_get_dens(dens)
    if(nread/=0d0)call ed_search_chemical_potential(xmu,sum(dens),converged)
    
    call end_loop
  enddo

  call dmft_kinetic_energy(Hk,Smats)

  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)

  call finalize_MPI()

  contains

    !-------------------------------------------------------------------------------------------
    !PURPOSE:  Hk model for the 2d square lattice
    !-------------------------------------------------------------------------------------------
    function hk_model(kpoint,N) result(hk)
      real(8),dimension(:) :: kpoint
      integer              :: N,ih
      real(8)              :: kx,ky
      complex(8)           :: hk(N,N)
      complex(8)           :: hk_t2g(6,6)
      kx=kpoint(1)
      ky=kpoint(2)
      Hk = zero
      hk_t2g = zero
      !
      !the base is ordered as [[yz, xz,xy]_up, [yz, xz,xy]_dw]
      !plus-minus x
      !xz,xz
      Hk_t2g(2,2) = Hk_t2g(2,2) - t2*exp(-xi*dot_product(kpoint,[1d0,0d0])) !right
      Hk_t2g(2,2) = Hk_t2g(2,2) - t2*exp(-xi*dot_product(kpoint,[-1d0,0d0])) !left
      !yz,yz
      Hk_t2g(1,1) = Hk_t2g(1,1) - t3*exp(-xi*dot_product(kpoint,[1d0,0d0])) !right
      Hk_t2g(1,1) = Hk_t2g(1,1) - t3*exp(-xi*dot_product(kpoint,[-1d0,0d0])) !left
      !yz,yz
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[1d0,0d0])) !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[-1d0,0d0])) !left
      !xz,yz
      Hk_t2g(1,2) = 0d0 !right and left
      Hk_t2g(2,1) = 0d0 !right and left
      !plus-minus y
      !xz,xz
      Hk_t2g(2,2) = Hk_t2g(2,2) - t3*exp(-xi*dot_product(kpoint,[0d0,1d0])) !right
      Hk_t2g(2,2) = Hk_t2g(2,2) - t3*exp(-xi*dot_product(kpoint,[0d0,-1d0])) !left
      !yz,yz
      Hk_t2g(1,1) = Hk_t2g(1,1) - t2*exp(-xi*dot_product(kpoint,[0d0,1d0])) !right
      Hk_t2g(1,1) = Hk_t2g(1,1) - t2*exp(-xi*dot_product(kpoint,[0d0,-1d0])) !left
      !yz,yz
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[0d0,1d0])) !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[0d0,-1d0])) !left
      !xz,yz
      Hk_t2g(1,2) = 0d0
      Hk_t2g(2,1) = 0d0
      !plusx plusy
      !xz,xz
      Hk_t2g(2,2) = 0d0
      Hk_t2g(2,2) = 0d0
      !yz,yz
      Hk_t2g(1,1) = 0d0
      Hk_t2g(1,1) = 0d0
      !yz,yz
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[1d0,1d0])) !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[-1d0,-1d0])) !left
      !xz,yz
      Hk_t2g(1,2) = Hk_t2g(1,2) - t4*exp(-xi*dot_product(kpoint,[1d0,1d0])) !right
      Hk_t2g(1,2) = Hk_t2g(1,2) - t4*exp(-xi*dot_product(kpoint,[-1d0,-1d0])) !left
      Hk_t2g(2,1) = Hk_t2g(2,1) - t4*exp(-xi*dot_product(kpoint,[1d0,1d0])) !right
      Hk_t2g(2,1) = Hk_t2g(2,1) - t4*exp(-xi*dot_product(kpoint,[-1d0,-1d0])) !left
      !plusx minusy
      !xz,xz
      Hk_t2g(2,2) = 0d0
      Hk_t2g(2,2) = 0d0
      !yz,yz
      Hk_t2g(1,1) = 0d0
      Hk_t2g(1,1) = 0d0
      !yz,yz
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[1d0,-1d0])) !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[-1d0,1d0])) !left
      !xz,yz
      Hk_t2g(1,2) = Hk_t2g(1,2) - t4*exp(-xi*dot_product(kpoint,[1d0,-1d0])) !right
      Hk_t2g(1,2) = Hk_t2g(1,2) - t4*exp(-xi*dot_product(kpoint,[-1d0,1d0])) !left
      Hk_t2g(2,1) = Hk_t2g(2,1) - t4*exp(-xi*dot_product(kpoint,[1d0,-1d0])) !right
      Hk_t2g(2,1) = Hk_t2g(2,1) - t4*exp(-xi*dot_product(kpoint,[-1d0,1d0])) !left
      !
      !
      Hk_t2g(4:6,4:6) = Hk_t2g(1:3,1:3)
      !here add SOC!
      Hk_t2g = Hk_t2g + spinorbit_matrix_t2g
      !
      if(NORB==2)then
        Hk(1:2,1:2) = Hk_t2g(1:3:2,1:3:2)
        Hk(3:4,3:4) = Hk_t2g(4:6:2,4:6:2)
      else
        Hk = Hk_t2g      
      endif
      !
    end function hk_model


    !---------------------------------------------------------------------
    !PURPOSE: GET Z
    !---------------------------------------------------------------------

    subroutine get_zmats()
      integer :: ii
      complex(8),dimension(Nso,Nso) :: s0so,Uinv,tmpmat,s0nh
      !
      !
      Zmats=zero
      s0so=nn2so(s0)
      Uinv=u_jbasis
      call inv(Uinv)
      if(ZJBASIS)then
         s0so=matmul(s0so,u_jbasis)
         s0so=matmul(Uinv,s0so)
      endif
      !
      s0nh = (s0so-conjg(transpose(s0so)))/(2*xi)
      tmpmat = zeye(6)-s0nh/(pi/beta)
      !tmpmat = zeye(6) - imag(s0so)/(pi/beta)
      call inv(tmpmat)
      Zmats=real(tmpmat)
      !
      if(master)then
        unit=free_unit()
        open(unit,file="z.dat",position='append')
        write(unit,*)Zmats(1,1),Zmats(2,2),Zmats(3,3),Zmats(4,4),Zmats(5,5),Zmats(6,6)
        close(unit)
        print*,"print Zmats"
        call print_a_matrix(Zmats)
      endif             
    end subroutine get_zmats



    !+---------------------------------------------------------------------------+
    !PURPOSE : reshape 
    !+---------------------------------------------------------------------------+
    function so2nn(Hlso) result(Hnnn)
      complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hlso
      complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnnn
      integer                                     :: iorb,jorb
      integer                                     :: ispin,jspin
      integer                                     :: is,js
      Hnnn=zero
      do ispin=1,Nspin
        do jspin=1,Nspin
          do iorb=1,Norb
            do jorb=1,Norb
              is = iorb + (ispin-1)*Norb
              js = jorb + (jspin-1)*Norb
              Hnnn(ispin,jspin,iorb,jorb) = Hlso(is,js)
            enddo
          enddo
        enddo
      enddo
    end function so2nn

    function nn2so(Hnnn) result(Hlso)
      complex(8),dimension(Nspin,Nspin,Norb,Norb) :: Hnnn
      complex(8),dimension(Nspin*Norb,Nspin*Norb) :: Hlso
      integer                                     :: iorb,jorb
      integer                                     :: ispin,jspin
      integer                                     :: is,js
      Hlso=zero
      do ispin=1,Nspin
        do jspin=1,Nspin
          do iorb=1,Norb
            do jorb=1,Norb
              is = iorb + (ispin-1)*Norb
              js = jorb + (jspin-1)*Norb
              Hlso(is,js) = Hnnn(ispin,jspin,iorb,jorb)
            enddo
          enddo
        enddo
      enddo
    end function nn2so
    

    !+---------------------------------------------------------------------------+
    !set hamiltonian or dos
    !+---------------------------------------------------------------------------+
    
    subroutine set_hamiltonian()
      !
      call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
      Lk = Nx*Nx
      allocate(Hk(Nso,Nso,Lk),Wt(Lk))
      call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
      Wt = 1d0/Lk
      if(master)call TB_write_hk(Hk(:,:,:),"Hk2d.dat",Nlat=1,Nspin=Nspin,Norb=Norb,Nkvec=[Nx,Nx])
      !Set Hloc of the lattice
      Hloc_latt   = zero
      Hloc_latt = sum(Hk,dim=3)/Lk
      where(abs(Hloc_latt)<1.d-10) Hloc_latt=0d0
      Hloc_latt_nn=so2nn(Hloc_latt)
      !
      Hloc_imp = Hloc_latt !If I'm working with HK: SOC already included
      Hloc_imp_nn = so2nn(Hloc_imp)
      !
      call ed_set_hloc(Hloc_imp) !Set the Hamiltonian of the impurity problem as Hloc_imp
      !
    end subroutine set_hamiltonian
    
    !+---------------------------------------------------------------------------+
    !set bath
    !+---------------------------------------------------------------------------+
    subroutine set_bath()
      !
      if(bath_type =="replica" .or. bath_type =="general")then
        select case(parametrization)
          case("symmetric")
            call parametrization_symmetric()
        end select
        !
        if(bath_type =="replica")then
          call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
        else
          call ed_set_Hgeneral(Hsym_basis,lambdasym_vectors)
        endif
        !
        Nb=ed_get_bath_dimension(NSymMats)
      else
        Nb=ed_get_bath_dimension()
      endif
      !
      !
      allocate(bath(Nb))
      allocate(bath_(Nb))
      bath_ = zero
    end subroutine set_bath

    !
    subroutine parametrization_symmetric()
       allocate(tmpmat(Nso,Nso))
       NsymMats=2
        !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      allocate(Hsym_t2g(Nspin,Nspin,3,3,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !
      !
      !#symmetry 1
      Hsym_t2g(:,:,:,:,1) = so2nn(zeye(Nso))
      Hsym_t2g(:,:,:,:,2) = spinorbit_matrix_t2g_nn
      if(abs(lambda)>0.0) Hsym_t2g(:,:,:,:,2) = Hsym_t2g(:,:,:,:,2)/(-0.5*lambda)
      !
      if (NORB==2)then
        Hsym_basis = Hsym_t2g(:,:,1:3:2,1:3:2,:)
      else
        Hsym_basis = Hsym_t2g
      endif
      !
      !
      write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)
      !
      !
      do irepl=1,Nbath
        onsite = irepl -1 - (Nbath-1)/2d0          ![-(Nbath-1)/2:(Nbath-1)/2]
        if(Nbath>1) onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite        !Multiplies the suitable identity
      enddo
      !
      !random values around lambda for off-diagonal components
      call random_init(.true., .true.)
      call RANDOM_NUMBER(lambdasym_vectors(:,2:NSymMats))
      lambdasym_vectors(:,2:NSymMats) = (lambdasym_vectors(:,2:NSymMats) - 0.5d0)*0.1 + lambda
      !
      if(mod(Nbath,2)==0 .and. Nbath>2)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
      endif    

    end subroutine parametrization_symmetric



    !+---------------------------------------------------------------------------+
    !set sigma
    !+---------------------------------------------------------------------------+
    
    subroutine set_sigma()

      call ed_get_sigma(Smats,axis='mats')
      call ed_get_sigma(Sreal,axis='real')
      s0=Smats(:,:,:,:,1)


    call dmft_write_gf(Smats,"Smats",axis='mats',iprint=6)

    end subroutine set_sigma 
 
    !+---------------------------------------------------------------------------+
    !get glocal
    !+---------------------------------------------------------------------------+
    subroutine get_glocal()
      !
      call dmft_get_gloc(Hk,Gmats,Smats,axis='mats')
      call dmft_get_gloc(Hk,Greal,Sreal,axis='real')
      !
      call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=6)
      call dmft_write_gf(Greal,"Gloc",axis='real',iprint=6)    
      !
    end subroutine get_glocal
    
    !+---------------------------------------------------------------------------+
    !get weiss field
    !+---------------------------------------------------------------------------+
    
    subroutine get_weissfield()
      !
      if(cg_scheme=='delta')then
        call dmft_self_consistency(Gmats,Smats,Weiss,Hloc_imp_nn)
        call dmft_write_gf(Weiss,"Delta",axis='mats',iprint=6)
      else
        call dmft_self_consistency(Gmats,Smats,Weiss)
        call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=6)
      endif          
    end subroutine get_weissfield


    !+---------------------------------------------------------------------------+
    !fit new bath
    !+---------------------------------------------------------------------------+
    
    subroutine fit_new_bath()
    
    !If mixing Weiss, do it before
    if(mixG0)then
      if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
      Weiss_=Weiss
    endif
    
    select case(bath_type)
    case("normal")
      call ed_chi2_fitgf(Weiss,bath)
      !all ed_ph_symmetrize_bath(bath,save=.true.)
      call ed_spin_symmetrize_bath(bath,save=.true.)
      call ed_orb_symmetrize_bath(bath,save=.true.)
    case("hybrid")
      call ed_chi2_fitgf(Weiss,bath)
      !call ed_ph_symmetrize_bath(bath,save=.true.)
      !call ed_spin_symmetrize_bath(bath,save=.true.)
    case("replica")
      call ed_chi2_fitgf(Weiss,bath)
    case("general")
      call ed_chi2_fitgf(Weiss,bath)
    end select


    !If mixing bath, do it after
    if(.not.mixG0)then
      if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
      Bath_=Bath
    endif
    end subroutine fit_new_bath


    !+---------------------------------------------------------------------------+
    !print zmats
    !+---------------------------------------------------------------------------+
    
  subroutine print_a_matrix(mat)
    real(8),dimension(Nspin*Norb,Nspin*Norb)               :: mat
    integer                                                :: is,js,Nso,unit
    character(len=32)                                      :: fmt
    !
    unit=LOGfile
    !
    Nso = Nspin*Norb
    do is=1,Nso
       write(unit,"(12(A1,F8.4,A1,2x))")&
            ('(',mat(is,js),')',js =1,Nso)
    enddo
    write(unit,*)""
    !
  end subroutine print_a_matrix
    
  subroutine print_c_matrix(mat)
    complex(8),dimension(Nspin*Norb,Nspin*Norb)            :: mat
    integer                                                :: is,js,Nso,unit
    character(len=32)                                      :: fmt
    !
    unit=LOGfile
    !
    Nso = Nspin*Norb
    do is=1,Nso
       write(unit,"(20(A1,F8.4,A1,F8.4,A1,2x))")&
            ('(',real(mat(is,js)),',',imag(mat(is,js)),')',js =1,Nso)
    enddo
    write(unit,*)""
    !
  end subroutine print_c_matrix


end program soc_test_model





