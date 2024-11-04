program soc_test_model
  USE SCIFOR
  USE DMFT_TOOLS
  USE EDIPACK2 
  USE MPI
  implicit none

  integer                                       :: iloop,Nb,Lk,Nx,Nso,ik,iorb,irepl,le,iii,jjj,NsymMats
  logical                                       :: converged
  real(8),dimension(:),allocatable              :: Wbethe,Dbethe
  real(8),dimension(:,:),allocatable            :: Dbands
  real(8),dimension(:,:),allocatable            :: Ebands
  real(8)                                       :: wmixing,onsite,lambda,orb_shift
  real(8),dimension(:),allocatable              :: ts,Dband
  real(8),dimension(:),allocatable              :: de,dens
  real(8),dimension(:),allocatable              :: Wband,wfreq,eps_array,wm,wr
  !Bath:
  real(8),allocatable                           :: Bath(:),Bath_(:)
  real(8),dimension(:,:),allocatable            :: lambdasym_vectors
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis
  complex(8),dimension(:,:,:,:),allocatable     :: spinorbit_matrix_d_nn,spinorbit_matrix_t2g_nn
  complex(8),dimension(:,:),allocatable         :: spinorbit_matrix_d,spinorbit_matrix_t2g, u_jbasis
  !The local hybridization function:
  complex(8),allocatable                        :: Hloc_imp(:,:),Hloc_latt(:,:),tmpmat(:,:)
  complex(8),allocatable                        :: Hloc_imp_nn(:,:,:,:),Hloc_latt_nn(:,:,:,:)
  real(8),dimension(:),allocatable              :: H0 !elements on the diagonal of Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats,OneoverGmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats,Smats_for_gloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Sreal,Sreal_for_gloc
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
  logical                                       :: betheSC,mixG0,symOrbs,DOSFLAG,ZJBASIS


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  

  call parse_cmd_variable(finput,"FINPUT",default='input.in')
  call parse_input_variable(model,"MODEL",finput,default='bethe')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(betheSC,"BETHESC",finput,default=.false.)
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
  !if model uses dos, use DOSFLAG=T
  if(model=="bethe" .or. model=="cosinedos")DOSFLAG=.true.
  !
  print*,"DOSFLAG IS ",DOSFLAG
  !
  Nso=Nspin*Norb
  allocate(TS(Nso))
  allocate(Dband(Nso))
  allocate(Dbethe(Nso))
  allocate(Wbethe(Nso))
  
  
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)
  

  !
  call parse_input_variable(ts,"TS",finput,default=(/(0.5d0,iii=1,size(TS))/),comment="hopping parameter")
  call parse_input_variable(Dband,"Dband",finput,default=(/(0.d0,iii=1,size(Dband))/),comment="cystal field splittig (bands shift)")
  call parse_input_variable(Wbethe,"WBETHE",finput,default=(/(1.d0,iii=1,size(Wbethe))/),comment="Bethe half bandwidth")
  call parse_input_variable(Dbethe,"DBETHE",finput,default=(/(0.d0,iii=1,size(Dbethe))/),comment="Bethe cystal field splittig (dos center)")

  if (Nspin>1)then
    do iii =1,Norb
      TS(iii+Norb)=TS(iii)
      Dband(iii+Norb)=Dband(iii)
      Dbethe(iii+Norb)=DBETHE(iii)
      Wbethe(iii+Norb)=Wbethe(iii)
    enddo
  endif
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
  allocate(OneoverGmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats_for_gloc(Nspin,Nspin,Norb,Norb,Lmats),Sreal_for_gloc(Nspin,Nspin,Norb,Norb,Lreal))
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
  OneoverGmats=zero
  Greal=zero
  Smats=zero
  Smats_for_gloc=zero
  Sreal=zero
  Sreal_for_gloc=zero
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
    !call test_weissfield()
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

  !Get kinetic energy:
  if(DOSFLAG)then
    select case(model)
      case("cosine")
        call dmft_kinetic_energy(Hk,Smats_for_gloc)
      case("cosinedos")
        call dmft_kinetic_energy(Ebands,Dbands,H0,Smats_for_gloc)
      case("bethe")
        call dmft_kinetic_energy(Ebands,Dbands,H0,Smats_for_gloc)
    end select
  endif

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
      kx=kpoint(1)
      ky=kpoint(2)
      Hk = zero
      do ih=1,N
        Hk(ih,ih) = -one*2d0*ts(ih)*(cos(kx)+cos(ky)) + Dband(ih)
      enddo
      !
      !here add SOC!
      if( .not. DOSFLAG)Hk = Hk + spinorbit_matrix_t2g
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
    !PURPOSE : square lattice H to dos
    !+---------------------------------------------------------------------------+
      
    subroutine square_lattice_h_to_dos(hk)
      complex(8),dimension(:,:,:)              :: hk
      real(8)                                  :: val
      integer                                  :: is,iw,ik
      !
      allocate(Ebands(Nso,Le));Ebands=zero
      allocate(Dbands(Nso,Le));Dbands=zero
      allocate(de(Nso))
      allocate(H0(Nso));H0=zero
      !
      do is=1,size(Hk,1)
        Ebands(is,:) = linspace(-4d0*ts(is) + Dband(is),4d0*ts(is) + Dband(is),Le,mesh=de(is)) !Ebands is the range of frequencies (energies)
        H0(is) = 0d0 !here i have to set it to 0, not Dband(is), because the center of the DOS is encoded in the line above
      enddo
      !
      do ik=1,size(Hk,3)
        do is=1,size(Hk,1)
          val=Hk(is,is,ik)
          do iw=1,size(ebands,2)
            if(ebands(is,iw)<val.AND.ebands(is,iw+1)>val)then
              dbands(is,iw)=dbands(is,iw)+1
            endif
          enddo
        enddo
      enddo
      !this normalizes to 1
      do is=1,size(Hk,1)
        dbands(is,:)=dbands(is,:)/(de(is)*size(Hk,3))
      enddo
      !
      !multiply again by de because the code wants so, cfr invocation of dens_bethe in text function
      do is=1,size(Hk,1)
        dbands(is,:)=dbands(is,:)*de(is)
      enddo
    end subroutine square_lattice_h_to_dos


    !+---------------------------------------------------------------------------+
    !set hamiltonian or dos
    !+---------------------------------------------------------------------------+
    
    subroutine set_hamiltonian()
      !
      select case(model)
        case("cosine")
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

        case("cosinedos")
          call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
          Lk = Nx*Nx
          allocate(Hk(Nso,Nso,Lk),Wt(Lk))
          call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
          if(master)call TB_write_hk(Hk(:,:,:),"Hk2d.dat",Nlat=1,Nspin=Nspin,Norb=Norb,Nkvec=[Nx,Nx])
          call square_lattice_h_to_dos(hk) 
          do iii=1,Nso
            unit=free_unit()
            open(unit,file="testdos_so"//reg(txtfy(iii))//".dat")
            do jjj=1,size(ebands,2)
              write(unit,*)ebands(iii,jjj),dbands(iii,jjj)/(de(iii))
            enddo
            close(unit)
          enddo
          !Set Hloc of the lattice
          Hloc_latt=diag(H0)
          Hloc_latt_nn=so2nn(Hloc_latt)
          
        case("bethe")
          allocate(Ebands(Nso,Le));Ebands=zero
          allocate(Dbands(Nso,Le));Dbands=zero
          allocate(Wband(Nso));Wband=zero
          allocate(H0(Nso));H0=zero
          allocate(de(Nso))
          Wband = Wbethe(:Nso)  !half bandwidth
          H0    = Dbethe(:Nso)  !crystal field (centri delle dos)
          !
          H0(Norb) = H0(Norb) + orb_shift  !shift xy
          if(Nspin>1) H0(Nso) = H0(Nso) + orb_shift  !shift xy
          do iorb=1,Nso
            Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb)) !Ebands is the range of frequencies (energies)
            Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb) !Dbands is the dos of the given Ebands
          enddo
          do iii=1,Nso
            unit=free_unit()
            open(unit,file="testdos_so"//reg(txtfy(iii))//".dat")
            do jjj=1,size(ebands,2)
              write(unit,*)ebands(iii,jjj),dbands(iii,jjj)
            enddo
            close(unit)
          enddo
          !Set Hloc of the lattice
          Hloc_latt=diag(H0)
          Hloc_latt_nn=so2nn(Hloc_latt)    
      end select
         
      !
      !Set Hloc of the impurity
      if(DOSFLAG)then
        Hloc_imp = Hloc_latt + spinorbit_matrix_t2g  !If I'm working with DOS: we need to add local SOC to hloc
        Hloc_imp_nn = so2nn(Hloc_imp)
      else
        Hloc_imp = Hloc_latt !If I'm working with HK: SOC already included
        Hloc_imp_nn = so2nn(Hloc_imp)
      endif
    
      call ed_set_hloc(Hloc_imp) !Set the Hamiltonian of the impurity problem as Hloc_imp
      
    end subroutine set_hamiltonian
    
    !+---------------------------------------------------------------------------+
    !set bath
    !+---------------------------------------------------------------------------+
    subroutine set_bath()
      !
      if(bath_type =="replica" .or. bath_type =="general")then
        select case(parametrization)
          case("diagonal")
            call parametrization_diagonal()
          case("orb_shifted")
            call parametrization_orb_shifted()
          case("general")
            call parametrization_general()
          case("jxjyjz_real")
            call parametrization_jxjyjz_real()
          case("symmetric")
            call parametrization_symmetric()
          case("orb_shifted_symmetric")
            call parametrization_orbshift_symmetric()
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

    subroutine parametrization_diagonal()
      NSymMats=1
      !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero !SOC matrix is Hermitian and has 12 nonzero entries. 6 of them are independent and in principle complex, so 12 parameters. +1 which is the identity
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !#symmetry 1
      Hsym_basis(:,:,:,:,1) = so2nn(zeye(Nso)) !identity
      !
      write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)
      !
      do irepl=1,Nbath
        onsite = irepl -1 - (Nbath-1)/2d0          ![-(Nbath-1)/2:(Nbath-1)/2]
        if(Nbath>1) onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite        !Multiplies the suitable identity
      enddo
      !
      if(mod(Nbath,2)==0 .and. Nbath>2)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
      endif   
    end subroutine parametrization_diagonal
    !
    subroutine parametrization_orb_shifted()
      NSymMats=2
      !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero !SOC matrix is Hermitian and has 12 nonzero entries. 6 of them are independent and in principle complex, so 12 parameters. +1 which is the identity
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !#symmetry 1
      tmpmat=zeye(Nso)
      tmpmat(3,3)=0.d0
      tmpmat(6,6)=0.d0
      Hsym_basis(:,:,:,:,1) = so2nn(tmpmat) !first two orbitals 
      !
      tmpmat=zeye(Nso)
      tmpmat(1,1)=0.d0
      tmpmat(2,2)=0.d0
      tmpmat(4,4)=0.d0
      tmpmat(5,5)=0.d0
      Hsym_basis(:,:,:,:,2) = so2nn(tmpmat) !third orbital
      !
      write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)
      !
      do irepl=1,Nbath
        onsite = irepl -1 - (Nbath-1)/2d0          ![-(Nbath-1)/2:(Nbath-1)/2]
        if(Nbath>1) onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite              !Multiplies the suitable identity
        lambdasym_vectors(irepl,2) = onsite + orb_shift  !Multiplies the suitable identity
      enddo
      !
      if(mod(Nbath,2)==0 .and. Nbath>2)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
      endif   
    end subroutine parametrization_orb_shifted
    !
    subroutine parametrization_general()
      allocate(tmpmat(Nso,Nso))
      NsymMats=13
      !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero !SOC matrix is Hermitian and has 12 nonzero entries. 6 of them are independent and in principle complex, so 12 parameters. +1 which is the identity
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !
      !
      !#symmetry 1
      Hsym_basis(:,:,:,:,1) = so2nn(zeye(Nso)) !identity
      !
      !#symmetry 2
      tmpmat=zero
      tmpmat(1,2)=1
      tmpmat(2,1)=1
      Hsym_basis(:,:,:,:,2) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 3
      tmpmat=zero
      tmpmat(1,2)=xi
      tmpmat(2,1)=-xi
      Hsym_basis(:,:,:,:,3) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 4
      tmpmat=zero
      tmpmat(1,6)=1
      tmpmat(6,1)=1
      Hsym_basis(:,:,:,:,4) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 5
      tmpmat=zero
      tmpmat(1,6)=xi
      tmpmat(6,1)=-xi
      Hsym_basis(:,:,:,:,5) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 6
      tmpmat=zero
      tmpmat(2,6)=1
      tmpmat(6,2)=1
      Hsym_basis(:,:,:,:,6) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 7
      tmpmat=zero
      tmpmat(2,6)=xi
      tmpmat(6,2)=-xi
      Hsym_basis(:,:,:,:,7) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 8
      tmpmat=zero
      tmpmat(3,4)=1
      tmpmat(4,3)=1
      Hsym_basis(:,:,:,:,8) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 9
      tmpmat=zero
      tmpmat(3,4)=xi
      tmpmat(4,3)=-xi
      Hsym_basis(:,:,:,:,9) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 10
      tmpmat=zero
      tmpmat(3,5)=1
      tmpmat(5,3)=1
      Hsym_basis(:,:,:,:,10) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 11
      tmpmat=zero
      tmpmat(3,5)=xi
      tmpmat(5,3)=-xi
      Hsym_basis(:,:,:,:,11) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 12
      tmpmat=zero
      tmpmat(4,5)=1
      tmpmat(5,4)=1
      Hsym_basis(:,:,:,:,12) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 13
      tmpmat=zero
      tmpmat(4,5)=xi
      tmpmat(5,4)=-xi
      Hsym_basis(:,:,:,:,13) = -0.5*so2nn(tmpmat)
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
      !call RANDOM_NUMBER(lambdasym_vectors(:,2:NSymMats:2))
      !lambdasym_vectors(:,2:NSymMats:2) = (lambdasym_vectors(:,2:NSymMats:2) - 0.5d0)*0.1 + lambda
      call RANDOM_NUMBER(lambdasym_vectors(:,2:NSymMats))
      lambdasym_vectors(:,2:NSymMats) = (lambdasym_vectors(:,2:NSymMats) - 0.5d0)*0.1 + lambda
      !
      if(mod(Nbath,2)==0 .and. Nbath>2)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
      endif    
    !
    end subroutine parametrization_general
    !
    !
    !
    subroutine parametrization_jxjyjz_real()
       allocate(tmpmat(Nso,Nso))
       NsymMats=4
        !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !
      !
      !#symmetry 1
      Hsym_basis(:,:,:,:,1) = so2nn(zeye(Nso)) !identity
      !
      !#symmetry 2
      tmpmat=zero
      tmpmat(2,6)=xi
      tmpmat(6,2)=-xi
      tmpmat(3,5)=-xi
      tmpmat(5,3)=xi
      Hsym_basis(:,:,:,:,2) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 3
      tmpmat=zero
      tmpmat(1,6)=-1
      tmpmat(6,1)=-1
      tmpmat(3,4)=1
      tmpmat(4,3)=1
      Hsym_basis(:,:,:,:,3) = -0.5*so2nn(tmpmat)
      !
      !#symmetry 4
      tmpmat=zero
      tmpmat(1,2)=xi
      tmpmat(2,1)=-xi
      tmpmat(4,5)=-xi
      tmpmat(5,4)=xi
      Hsym_basis(:,:,:,:,4) = -0.5*so2nn(tmpmat)
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
    !    
    end subroutine parametrization_jxjyjz_real
    !
    subroutine parametrization_symmetric()
       allocate(tmpmat(Nso,Nso))
       NsymMats=2
        !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !
      !
      !#symmetry 1
      Hsym_basis(:,:,:,:,1) = so2nn(zeye(Nso))
      Hsym_basis(:,:,:,:,2) = spinorbit_matrix_t2g_nn
      if(abs(lambda)>0.0) Hsym_basis(:,:,:,:,2) = Hsym_basis(:,:,:,:,2)/(-0.5*lambda)
      !
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

    subroutine parametrization_orbshift_symmetric()
       allocate(tmpmat(Nso,Nso))
       NsymMats=3
        !
      allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
      !
      !define all the matrices
      !
      !
      !
      !define all the matrices
      !#symmetry 1
      tmpmat=zeye(Nso)
      tmpmat(3,3)=0.d0
      tmpmat(6,6)=0.d0
      Hsym_basis(:,:,:,:,1) = so2nn(tmpmat) !first two orbitals 
      !#symmetry 2
      tmpmat=zeye(Nso)
      tmpmat(1,1)=0.d0
      tmpmat(2,2)=0.d0
      tmpmat(4,4)=0.d0
      tmpmat(5,5)=0.d0
      Hsym_basis(:,:,:,:,2) = so2nn(tmpmat) !third orbital
      !#symmetry 3
      Hsym_basis(:,:,:,:,3) = spinorbit_matrix_t2g_nn
      if(abs(lambda)>0.0) Hsym_basis(:,:,:,:,3) = Hsym_basis(:,:,:,:,3)/(-0.5*lambda)
      !
      write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)
      !
      do irepl=1,Nbath
        onsite = irepl -1 - (Nbath-1)/2d0          ![-(Nbath-1)/2:(Nbath-1)/2]
        if(Nbath>1) onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite              !Multiplies the suitable identity
        lambdasym_vectors(irepl,2) = onsite + orb_shift  !Multiplies the suitable identity
      enddo
      !
      !random values around lambda for off-diagonal components
      call random_init(.true., .true.)
      call RANDOM_NUMBER(lambdasym_vectors(:,3:NSymMats))
      lambdasym_vectors(:,3:NSymMats) = (lambdasym_vectors(:,3:NSymMats) - 0.5d0)*0.1 + lambda
      !
      if(mod(Nbath,2)==0 .and. Nbath>2)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
      endif    

    end subroutine parametrization_orbshift_symmetric

    !+---------------------------------------------------------------------------+
    !set sigma
    !+---------------------------------------------------------------------------+
    
    subroutine set_sigma()

      call ed_get_sigma(Smats,axis='mats')
      call ed_get_sigma(Sreal,axis='real')
      s0=Smats(:,:,:,:,1)
      !DOSFLAG: use the dos without the Hlambda, put Hlambda in the Sigma
      
      if(DOSFLAG)then
        do iii=1,Lmats
          Smats_for_gloc(:,:,:,:,iii)=Smats(:,:,:,:,iii) + so2nn(spinorbit_matrix_t2g) !Add Hlambda for Gloc calculation
        enddo
        do iii=1,Lreal
          Sreal_for_gloc(:,:,:,:,iii)=Sreal(:,:,:,:,iii) + so2nn(spinorbit_matrix_t2g) !Add Hlambda for Gloc calculation
        enddo
      else
        do iii=1,Lmats
          Smats_for_gloc(:,:,:,:,iii)=Smats(:,:,:,:,iii)
        enddo
        do iii=1,Lreal
          Sreal_for_gloc(:,:,:,:,iii)=Sreal(:,:,:,:,iii)
        enddo
      endif
      call dmft_write_gf(Smats_for_gloc,"Smats_for_gloc",axis='mats',iprint=6)

    end subroutine set_sigma 
 
    !+---------------------------------------------------------------------------+
    !get glocal
    !+---------------------------------------------------------------------------+
    subroutine get_glocal()
      !
      allocate(tmpmat_so(Nso,Nso))
      if(DOSFLAG)then
        print*,"Outside the function:"
        print*,"size(ebands,1)=",size(ebands,1)
        call dmft_get_gloc(Ebands,Dbands,H0,Gmats,Smats_for_gloc,axis='mats',diagonal=.false.)
        !call dmft_get_gloc(Ebands,Dbands,H0,Greal,Sreal_for_gloc,axis='real',diagonal=.false.)
      else
        call dmft_get_gloc(Hk,Gmats,Smats_for_gloc,axis='mats')
        !call dmft_get_gloc(Hk,Greal,Sreal_for_gloc,axis='real')
      endif
      !
      call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=6)
      !call dmft_write_gf(Greal,"Gloc",axis='real',iprint=6)    
      !
      deallocate(tmpmat_so)
    end subroutine get_glocal
    
    !+---------------------------------------------------------------------------+
    !get weiss field
    !+---------------------------------------------------------------------------+
    
    subroutine get_weissfield()
      !
      if(cg_scheme=='delta')then
        if(model=='bethe' .and. betheSC)then
          call dmft_get_delta_normal_bethe(Gmats,Weiss,Hloc_imp_nn,Wband)
        else
          call dmft_self_consistency(Gmats,Smats,Weiss,Hloc_imp_nn)
        endif
        call dmft_write_gf(Weiss,"Delta",axis='mats',iprint=6)
      else
        if(model=='bethe' .and. betheSC)then
          call dmft_get_weiss_normal_bethe(Gmats,Weiss,Hloc_imp_nn,Wband)
        else
          call dmft_self_consistency(Gmats,Smats,Weiss)
        endif
        call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=6)
      endif          
    end subroutine get_weissfield


    !Test function to see that Weiss and Delta are consistent

    subroutine test_weissfield()
      !
      allocate(tmpmat_so(Nso,Nso))
      !Get Weiss from library (the object is called Weiss_)
      call dmft_self_consistency(Gmats,Smats,Weiss_)
      call dmft_write_gf(Weiss_,"Weiss",axis='mats',iprint=6)
      Weiss_=zero
      !Get Delta from library (the object is called Weiss)
      call dmft_self_consistency(Gmats,Smats,Weiss,Hloc_imp_nn)
      call dmft_write_gf(Weiss,"Delta",axis='mats',iprint=6)
      !Get Weiss from here
      do iii=1,Lmats
        tmpmat_so=xi*wm(iii)*eye(Nso)-Hloc_imp-nn2so(Weiss(:,:,:,:,iii))
        call inv(tmpmat_so)
        Weiss_(:,:,:,:,iii)=so2nn(tmpmat_so)
      enddo
      !Print Weiss from here
      call dmft_write_gf(Weiss_,"Weiss_fromhere",axis='mats',iprint=6)
      deallocate(tmpmat_so)
      !
    end subroutine test_weissfield
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




  subroutine dmft_get_weiss_normal_bethe(Gloc,Weiss,Hloc,Wbands)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Weiss ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    real(8),dimension(:),intent(in)               :: Wbands ![Nspin*Norb]
    !aux
    complex(8),dimension(size(Gloc,5))            :: invWeiss ![Lmats]
    integer                                       :: Nspin,Norb,Nso,Lmats
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !
    !
    if(master)then
       !Testing part:
       Nspin = size(Gloc,1)
       Norb  = size(Gloc,3)
       Lmats = size(Gloc,5)
       Nso   = Nspin*Norb
       !
       !
       !\calG0^{-1}_aa = iw + mu - H_0 - d**2/4*Gmats
       Weiss=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             invWeiss = (xi*wm(:)+xmu) - Hloc(ispin,ispin,iorb,iorb) - &
                  0.25d0*Wbands(iorb+(ispin-1)*Norb)**2*Gloc(ispin,ispin,iorb,iorb,:)
             Weiss(ispin,ispin,iorb,iorb,:) = one/invWeiss
          enddo
       enddo
       !
       !
    endif
  end subroutine dmft_get_weiss_normal_bethe
  
  
  subroutine dmft_get_delta_normal_bethe(Gloc,Delta,Hloc,Wbands)
    complex(8),dimension(:,:,:,:,:),intent(in)    :: Gloc  ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:,:),intent(inout) :: Delta ! [Nspin][Nspin][Norb][Norb][Lmats]
    complex(8),dimension(:,:,:,:),intent(in)      :: Hloc  ! [Nspin][Nspin][Norb][Norb]
    real(8),dimension(:),intent(in)               :: Wbands ![Nspin*Norb]
    !aux
    integer                                       :: Nspin,Norb,Nso,Lmats
    integer                                       :: i,iorb,jorb,ispin,jspin,io,jo
    !
    !
    !
    if(master)then
       !Testing part:
       Nspin = size(Gloc,1)
       Norb  = size(Gloc,3)
       Lmats = size(Gloc,5)
       Nso   = Nspin*Norb
       !
       !
       !\calG0^{-1}_aa = d**2/4*Gmats
       Delta=zero
       do ispin=1,Nspin
          do iorb=1,Norb
             Delta(ispin,ispin,iorb,iorb,:) = 0.25d0*Wbands(iorb+(ispin-1)*Norb)**2*Gloc(ispin,ispin,iorb,iorb,:)
          enddo
       enddo
       !
       !
    endif
  end subroutine dmft_get_delta_normal_bethe







end program soc_test_model





