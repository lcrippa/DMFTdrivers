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
  real(8)                                       :: wmixing,onsite,lambda
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
  complex(8),allocatable                        :: Hloc(:,:),tmpmat(:,:)
  complex(8),allocatable                        :: Hloc_nn(:,:,:,:)
  real(8),dimension(:),allocatable              :: H0 !elements on the diagonal of Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats,Smats_for_gloc
  complex(8),allocatable,dimension(:,:,:,:)     :: S0
  real(8),allocatable,dimension(:,:)            :: Zmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gkmats
  complex(8),allocatable,dimension(:)           :: Gtest
  !
  character(len=16)                             :: finput,foutput,model,PARAMETRIZATION
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
  !
  !
  call ed_read_input(trim(finput))
  !
  !
  !if model uses dos, use DOSFLAG=T
  if(model=="bethe" .or. model=="cosinedos")DOSFLAG=.true.
  !
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
  call parse_input_variable(Wbethe,"WBETHE",finput,default=(/(0.d0,iii=1,size(Wbethe))/),comment="Bethe half bandwidth")
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
    spinorbit_matrix_t2g(1,:) = (-1/2)*[  zero,   -xi,   zero,   zero,  zero,   one]
    spinorbit_matrix_t2g(2,:) = (-1/2)*[    xi,  zero,   zero,   zero,  zero,   -xi]
    spinorbit_matrix_t2g(3,:) = (-1/2)*[  zero,  zero,   zero,   -one,    xi,  zero]
    spinorbit_matrix_t2g(4,:) = (-1/2)*[  zero,  zero,   -one,   zero,    xi,  zero]
    spinorbit_matrix_t2g(5,:) = (-1/2)*[  zero,  zero,    -xi,    -xi,  zero,  zero]  
    spinorbit_matrix_t2g(6,:) = (-1/2)*[   one,   xi,    zero,   zero,  zero,  zero] 
    spinorbit_matrix_t2g_nn = so2nn(spinorbit_matrix_t2g)
    !
    u_jbasis(1,:) = (1/sqrt(6d0))*[ -one*sqrt(3d0),     one,          zero,  zero,           zero, -one*sqrt(2d0)]
    u_jbasis(2,:) = (1/sqrt(6d0))*[  -xi*sqrt(3d0),     -xi,          zero,  zero,           zero,   xi*sqrt(2d0)]
    u_jbasis(3,:) = (1/sqrt(6d0))*[           zero,    zero,          zero, 2*one, -one*sqrt(2d0),           zero]
    u_jbasis(4,:) = (1/sqrt(6d0))*[           zero,    zero, one*sqrt(3d0),  -one, -one*sqrt(2d0),           zero]
    u_jbasis(5,:) = (1/sqrt(6d0))*[           zero,    zero, -xi*sqrt(3d0),   -xi,  -xi*sqrt(2d0),           zero]  
    u_jbasis(6,:) = (1/sqrt(6d0))*[           zero,   2*one,          zero,  zero,           zero,  one*sqrt(2d0)] 
  endif


  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_for_gloc(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Hloc(Nso,Nso))
  allocate(Zmats(Nso,Nso))
  allocate(Hloc_nn(Nspin,Nspin,Norb,Norb))
  allocate(S0(Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  allocate(Gtest(Lmats))

  Weiss=zero
  Weiss_=zero
  Gmats=zero
  Greal=zero
  Smats=zero
  Smats_for_gloc=zero
  Hloc=zero
  Hloc_nn=zero
  dens=zero
  Gtest=zero


  !Build Hk and Hloc
  allocate(Ebands(Nso,Le));Ebands=zero
  allocate(Dbands(Nso,Le));Dbands=zero
  allocate(Wband(Nso));Wband=zero
  allocate(H0(Nso));H0=zero
  allocate(de(Nso))
  Wband = Wbethe(:Nso)  !half bandwidth
  H0    = Dbethe(:Nso)  !crystal field (centri delle dos)
  do iorb=1,Nso
    Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb)) !Ebands is the range of frequencies (energies)
    Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb) !Dbands is the dos of the given Ebands
  enddo
  Hloc=diag(H0)
  do iii=1,Nso
    unit=free_unit()
    open(unit,file="testdos_so"//reg(txtfy(iii))//".dat")
    do jjj=1,size(ebands,2)
      write(unit,*)ebands(iii,jjj),dbands(iii,jjj)
    enddo
    close(unit)
  enddo
  Hloc_nn=so2nn(Hloc)   


  !Setup Bath
  call set_bath()

  call ed_set_hloc(Hloc-(lambda/2)*spinorbit_matrix_t2g) !this sets Hloc and Hloc plus one-body terms that I want to be in the interaction part (the SOC)


  !Setup solver
  call ed_init_solver(bath)


  

  !Fit the bath
  if(.false.)then
    call read_weissdelta_matsubara(Weiss)
  else
    call read_sigma_matsubara(Smats)
    do iii=1,Lmats
      Smats_for_gloc(:,:,:,:,iii)=Smats(:,:,:,:,iii)-so2nn((lambda/2)*spinorbit_matrix_t2g)
    enddo
    call dmft_get_gloc(Ebands,Dbands,H0,Gmats,Smats_for_gloc,axis='mats',diagonal=.false.)
    call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=6)
    !
    if(cg_scheme=='delta')then
      call dmft_self_consistency(Gmats,Smats,Weiss,Hloc_nn-so2nn((lambda/2)*spinorbit_matrix_t2g))
      call dmft_write_gf(Weiss,"Delta",axis='mats',iprint=6)
    else
      call dmft_self_consistency(Gmats,Smats,Weiss)
      call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=6)       
    endif
  endif
   
  select case(bath_type)
  case("normal")
    call ed_chi2_fitgf(Weiss,bath)
    call ed_ph_symmetrize_bath(bath,save=.true.)
    call ed_spin_symmetrize_bath(bath,save=.true.)
  case("hybrid")
    call ed_chi2_fitgf(Weiss,bath)
    call ed_ph_symmetrize_bath(bath,save=.true.)
    call ed_spin_symmetrize_bath(bath,save=.true.)
  case("replica")
    call ed_chi2_fitgf(Weiss,bath)
  case("general")
    call ed_chi2_fitgf(Weiss,bath)
  end select

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
      if( .not. DOSFLAG)Hk = Hk + lambda*spinorbit_matrix_t2g
      !
    end function hk_model

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
    !set bath
    !+---------------------------------------------------------------------------+
    subroutine set_bath()
      !
      if(bath_type =="replica" .or. bath_type =="general")then
        select case(parametrization)
          case("diagonal")
            call parametrization_diagonal()
          case("general")
            call parametrization_general()
          case("jxjyjz_real")
            call parametrization_jxjyjz_real()
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
        onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite        !Multiplies the suitable identity
      enddo
      !
      if(mod(Nbath,2)==0 .and. Nbath>2)then
        lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
        lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
      endif   
    end subroutine parametrization_diagonal
    !
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
        onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
        lambdasym_vectors(irepl,1) = onsite        !Multiplies the suitable identity
      enddo
      !
      !random values around lambda for off-diagonal components
      call random_init(.true., .true.)
      call RANDOM_NUMBER(lambdasym_vectors(:,2:NSymMats:2))
      lambdasym_vectors(:,2:NSymMats:2) = (lambdasym_vectors(:,2:NSymMats:2) - 0.5d0)*1 + lambda
      lambdasym_vectors(:,1:NSymMats:2) = (lambdasym_vectors(:,1:NSymMats:2) - 0.5d0)*1 + lambda
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
        onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
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

  !-----------------------------------------------------------------------------!
  ! purpose: read quantities from disk
  !-----------------------------------------------------------------------------!

  subroutine read_weissdelta_matsubara(Weiss)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)      :: Weiss
    character(len=30)                                      :: suffix,filename
    integer                                                :: ispin,jspin,iorb,jorb
    !

    !
    if(master)then
      do ispin=1,Nspin
         do jspin=1,Nspin
            do iorb=1,Norb
                do jorb=1,Norb
                 suffix="_l"//reg(txtfy(iorb))//"m"//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.dat"
                 if(cg_scheme=="delta")then
                     call sread("../weiss_components/Delta"//trim(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))
                 else
                     call sread("../weiss_components/Delta"//trim(suffix),wm,Weiss(ispin,jspin,iorb,jorb,:))                 
                 endif
                enddo
            enddo
         enddo
      enddo
    endif

  end subroutine read_weissdelta_matsubara
  
    subroutine read_sigma_matsubara(Selfenergy)
    complex(8),dimension(Nspin,Nspin,Norb,Norb,Lmats)      :: Selfenergy
    character(len=30)                                      :: suffix
    integer                                                :: ispin,jspin,iorb,jorb
    !
    if(master)then
       do ispin=1,Nspin
        do jspin=1,Nspin
          do iorb=1,Norb
              do jorb=1,Norb
               suffix="_l"//reg(txtfy(iorb))//reg(txtfy(jorb))//"_s"//reg(txtfy(ispin))//reg(txtfy(jspin))//"_iw.ed"
               call sread("../sigma_components/impSigma"//trim(suffix),wm,Selfenergy(ispin,jspin,iorb,jorb,:))
              enddo
            enddo
          enddo
       enddo
    endif
  end subroutine read_sigma_matsubara

end program soc_test_model





