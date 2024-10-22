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
  complex(8),allocatable                        :: Hloc_imp(:,:),Hloc_latt(:,:),tmpmat(:,:)
  complex(8),allocatable                        :: Hloc_imp_nn(:,:,:,:),Hloc_latt_nn(:,:,:,:)
  real(8),dimension(:),allocatable              :: H0 !elements on the diagonal of Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
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
  character(len=16)                             :: finput,foutput,model,PARAMETRIZATION
  complex(8),allocatable                        :: Hk(:,:,:)
  real(8),allocatable                           :: Wt(:)
  !
  integer                                       :: comm,rank,unit
  logical                                       :: master
  logical                                       :: betheSC,mixG0,symOrbs,DOSFLAG,ZJBASIS
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  
  call parse_cmd_variable(finput,"FINPUT",default='input.in')
  call ed_read_input(trim(finput))
  !
  !
  Nso=Nspin*Norb
  !
  !
  ! 
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)
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


  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc_imp(Nso,Nso))
  allocate(Hloc_imp_nn(Nspin,Nspin,Norb,Norb))


  Hloc_imp=zero
  call ed_set_hloc(Hloc_imp)

  !Setup Bath
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb))
  bath = zero

  !Setup solver
  call ed_init_solver(bath)

  !Solve once
  call ed_solve(bath,sflag=.false.) 

end program soc_test_model





