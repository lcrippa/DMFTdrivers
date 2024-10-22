program ed_bilayer
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                     :: iloop,Lk,Nso,Nlso,Nlat
  logical                                     :: converged,conv_dens,converged0
  integer :: ilat
  !Bath:
  integer                                     :: Nb
  real(8),allocatable                         :: Bath(:,:),Bath_prev(:,:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:,:)
  complex(8),allocatable                      :: Sreal(:,:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:,:)
  complex(8),allocatable                      :: Greal(:,:,:,:,:,:)
  complex(8),allocatable,dimension(:,:)         :: Gtest
  real(8),allocatable,dimension(:)            :: dens,dens_prev,docc
  !hamiltonian input:
  !variables for the model:
  complex(8),allocatable                      :: Hk(:,:,:),Hk_loop(:,:,:)
  complex(8),allocatable                      :: Hloc(:,:,:,:,:)
  complex(8),allocatable                      :: SigmaHk(:,:),Zmats(:,:)
  integer                                     :: Nk,Nkpath
  real(8)                                     :: ts,tperp,Vel,lambda,wmixing,ntarget,wmixing_dens
  character(len=16)                           :: finput
  !MPI Vars:
  integer                                     :: comm,rank,mpierr,mpiSize
  logical                                     :: master
  !
  !
  real(8),dimension(2)                        :: e1,e2   !real-space lattice basis
  real(8),dimension(2)                        :: bk1,bk2 !reciprocal space lattice basis
  real(8),dimension(2)                        :: bklen
  !
  real(8)                                     :: Vnn  !+- nearest neighbourg interaction
  real(8)                                     :: ntop,nbot,dens_check
  integer                                     :: uio,uio_loop,imu
  real(8)                                     :: ntest1,ntest2,xmu_imp,xmu1,xmu2,alpha_mu,dens_error
  logical                                     :: fix_mu,flag_mpi,printG
  logical                                     :: fixing_newton
  integer                                     :: test_strides
  !
  real(8),dimension(:,:),allocatable :: kgrid
  integer,dimension(:,:),allocatable :: ik2ij,ik_diff
  integer,dimension(:,:),allocatable :: ij2ik
  !
  !+- exciton-polaron problem
  real(8)              :: thop_mev,am_nm,thop_empty,gap_empty,beta_ancillary,dos_ancillary
  real(8)              :: wX,epsX,winiX,wfinX,mX,h0,delta_wgrid

  logical              :: empty_quadratic
  

  integer              :: LwX,Lx,iw,ik,jL,jR
  real(8),dimension(:),allocatable :: wrX,wr_tmp,wr
  complex(8)           :: wcmplx,DX_tmp
  complex(8),dimension(:),allocatable :: DX
  complex(8),dimension(:,:),allocatable :: DX_dressed

  character(len=6)     :: fix_mu_scheme           !newton/f_zero
  real(8),dimension(:,:),allocatable :: Aloc,Ak_test
    
  real(8),dimension(2)               :: pointK,pointKp,pointM,kdiff
  real(8),dimension(:,:),allocatable :: KPath

  complex(8),dimension(:,:,:),allocatable :: pi_irreducible,Lambda_vertex
  complex(8),dimension(:,:),allocatable :: Sigma_exciton
  !
  real(8),dimension(:,:),allocatable :: lambda_imag_tmp
  real(8) :: Vex
  real(8) :: dw1,dw2,dw3
  logical :: tmp_debug_vertex
  !
  !
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  mpiSize = get_size_MPI(comm)
  !

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BL.in')  
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=100)
  call parse_input_variable(ts,"TS",finput,default=1d0)
  call parse_input_variable(tperp,"TPERP",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(Vel,"Vel",finput,default=0.d0)

  call parse_input_variable(Vnn,"Vnn",finput,default=0.d0)
  call parse_input_variable(ntop,"ntop",finput,default=1.d0)
  call parse_input_variable(nbot,"nbot",finput,default=1.d0)
  !
  call parse_input_variable(alpha_mu,"alpha_mu",finput,default=0.5d0)
  call parse_input_variable(fix_mu,"fix_mu",finput,default=.false.)
  call parse_input_variable(flag_mpi,"flag_mpi",finput,default=.false.)
  call parse_input_variable(conv_dens,"conv_dens",finput,default=.false.)  
  call parse_input_variable(dens_error,"dens_error",finput,default=1.d-5)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.5d0)
  call parse_input_variable(wmixing_dens,"WMIXING_DENS",finput,default=0.5d0)
  call parse_input_variable(ntarget,"NTARGET",finput,default=2d0)
  call parse_input_variable(fix_mu_scheme,"fix_mu_scheme",finput,default='f_zero')
  call parse_input_variable(printG,"printG",finput,default=.true.)
  !
  call parse_input_variable(thop_mev,"thop_mev",finput,default=0.5d0,comment='hopping in meV')
  call parse_input_variable(thop_empty,"thop_empty",finput,default=2.d0,comment='thop in the the 2nd moire band; dimensionless')
  call parse_input_variable(gap_empty,"gap_empty",finput,default=30.d0,comment='gap btw 1st and 2nd moire bands; dimensionless')


  call parse_input_variable(empty_quadratic,"empty_quadratic",finput,default=.false.,comment='gap btw 1st and 2nd moire bands; dimensionless')
  call parse_input_variable(beta_Ancillary,"beta_ancillary",finput,default=1.d0,comment='smooth the ancillary dos')
  call parse_input_variable(dos_ancillary,"dos_ancillary",finput,default=0.1d0,comment='smooth the ancillary dos')

  call parse_input_variable(tmp_debug_vertex,"tmp_debug_vertex",finput,default=.false.,comment='tmp debug')

  !  
  call parse_input_variable(am_nm,"am_nm",finput,default=25d0,comment='moire periodicity in nm')
  !
  call parse_input_variable(wX,"wX",finput,default=1.6d0,comment='bare exciton energy in eV')

  call parse_input_variable(mX,"mX",finput,default=1.3d0,comment='bare exciton mass [units of me]')

  call parse_input_variable(epsX,"epsX",finput,default=1d-3,comment='exciton energy lifetime in eV')
  !
  call parse_input_variable(winiX,"winiX",finput,default=-5.d0,comment='minimum exciton energy in eV')
  call parse_input_variable(wfinX,"wfinX",finput,default= 5.d0,comment='maximum exciton energy in eV')
  call parse_input_variable(Vex,"Vex",finput,default=-1.d0,comment='e-X interaction')



  call parse_input_variable(delta_wgrid,"dwx",finput,default=0.05d0,comment='range for the finer grid eV')

  call parse_input_variable(LwX,"LwX",finput,default= 1000,comment='frequency discretizion for the exciton')

  !
  call parse_input_variable(test_strides,"test_strides",finput,default=0)

  !
  call ed_read_input(trim(finput),comm)
  !
  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")
  !
  !
  h0=7.62d-02 ! hbar^2/m [eV x nm^2]
  

  Nlat=2
  if(Nspin/=1.OR.Nlat/=2)stop "Wrong setup from input file: Nspin=1, Norb=2 -> 2Spin-Orbitals"
  Nso=Nspin*Norb
  Nlso=Nlat*Nso
  
  !Allocate Weiss Field:

  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gtest(Nlat,Lmats))
  allocate(dens(Norb)); allocate(dens_prev(Norb))
  allocate(docc(Nlat)); 
  allocate(SigmaHk(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  !Buil the Hamiltonian on a grid or on  path
  call set_SigmaHk()
  
  !+- a=1 is the hexagon edge
  !+- a' = sqrt(3.0)*a is the modulus of the primitive vector of the triangular lattice
  !+- e_1 = sqrt(3.0)/2.0 * a [-1, sqrt(3.0)]
  !+- e_2 = sqrt(3.0)/2.0 * a [ 1, sqrt(3.0)]
  e1 = sqrt(3d0)/2d0*[-1d0, sqrt(3d0)]
  e2 = sqrt(3d0)/2d0*[ 1d0, sqrt(3d0)]
  
  !RECIPROCAL LATTICE VECTORS:
  bklen=2d0*pi/3d0
  bk1=bklen*[ -sqrt(3d0), 1d0]
  bk2=bklen*[  sqrt(3d0), 1d0]
  call TB_set_bk(bkx=bk1,bky=bk2) 
  call build_hk_honeycomb()
  !
  

  !+- now do a single DMFT-loop
  if(nloop.gt.1) then
     if(master) write(*,*) 'nloop set to 1'
     nloop=1

     if(master) write(*,*) 'fix_mu set to F'
     fix_mu=.false.

  end if
  !
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nlat,Nb))
  allocate(Bath_prev(Nlat,Nb))
  call ed_init_solver(comm,bath)
  Bath_prev=Bath
  !DMFT loop
  iloop=0;converged=.false.;converged0=.false.
  dens=[ntop,nbot]
  dens_prev=dens
  uio=free_unit()
  if(master) then
     open(uio,file='init_ndens_hf.out')
     write(uio,'(3F18.10)') dens
     close(uio)
  end if
  !
  uio=free_unit()
  if(master) then
     open(uio,file='iloop_observables.out')
     close(uio)
  end if
  call set_Hloc(dens)    
  !
  !
  xmu1=0.d0
  xmu2=xmu1


  xmu_imp=xmu
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     
     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call set_Hloc(dens,'Hloc_set'//str(iloop,Npad=3)//'.dat')     
     
     !+- set fix_mu=F and skip all this part
     if(fix_mu) then
        !   bracketing the xmu  !

        select case(fix_mu_scheme)
        case('newton')
           call newton(get_delta_dens_imp_,xmu,eps=1d-4)
        case('f_zero')
           xmu1= -50
           xmu2=  50
           call fzero(get_delta_dens_imp,xmu1,xmu2,imu,tol_rel=1.d-8,tol_abs=1.d-6)
           xmu=xmu1
        case default           
           xmu1= -50
           xmu2=  50
           call fzero(get_delta_dens_imp,xmu1,xmu2,imu,tol_rel=1.d-8,tol_abs=1.d-6)
           xmu=xmu1
        end select
        !
     end if
     !
     call ed_solve(comm,bath,Hloc,mpi_lanc=flag_mpi)     
     call ed_get_sigma_matsubara(Smats,Nlat)
     call ed_get_sigma_realaxis(Sreal,Nlat)

     !
     !
     call ed_get_dens(dens,Nlat,iorb=1)
     call ed_get_docc(docc,Nlat,iorb=1)
     !
     uio=free_unit()
     if(master) then
        open(uio,file='iloop_observables.out',status='old',position='append')
        write(uio,'(10F18.10)') dens,docc,xmu
        close(uio)
     end if
     
     ! !Get GLOC:
     
     !+- here I should use Hk_loop !!!!
     if(.not.allocated(hk_loop)) then
        call mpi_barrier(comm,mpiERR)
        stop
     end if
     call dmft_gloc_matsubara(Hk_loop,Gmats,Smats)
     call dmft_gloc_realaxis(Hk_loop,Greal,Sreal)
     
     ! !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
     !
     if(printG) then
        call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=4,ineq_pad=2)
        call dmft_print_gf_realaxis(Greal,"Gloc",iprint=4,ineq_pad=2)
        call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=1,ineq_pad=2)
     end if               
     !
     call ed_chi2_fitgf(bath,Weiss,Hloc,ispin=1)  !+- perche qui si mangia anche Hloc?
     !
     !+- mix things
     Bath = wmixing*Bath + (1.d0-wmixing)*Bath_prev
     dens = wmixing_dens*dens + (1.d0-wmixing_dens)*dens_prev

     Gtest=Weiss(:,1,1,1,1,:)
     if(.not.conv_dens) then
        converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     else
        converged = check_convergence_local(dens,dens_error,nsuccess,nloop)
     end if
     !
     Bath_prev = Bath
     dens_prev = dens
     xmu_imp   = xmu
     !
     call end_loop
     !
  enddo
  !
  call set_Hloc(dens,'Hloc_last.dat')     
  !
  ! call dmft_gloc_realaxis(Hk_loop,Greal,Sreal)
  ! call dmft_kinetic_energy(Hk_loop,Smats)
  !


  !+- THE EXCITON-POLARON PROBLEM -+!
  allocate(wr(Lreal)); wr=0.d0
  wr=linspace(wini,wfin,Lreal) 

  !+- test Akw DONE
  allocate(Aloc(2,Lreal)); Aloc=0.d0
  allocate(Ak_test(2,Lreal)); Ak_test=0.d0
  do ik=1,Lk
     call get_Akw_serial(ik,Sreal(:,1,1,1,1,:),Ak_test)
     Aloc = Aloc + Ak_test/dble(Lk)
     !     
     if(thop_empty.gt.0d0) then
        call get_Akw_empty(ik,Ak_test)
        Aloc = Aloc + Ak_test/dble(Lk)
     end if
     !
  end do
  if(master) then
     uio=free_unit()
     open(uio,file='test_get_akw.out')
     do iw=1,Lreal
        write(uio,'(5F18.10)') wr(iw),Aloc(:,iw)
     end do
     close(uio)
  end if
  !call mpi_stop
  !


  !+- set the eXciton frequency grid
  Lx=2*LwX+Lreal
  allocate(wrX(Lx)); wrX=0.d0    
  !+- coarse grid (negative frequency)
  allocate(wr_tmp(LwX))
  wr_tmp=linspace(winiX,wX-delta_wgrid,LwX,mesh=dw1,iend=.false.)
  wrX(1:Lwx) = wr_tmp 
  deallocate(wr_tmp)
  !+- fine grid (around exciton)
  allocate(wr_tmp(Lreal))
  wr_tmp=linspace(wX-delta_wgrid,wX+delta_wgrid,Lreal,mesh=dw2,iend=.false.)
  wrX(LwX+1:LwX+Lreal) = wr_tmp 
  deallocate(wr_tmp)
  !+- coarse grid (positive frequency)
  allocate(wr_tmp(LwX))
  wr_tmp=linspace(wX+delta_wgrid,wfinX,LwX,mesh=dw3)
  wrX(LwX+Lreal+1:2*LwX+Lreal) = wr_tmp 
  deallocate(wr_tmp)  
  !
  !+- the bare exciton propagator
  allocate(DX(LX)); DX=0.d0
  do iw=1,LX
     wcmplx=wrX(iw)+xi*epsX
     DX(iw) = 1d0/(wcmplx-wX)
  end do
  uio=free_unit()
  if(master) then
     open(uio,file='bare_exciton.out')
     do iw=1,LX
        write(uio,'(5F18.10)') wrX(iw),-1.d0*dimag(DX(iw)),DX(iw)
     end do
     close(uio)
  end if
  !
  !
  !
  !+- this has been tested
  ! allocate(wr_tmp(437))
  ! wr_tmp=linspace(wX-2*delta_wgrid,wX+2*delta_wgrid,437)
  ! do iw=1,437     
  !    call get_LR_freq(wr_tmp(iw),jL,jR)
  !    DX_tmp = DX(jL) + (DX(jR)-DX(jL))/(wrX(jR)-wrX(jL))*(wr_tmp(iw)-wrX(jL))
  !    if(master) write(800,*) wr_tmp(iw),wrX(jL),wrX(jR)     
  !    if(master) write(900,*) wr_tmp(iw),-1.d0*dimag(DX_tmp)
  ! end do
  ! deallocate(wr_tmp)  
  ! call mpi_stop
  

  call get_pi_irreducible(pi_irreducible)
  !
  !
  !+- here plot Pi_irreducible for q=0 -+!
  uio=free_unit()
  if(master) then
     open(uio,file='q0_pi_irreducible.out')
     do iw=1,LX
        write(uio,'(10F18.10)') wrX(iw),pi_irreducible(1:2,1,iw)
     end do
     close(uio)
  end if
  !  
  if(allocated(Lambda_vertex)) deallocate(Lambda_vertex)
  allocate(Lambda_vertex(2,Lk,LX)); Lambda_vertex=0d0
  if(tmp_debug_vertex) then
     Pi_irreducible = -1.d0*Pi_irreducible
     Lambda_vertex = -Vex*Vex*Pi_irreducible/(1.0-Vex*Pi_irreducible)
  else
     Lambda_vertex=Vex*Vex*Pi_irreducible/(1.0+Vex*Pi_irreducible)
  end if
  !
  !+- here compute the vertex
  !
  uio=free_unit()
  if(master) then
     if(tmp_debug_vertex) then
        do ik=1,Lk
           open(uio,file='lambda_vertex_ik'//str(ik,Npad=4)//'.out')
           do iw=1,LX
              write(uio,'(10F18.10)') wrX(iw),Lambda_vertex(1:2,ik,iw)
           end do
           close(uio)           
        end do
     else
        do ik=1,1
           open(uio,file='lambda_vertex_ik'//str(ik,Npad=4)//'.out')
           do iw=1,LX
              write(uio,'(10F18.10)') wrX(iw),Lambda_vertex(1:2,ik,iw)
           end do
           close(uio)           
        end do
     end if
  end if
  !
  !call mpi_stop
  !+- here compute the SigmaX
  call get_sigma_exciton(Sigma_exciton,Lambda_vertex)
  !

  uio=free_unit()
  if(master) then
     open(uio,file='SigmaX.out')
     do iw=1,LX
        write(uio,'(10F18.10)') wrX(iw),Sigma_exciton(1:2,iw)
     end do
     close(uio)
  end if

  !+- here do the dyson equation for the excitons
  allocate(DX_dressed(2,LX));DX_dressed=0d0
  DX_dressed(1,:) = 1d0/(DX**(-1d0)-sigma_exciton(1,:))
  DX_dressed(2,:) = 1d0/(DX**(-1d0)-sigma_exciton(2,:))

  uio=free_unit()
  if(master) then
     open(uio,file='dressedX.out')
     do iw=1,LX
        write(uio,'(10F18.10)') wrX(iw),DX_dressed(1:2,iw)
     end do
     close(uio)
  end if


  DX_dressed=0d0
  DX_dressed(1,:) = 1d0/(DX**(-1d0)-sigma_exciton(1,:)-Vex*dens(1))
  DX_dressed(2,:) = 1d0/(DX**(-1d0)-sigma_exciton(2,:)-Vex*dens(2))

  uio=free_unit()
  if(master) then
     open(uio,file='hartree_dressedX.out')
     do iw=1,LX
        write(uio,'(10F18.10)') wrX(iw),DX_dressed(1:2,iw)
     end do
     close(uio)
  end if

  !
  !
  ! uio=free_unit()
  ! if(master) then
  !    open(uio,file='q_pi_irreducible.out')
  !    do iw=1,LX
  !       write(uio,'(10F18.10)') wrX(iw),pi_irreducible(1:2,1,iw)
  !    end do
  !    close(uio)
  ! end if
  !
  call finalize_MPI()
  !
contains


  !+- first to test: get_Akw

  subroutine get_Akw(ik,Sigma,Akw)
    integer,intent(in) :: ik 
    complex(8),dimension(2,Lreal),intent(in) :: Sigma
    real(8),dimension(:,:),allocatable,intent(out) :: Akw
    complex(8),dimension(2,2,Lreal) :: Gkw,Gkw_tmp
    !
    !
    !+ 
    Gkw_tmp=0.d0
    do iw=1+rank,Lreal,mpiSize
       Gkw_tmp(:,:,iw) = (wr(iw)+xi*eps+xmu)*eye(2) - Hk_loop(:,:,ik)
       Gkw_tmp(1,1,iw) = Gkw_tmp(1,1,iw) - Sigma(1,iw)
       Gkw_tmp(2,2,iw) = Gkw_tmp(2,2,iw) - Sigma(2,iw)
       call inv(Gkw_tmp(:,:,iw))
    end do
    call mpi_allreduce(Gkw_tmp,Gkw,lreal*2*2,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)    

    if(allocated(Akw)) deallocate(Akw)
    allocate(Akw(2,Lreal)); Akw=0.d0
    !
    Akw(1,:) = -1.d0/pi*dimag(Gkw(1,1,:))
    Akw(2,:) = -1.d0/pi*dimag(Gkw(2,2,:))
    !
  end subroutine get_Akw



  subroutine get_Akw_serial(ik,Sigma,Akw)
    integer,intent(in) :: ik 
    complex(8),dimension(2,Lreal),intent(in) :: Sigma
    real(8),dimension(:,:),allocatable,intent(out) :: Akw
    complex(8),dimension(2,2,Lreal) :: Gkw
    !
    !
    !+ 
    Gkw=0.d0
    do iw=1,Lreal
       Gkw(:,:,iw) = (wr(iw)+xi*eps+xmu)*eye(2) - Hk_loop(:,:,ik)       
       Gkw(1,1,iw) = Gkw(1,1,iw) - Sigma(1,iw)
       Gkw(2,2,iw) = Gkw(2,2,iw) - Sigma(2,iw)
       call inv(Gkw(:,:,iw))
    end do
    if(allocated(Akw)) deallocate(Akw)
    allocate(Akw(2,Lreal)); Akw=0.d0
    !
    Akw(1,:) = -1.d0/pi*dimag(Gkw(1,1,:))
    Akw(2,:) = -1.d0/pi*dimag(Gkw(2,2,:))
    !
  end subroutine get_Akw_serial



  subroutine get_Akw_empty(ik,Akw)
    integer,intent(in) :: ik 
    real(8),dimension(:,:),allocatable,intent(out) :: Akw
    complex(8),dimension(2,2,Lreal) :: Gkw
    real(8) :: tmpWband
    !
    !
    !+ 
    if(allocated(Akw)) deallocate(Akw)
    allocate(Akw(2,Lreal)); Akw=0.d0
       
    if(.not.empty_quadratic) then
       Gkw=0.d0
       do iw=1,Lreal
          Gkw(:,:,iw) = (wr(iw)+xi*eps+xmu)*eye(2) 
          !+- add the ancillary empty band 
          Gkw(:,:,iw) = Gkw(:,:,iw) -  thop_empty*(Hk(:,:,ik) - pauli_z*Vel) - pauli_z*Vel - gap_empty*eye(2)                        
          call inv(Gkw(:,:,iw))
       end do
       !
       Akw(1,:) = -1.d0/pi*dimag(Gkw(1,1,:))
       Akw(2,:) = -1.d0/pi*dimag(Gkw(2,2,:))

    else
       ! just add a broad background
       tmpWband=10.d0*thop_empty


       do iw=1,Lreal
          Akw(1,iw) = fermi(-wr(iw)+gap_empty-tmpWband*0.5d0+Vel,beta_ancillary)*fermi(wr(iw)-gap_empty-tmpWband*0.5d0-Vel,beta_ancillary)*dos_ancillary/thop_empty
          Akw(2,iw) = fermi(-wr(iw)+gap_empty-tmpWband*0.5d0-Vel,beta_ancillary)*fermi(wr(iw)-gap_empty-tmpWband*0.5d0+Vel,beta_ancillary)*dos_ancillary/thop_empty
       end do
          !Gkw(:,:,iw) = Gkw(:,:,iw) -  0.5*h0*dot_product(kgrid(ik,:),kgrid(ik,:))/am_nm/am_nm/(thop_mev*1d-3) - pauli_z*Vel - gap_empty*eye(2) 
       
    end if

    !
  end subroutine get_Akw_empty




  subroutine get_Akw_w(ik,w_in,Sigma_in,Akw)
    integer,intent(in) :: ik 
    real(8),intent(in)    :: w_in
    complex(8),dimension(2),intent(in) :: Sigma_in
    real(8),dimension(2),intent(out) :: Akw
    complex(8),dimension(2,2) :: Gkw
    !
    !
    Gkw = (w_in+xi*eps+xmu)*eye(2) - Hk_loop(:,:,ik)
    Gkw(1,1) = Gkw(1,1) - Sigma_in(1)
    Gkw(2,2) = Gkw(2,2) - Sigma_in(2)
    call inv(Gkw)
    !
    Akw(1) = -1.d0/pi*dimag(Gkw(1,1))
    Akw(2) = -1.d0/pi*dimag(Gkw(2,2))
    !
  end subroutine get_Akw_w
  
  
  !+- this is the routine to be tested -+!
  subroutine get_pi_irreducible(pi_irreducible)
    complex(8),dimension(:,:,:),allocatable,intent(out) :: pi_irreducible
    complex(8),dimension(:,:,:),allocatable :: pi_irreducible_tmp
    real(8),dimension(:,:,:),allocatable :: re_pi_tmp,im_pi_tmp
    real(8) :: wtmp,wpX,wpX_,wtmp_rescale
    real(8),dimension(:,:),allocatable :: Akw,Akw_empty
    real(8),dimension(2) :: Akw_tmp
    integer :: ik,jk,iik,iw,jw
    complex(8),dimension(2) :: sigma_tmp
    complex(8),dimension(:,:),allocatable :: int_array
    real(8) :: win_tmp
    complex(8) :: wcmplx,wcmplx_
    !
    allocate(re_pi_tmp(2,Lk,Lx)); re_pi_tmp = 0.d0
    allocate(im_pi_tmp(2,Lk,Lx)); im_pi_tmp = 0.d0
    !
    if(allocated(pi_irreducible)) deallocate(pi_irreducible)
    allocate(pi_irreducible(2,Lk,Lx)); pi_irreducible = 0.d0
    allocate(pi_irreducible_tmp(2,Lk,Lx)); pi_irreducible_tmp = 0.d0
    !
    do ik=1+rank,Lk,mpiSize  !+- this is the external momentum
       !
       if(master) write(*,*) 'Pi-irreducible',ik,Lk       
       do jk=1,Lk   !+- this is the momentum over which the integration is performed
          !
          iik=ik_diff(ik,jk)   !+- this is the grid difference
          !
          call get_Akw_serial(iik,Sreal(:,1,1,1,1,:),Akw)  !+- here set the array of fermionic greens function
          if(thop_empty.gt.0d0) then
             call get_Akw_empty(iik,Akw_empty)  !+- here set the array of fermionic greens function
             Akw=Akw+Akw_empty
          end if
          !
          do iw=1,Lx          
             wcmplx = wrX(iw) + xi*epsX
             !
             wpX = wX + 0.5*h0*dot_product(kgrid(jk,:),kgrid(jk,:))/am_nm/am_nm !+- exciton energy
             !
             if(allocated(int_array)) deallocate(int_array)
             allocate(int_array(2,Lreal)); int_array=0.d0
             
             !+- do the integral over the dimensionless ferimionic frequencies
             !   at the end divide the result by thop
             wcmplx_ = wcmplx/thop_meV/1d-3
             wpX_ = wpX/thop_meV/1d-3
             !
             int_array(1,:) = (fermi(wr(:),beta)-1.d0)/(wcmplx_-wr(:)-wpX_)*Akw(1,:)
             int_array(2,:) = (fermi(wr(:),beta)-1.d0)/(wcmplx_-wr(:)-wpX_)*Akw(2,:)
             !
             !+- add to the integral
             pi_irreducible_tmp(1,ik,iw) = pi_irreducible_tmp(1,ik,iw) + & 
                  trapz(int_array(1,:),wr(1),wr(Lreal))/dble(Lk)/thop_meV/1d-3
             pi_irreducible_tmp(2,ik,iw) = pi_irreducible_tmp(2,ik,iw) + & 
                  trapz(int_array(2,:),wr(1),wr(Lreal))/dble(Lk)/thop_meV/1d-3
             !
          end do
       end do
       !
    end do
    CALL MPI_ALLREDUCE(pi_irreducible_tmp,pi_irreducible,2*Lk*Lx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)    
  end subroutine get_pi_irreducible



  subroutine get_sigma_exciton(SigmaX,L_vertex)
    complex(8),dimension(:,:),allocatable,intent(out) :: SigmaX
    complex(8),dimension(:,:,:),allocatable,intent(in) :: L_vertex    
    real(8),dimension(:,:),allocatable :: im_sigmaX_tmp,im_sigmaX
    real(8),dimension(2) :: ImL,ImL_l,ImL_r
    real(8),dimension(:,:), allocatable :: int_array
    real(8) :: wtmp
    integer :: ik,iw,jw,jwL,jwR
    real(8),dimension(:,:),allocatable :: Akw,Akw_empty
    
    if(.not.allocated(L_vertex)) call mpi_stop('if(.not.allocated(L_vertex))')
    if(size(L_vertex,1).ne.2) call mpi_stop('if(size(L_vertex,1).ne.2)')
    if(size(L_vertex,2).ne.Lk) call mpi_stop('if(size(L_vertex,2).ne.Lk)')
    if(size(L_vertex,3).ne.LX) call mpi_stop('if(size(L_vertex,3).ne.LX)')
    !
    allocate(im_sigmaX(2,Lx));     im_sigmaX=0.d0
    allocate(im_sigmaX_tmp(2,Lx)); im_sigmaX_tmp=0.d0    
    !allocate(int_array(2,Lx));int_array=0.d0
    !
    !
    do iw=1+rank,LX,mpiSize
       !
       im_sigmaX_tmp(:,iw) = 0.d0
       !
       if(master) write(*,*) 'sigmaX',iw,LX
       !+- integrate over momenta
       do ik=1,Lk
          !
          call get_Akw_serial(ik,Sreal(:,1,1,1,1,:),Akw)  !+- here set the array of fermionic greens function
          call get_Akw_empty(ik,Akw_empty)  !+- here set the array of fermionic greens function
          Akw = Akw + Akw_empty
          !
          if(allocated(int_array)) deallocate(int_array)
          allocate(int_array(2,Lreal)); int_array=0d0
          do jw=1,Lreal
             !
             wtmp=wrX(iw)+wr(jw)*thop_mev*1d-3
             if(wtmp.gt.winiX.and.wtmp.lt.wfinX) then
                ImL=0d0
                !
                !+- find the two closest left and right points on the grid and take the linear interpolation
                ! call get_LR_freq(wtmp,jwL,jwR)
                ! ImL_l = dimag(L_vertex(1:2,ik,jwL))
                ! ImL_r = dimag(L_vertex(1:2,ik,jwR))
                ! ImL(1:2) = ImL_l + (ImL_r-ImL_l)/(wrX(jwR)-wrX(jwL))*(wtmp-wrX(jwL))
                call linear_spline(wrX(:),dimag(L_vertex(1,ik,:)),wtmp,ImL(1))
                call linear_spline(wrX(:),dimag(L_vertex(2,ik,:)),wtmp,ImL(2))
                !
                int_array(1,jw) = (fermi(wr(jw),beta) - fermi(wtmp,beta/thop_mev/1d-3))*ImL(1)*Akw(1,jw)
                int_array(2,jw) = (fermi(wr(jw),beta) - fermi(wtmp,beta/thop_mev/1d-3))*ImL(2)*Akw(2,jw)

                ! !+- test
                ! if(abs(wrX(iw)-wX).lt.0.01) then
                !    write(500,*) wr(jw),(fermi(wr(jw),beta) - fermi(wtmp,beta/thop_mev/1d-3)),wtmp/thop_mev/1d-3
                ! end if
                !
             end if
             !
          end do
          !if(abs(wrX(iw)-wX).lt.0.01) call mpi_stop
          
          !
          ! if(tmp_debug_vertex) then
          !    im_sigmaX_tmp(1,iw) = im_sigmaX_tmp(1,iw) + trapz(int_array(1,:),wr(1),wr(Lreal))/dble(Lk)
          !    im_sigmaX_tmp(2,iw) = im_sigmaX_tmp(2,iw) + trapz(int_array(2,:),wr(1),wr(Lreal))/dble(Lk)
          ! else
          !    im_sigmaX_tmp(1,iw) = im_sigmaX_tmp(1,iw) - trapz(int_array(1,:),wr(1),wr(Lreal))/dble(Lk)
          !    im_sigmaX_tmp(2,iw) = im_sigmaX_tmp(2,iw) - trapz(int_array(2,:),wr(1),wr(Lreal))/dble(Lk)
          ! end if
          im_sigmaX_tmp(1,iw) = im_sigmaX_tmp(1,iw) - trapz(int_array(1,:),wr(1),wr(Lreal))/dble(Lk)
          im_sigmaX_tmp(2,iw) = im_sigmaX_tmp(2,iw) - trapz(int_array(2,:),wr(1),wr(Lreal))/dble(Lk)

          !
       end do
       !
    end do
    call mpi_allreduce(im_sigmaX_tmp,im_sigmaX,2*LX,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    im_sigmaX_tmp=0d0
    call get_kkt(im_sigmaX(1,:),im_sigmaX_tmp(1,:),wrX,'IR')    
    call get_kkt(im_sigmaX(2,:),im_sigmaX_tmp(2,:),wrX,'IR')
    if(allocated(SigmaX)) deallocate(SigmaX)
    allocate(SigmaX(2,LX));SigmaX=0d0    
    SigmaX=im_sigmaX_tmp+xi*im_sigmaX
    !
  end subroutine get_sigma_exciton


  subroutine get_LR_freq(wtmp,jw_tmp,jw_tmp_)
    real(8),intent(in) :: wtmp
    integer,intent(out) :: jw_tmp,jw_tmp_
    
    !+- find the two points on wrX grid  closest to wtmp -+!
    ! determine the sector of the grid 
    if(wtmp.le.wrX(LwX)) then !+- the negative sector
       !+-  the closest point to the left on the grid 
       jw_tmp  = 1 + floor((wtmp-winiX)/dw1)
       jw_tmp_ = 1 + ceiling((wtmp-winiX)/dw1)
       !
    else
       if(wtmp.le.wrX(LwX+Lreal)) then  !+- the exciton sector -+!
          !
          if(wtmp.lt.wrX(LwX+1)) then
             !+- between the grid 1 and 2
             jw_tmp  = LwX
             jw_tmp_ = LwX+1
          else
             !+- get the closest point on grid
             jw_tmp  = 1 + floor((wtmp-wrX(LwX+1))/dw2)                      
             jw_tmp_ = 1 + ceiling((wtmp-wrX(LwX+1))/dw2)                      
             !
             jw_tmp  = jw_tmp  + LwX
             jw_tmp_ = jw_tmp_ + LwX
             !
          end if
       else !+- the positive sector
          !
          if(wtmp.lt.wrX(LwX+Lreal+1)) then
             !+- between the grid 2 and 3
             jw_tmp  = LwX+Lreal
             jw_tmp_ = LwX+Lreal+1                         
          else                         
             jw_tmp  = 1 + floor((wtmp-wrX(LwX+Lreal+1))/dw3)
             jw_tmp_ = 1 + ceiling((wtmp-wrX(LwX+Lreal+1))/dw3)
             !
             jw_tmp  = jw_tmp  + LwX + Lreal
             jw_tmp_ = jw_tmp_ + LwX + Lreal
             !
          end if
          !
       end if
    end if
    

  end subroutine get_LR_freq

  

  subroutine get_KKT(ReS,ImS,wkkt,mode_)
    real(8),dimension(:) :: ReS,ImS
    real(8),dimension(:) :: wkkt
    character(len=2),optional :: mode_
    character(len=2) :: mode
    real(8) :: A,B
    integer :: iv,iw,Lw
    real(8),dimension(:),allocatable :: ImS_tmp
    !
    Lw = size(wkkt)
    if(size(ReS).ne.Lw) then
       if(rank==0) write(*,*) 'size(ReS).ne.Lw'
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)       
       stop
    end if
    !
    if(size(ImS).ne.Lw) then
       if(rank==0) write(*,*) 'size(ImS).ne.Lw'
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)       
       stop
    end if
    !
    mode='RI'
    if(present(mode_)) mode=mode_
    if(mode.ne.'IR'.and.mode.ne.'RI') stop "wrong mode KKT"
    !
    allocate(ImS_tmp(Lw))
    ImS_tmp=0.d0
    ImS=0.d0
    do iw=1+rank,Lw,mpiSize
       do iv=1,Lw-1          
          A = ReS(iv) -wkkt(iv)*(ReS(iv)-ReS(iv+1))/(wkkt(iv)-wkkt(iv+1))
          B = (ReS(iv)-ReS(iv+1))/(wkkt(iv)-wkkt(iv+1))          
          ImS_tmp(iw) =ImS_tmp(iw)  -B*(wkkt(iv+1)-wkkt(iv))          
          if(iv+1.ne.iw) ImS_tmp(iw) = ImS_tmp(iw) - (A+B*wkkt(iw))*log(abs(wkkt(iw)-wkkt(iv+1)))
          if(iv.ne.iw)   ImS_tmp(iw) = ImS_tmp(iw) + (A+B*wkkt(iw))*log(abs(wkkt(iw)-wkkt(iv)))
       end do
       ImS_tmp(iw) = ImS_tmp(iw)/pi
       if(mode.eq.'IR') ImS_tmp(iw)=-ImS_tmp(iw)
    end do
    CALL MPI_ALLREDUCE(ImS_tmp,ImS,Lw,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)
    !
  end subroutine get_KKT




  subroutine get_KKT_serial(ReS,ImS,wkkt,mode_)
    real(8),dimension(:) :: ReS,ImS
    real(8),dimension(:) :: wkkt
    character(len=2),optional :: mode_
    character(len=2) :: mode
    real(8) :: A,B
    integer :: iv,iw,Lw
    real(8),dimension(:),allocatable :: ImS_tmp
    !
    Lw = size(wkkt)
    if(size(ReS).ne.Lw) then
       if(rank==0) write(*,*) 'size(ReS).ne.Lw'
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)       
       stop
    end if
    !
    if(size(ImS).ne.Lw) then
       if(rank==0) write(*,*) 'size(ImS).ne.Lw'
       CALL MPI_BARRIER(MPI_COMM_WORLD,MPIerr)       
       stop
    end if
    !
    mode='RI'
    if(present(mode_)) mode=mode_
    if(mode.ne.'IR'.and.mode.ne.'RI') stop "wrong mode KKT"
    !
    ImS=0.d0
    !do iw=1+rank,Lw,mpiSize
    do iw=1,Lw
       do iv=1,Lw-1          
          A = ReS(iv) -wkkt(iv)*(ReS(iv)-ReS(iv+1))/(wkkt(iv)-wkkt(iv+1))
          B = (ReS(iv)-ReS(iv+1))/(wkkt(iv)-wkkt(iv+1))          
          ImS(iw) =ImS(iw)  -B*(wkkt(iv+1)-wkkt(iv))          
          if(iv+1.ne.iw) ImS(iw) = ImS(iw) - (A+B*wkkt(iw))*log(abs(wkkt(iw)-wkkt(iv+1)))
          if(iv.ne.iw)   ImS(iw) = ImS(iw) + (A+B*wkkt(iw))*log(abs(wkkt(iw)-wkkt(iv)))
       end do
       ImS(iw) = ImS(iw)/pi
       if(mode.eq.'IR') ImS(iw)=-ImS(iw)
    end do
    !
  end subroutine get_KKT_serial



  subroutine mpi_stop(msg)
    character(len=*),optional :: msg
    if(master) then
       write(*,*) "STOP!"
       if(present(msg)) write(*,*) msg       
    end if
    call mpi_barrier(comm,mpiERR)
    stop
  end subroutine mpi_stop
  
  function get_delta_dens_imp(xmu_in) result(dens_out)
    real(8),intent(in) :: xmu_in
    real(8) :: dens_out
    character(len=20) :: suffix
    !
    xmu=xmu_in
    !
    call ed_solve(comm,bath,Hloc,mpi_lanc=flag_mpi,iflag=.false.)
    call ed_get_dens(dens,Nlat,iorb=1)
    !
    call ed_reset_suffix    
    dens_out=dens(1)+dens(2)-ntarget        
    !
  end function get_delta_dens_imp


  function get_delta_dens_imp_(xmu_in) result(dens_out)
    real(8) :: xmu_in
    real(8) :: dens_out
    character(len=20) :: suffix
    !
    xmu=xmu_in
    !
    call ed_solve(comm,bath,Hloc,iflag=.false.)
    call ed_get_dens(dens,Nlat,iorb=1)
    !
    call ed_reset_suffix    
    dens_out=(dens(1)+dens(2)-ntarget)**2.d0        
    !
  end function get_delta_dens_imp_
    
  

  subroutine build_hk_honeycomb(file)
    character(len=*),optional          :: file
    !real(8),dimension(2)               :: pointK,pointKp,pointM,kdiff
    real(8),dimension(:,:),allocatable :: kgrid_tmp
    real(8),dimension(:),allocatable   :: gridx,gridy
    integer                            :: i,j,ik,unit_io,jk
    integer :: idiff,jdiff
    real(8) :: itmp,jtmp
    
    !
    Lk= Nk*Nk
    if(master)write(*,*)"Build H(k) for the honeycomb lattice:",Lk
    if(master)write(*,*)"# of SO-bands     :",Nlso
    !
    if(allocated(Hk))deallocate(Hk)
    !
    allocate(Hk(Nlso,Nlso,Lk));Hk=zero
    !
    call TB_build_model(Hk,hk_honeycomb,Nlso,[Nk,Nk],iprint=.true.,kgrid_out=kgrid_tmp)
    if(present(file).and.master) call TB_write_Hk(Hk,file,1,Nspin,Norb,[Nk,Nk])
    !
    pointK  = 1d0/3d0*bk1 + 2d0/3d0*bk2
    pointKp = 2d0/3d0*bk1 + 1d0/3d0*bk2
    pointM= 0.5d0*pointK+0.5d0*pointKp
    write(*,*) pointK
    write(*,*) pointKp
    allocate(Kpath(4,2))
    KPath(1,:)=[0,0]
    KPath(2,:)=pointK
    Kpath(3,:)=pointM
    KPath(4,:)=[0d0,0d0]
    !
    call TB_Solve_model(hk_honeycomb,Nlso,KPath,Nkpath,&
         colors_name=[red1,blue1,red1,blue1],&
         points_name=[character(len=10) :: "G","K","K`","G"],&
         file="Eigenbands.nint",iproject=.false.)
    !
    
    allocate(gridx(Nk)); gridx=linspace(0d0,1d0,Nk,iend=.false.)
    allocate(gridy(Nk)); gridy=linspace(0d0,1d0,Nk,iend=.false.)

    allocate(kgrid(Lk,2)); kgrid=0.d0
    allocate(ik2ij(Lk,2)); ik2ij=0
    allocate(ij2ik(Nk,Nk)); ij2ik=0

    ik=0
    do j=1,Nk
       do i=1,Nk
          ik=ik+1
          kgrid(ik,:) = gridx(i)*bk1+gridx(j)*bk2
          ik2ij(ik,1) = i
          ik2ij(ik,1) = j
          ij2ik(i,j)  = ik
       end do       
    end do
    !
    if(master) then
       unit_io=free_unit()
       open(unit_io,file='kgrid_driver.out')
       do ik=1,Lk
          write(unit_io,'(5F18.10)') kgrid_tmp(ik,:),kgrid(ik,:)
       end do
       close(unit_io)
    end if
    !
    !
    ! if(test_strides.gt.Lk) test_strides=0
    ! if(test_strides.lt.1) test_strides=0

    allocate(ik_diff(Lk,Lk)) ; ik_diff=0
    do ik=1,Lk
       do jk=1,Lk
          ! 
          kdiff = kgrid(ik,:) - kgrid(jk,:)
          !
          itmp = dot_product(e1,kdiff)/2.d0/pi
          jtmp = dot_product(e2,kdiff)/2.d0/pi
          !
          idiff = nint(Nk*itmp)+1
          jdiff = nint(Nk*jtmp)+1          

          !
          do while(idiff.gt.Nk) 
             idiff = idiff - Nk
          end do
          do while(idiff.lt.1) 
             idiff = idiff + Nk
          end do
          !
          do while(jdiff.gt.Nk) 
             jdiff = jdiff - Nk
          end do
          do while(jdiff.lt.1) 
             jdiff = jdiff + Nk
          end do
          !
          ! write(700,*) jtmp
          ! do while(jtmp.ge.1.d0) 
          !    jtmp = jtmp - 1.d0
          ! end do
          ! do while(jtmp.lt.0.d0) 
          !    jtmp = jtmp + 1.d0
          ! end do
          ! write(701,*) jtmp


          if(idiff.lt.1.or.idiff.gt.Nk) then
             if(master) write(*,*) idiff,itmp,ik,jk
             call mpi_stop('idiff.lt.1.or.idiff.gt.Nk')
          end if
          if(jdiff.lt.1.or.jdiff.gt.Nk) then
             if(master) write(*,*) jdiff,jtmp,ik,jk
             call mpi_stop('jdiff.lt.1.or.jdiff.gt.Nk')
          end if
          ik_diff(ik,jk) = ij2ik(idiff,jdiff)
          if(ik.eq.test_strides) then             
             write(800,*) kgrid(ik,:),kgrid(jk,:),kdiff,kgrid(ik_diff(ik,jk),:)
             !
             write(900,*) kgrid(ik,:)
             write(900,*) kgrid(jk,:)
             write(900,*) kgrid(ik_diff(ik,jk),:)
             write(900,*)
             write(900,*)
             !
             write(901,*) kdiff
             write(901,*)
             write(901,*)
             !
          end if
       end do
    end do
    !
    !
    !
  end subroutine build_hk_honeycomb

  subroutine set_hloc(nave,file_out)
    real(8),dimension(2) :: nave
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hloc_tmp
    character(len=*),optional :: file_out
    integer :: unit,iso,ik
    !
    Hloc_tmp = zero
    !+- sum the k Hk -+!
    Hloc_tmp = sum(Hk,dim=3)/Lk    
    where(abs(Hloc_tmp)<1d-6) Hloc_tmp=zero
    !
    if(allocated(Hloc)) deallocate(Hloc)    
    allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))    
    !
    if(allocated(Hk_loop)) deallocate(Hk_loop)
    allocate(Hk_loop(Nlso,Nlso,Lk)) ; Hk_loop=Hk
    !
    iso=lso2j_index(1,1,1)
    Hloc_tmp(iso,iso) = Hloc_tmp(iso,iso) + 3.d0*Vnn*nave(2) !
    !
    iso=lso2j_index(2,1,1)
    Hloc_tmp(iso,iso) = Hloc_tmp(iso,iso) + 3.d0*Vnn*nave(1) !
    !
    do ik=1,Lk
       iso=lso2j_index(1,1,1)
       Hk_loop(iso,iso,ik) = Hk_loop(iso,iso,ik) + 3.d0*Vnn*nave(2) !
       !
       iso=lso2j_index(2,1,1)
       Hk_loop(iso,iso,ik) = Hk_loop(iso,iso,ik) + 3.d0*Vnn*nave(1) !
    end do
    !
    Hloc = lso2nnn_reshape(Hloc_tmp,Nlat,Nspin,Norb)
    !
    if(present(file_out).and.master) then
       unit=free_unit()
       open(unit,file=trim(file_out))
       do ilat=1,Nlat
          write(unit,'(10F18.10)') Hloc(ilat,:,:,:,:)
       end do
       close(unit)
    end if
    !
  end subroutine set_hloc

  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_SigmaHk(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    SigmaHk = zero;if(present(sigma))SigmaHk=sigma
  end subroutine set_SigmaHk
  !
  function hk_honeycomb(kpoint,Nlso) result(hk)
    real(8),dimension(:)            :: kpoint
    integer                         :: Nlso
    complex(8),dimension(2,2)       :: hk11,hk22
    complex(8),dimension(Nlso,Nlso) :: hk
    real(8)                         :: h0,hx,hy,hz, kdote1, kdote2
    !
    kdote1 = dot_product(kpoint,e1)
    kdote2 = dot_product(kpoint,e2)
    !
    h0 = -2*ts*( cos(kdote1) + cos(kdote2) + cos(kdote1+kdote2) )
    hx = -tperp*( cos(kdote1) + cos(kdote2) + 1)
    hy = -tperp*( sin(kdote1) + sin(kdote2) )
    hz = Vel
    !
    hk = h0*pauli_0 + hx*pauli_x + hy*pauli_y + hz*pauli_z 
    !
  end function hk_honeycomb
  



  function lso2nnn_reshape(Fin,Nlat,Nspin,Norb) result(Fout)
    integer                                               :: Nlat,Nspin,Norb
    complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Fin
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fout
    integer                                               :: iorb,ispin,ilat,is
    integer                                               :: jorb,jspin,js
    Fout=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Fout(ilat,ispin,jspin,iorb,jorb) = Fin(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape



  function lso2j_index(ilat,ispin,iorb) result(isporb)
    integer :: ispin,iorb,ilat
    integer :: isporb
    if(iorb>Norb)stop "error so2j_index: iorb>Norb"
    if(ispin>Nspin)stop "error so2j_index: ispin>Nspin"
    isporb=(ilat-1)*Norb*Nspin + (ispin-1)*Nspin + iorb
  end function lso2j_index


  function lso2j(fg,Nlat) result(g)
    integer :: Nlat
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: fg
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: g
    integer                                     :: i,j,iorb,jorb,ispin,jspin,ilat
    g=zero
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   i=lso2j_index(ilat,ispin,iorb)
                   j=lso2j_index(ilat,jspin,jorb)
                   g(ilat,i,j) = fg(ilat,ispin,jspin,iorb,jorb)
                enddo
             enddo
          enddo
       enddo
    end do
  end function lso2j

  function j2lso(fg,Nlat) result(g)
    integer :: Nlat
    complex(8),dimension(Nlat,Nspin*Norb,Nspin*Norb) :: fg
    complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb) :: g
    integer                                          :: i,j,iorb,jorb,ispin,jspin,ilat
    do ilat=1,Nlat
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   i=lso2j_index(ilat,ispin,iorb)
                   j=lso2j_index(ilat,jspin,jorb)
                   g(ilat,ispin,jspin,iorb,jorb)  = fg(ilat,i,j)
                enddo
             enddo
          enddo
       enddo
    end do
  end function j2lso


end program
        ! uio=free_unit()
        ! if(master) open(uio,file='bracket_xmu'//str(iloop,Npad=3)//'.dat')
        ! if(iloop>1) then
        !    xmu1=xmu_imp
        !    xmu2=xmu1
        ! end if
        ! !+- check where to move - +!
        ! ntest1=get_delta_dens_imp(xmu1)       
        ! ntest2=ntest1
        ! alpha_mu= abs(alpha_mu)  !+- assuming starting from a system short of particles
        ! if(ntest1.gt.0) alpha_mu=-1.d0*abs(alpha_mu)  !+- if not adjust
        ! imu=0
        ! do while(ntest2*ntest1.ge.0.d0.and.imu.le.100)
        !    xmu1=xmu2
        !    ntest1=ntest2
        !    !
        !    xmu2=xmu2+alpha_mu
        !    ntest2=get_delta_dens_imp(xmu2)  
        !    if(master) write(uio,'(5F18.10)') xmu2,ntest2,ntest1
        !    imu=imu+1
        ! end do
        ! if(master) close(uio)
        ! xmu_imp = brentq(get_delta_dens_imp,xmu1,xmu2)            
        ! xmu=xmu_imp
        ! xmu_imp=50
        ! do imu=1,100
        !    ntest1=get_delta_dens_imp_(xmu_imp)
        !    if(master) write(450+iloop,*) xmu_imp,ntest1
        !    xmu_imp=xmu_imp-1.0
        ! end do
        ! if(iloop.eq.2) then
        !    call mpi_barrier(comm,mpiERR)
        !    stop
        ! end if
        ! xmu=1
  ! subroutine get_pi_irreducible(pi_irreducible)
  !   complex(8),dimension(:,:,:),allocatable,intent(out) :: pi_irreducible
  !   complex(8),dimension(:,:,:),allocatable :: pi_irreducible_tmp
  !   real(8),dimension(:,:,:),allocatable :: re_pi_tmp,im_pi_tmp
  !   real(8) :: wtmp,wpX,wtmp_rescale
  !   real(8),dimension(:,:),allocatable :: Akw
  !   real(8),dimension(2) :: Akw_tmp
  !   integer :: ik,jk,iik,iw,jw
  !   complex(8),dimension(2) :: sigma_tmp
  !   complex(8),dimension(:,:),allocatable :: int_array
  !   real(8) :: win_tmp
  !   complex(8) :: wcmplx
  !   !
  !   allocate(re_pi_tmp(2,Lk,Lx)); re_pi_tmp = 0.d0
  !   allocate(im_pi_tmp(2,Lk,Lx)); im_pi_tmp = 0.d0
  !   !
  !   if(allocated(pi_irreducible)) deallocate(pi_irreducible)
  !   allocate(pi_irreducible(2,Lk,Lx)); pi_irreducible = 0.d0
  !   allocate(pi_irreducible_tmp(2,Lk,Lx)); pi_irreducible_tmp = 0.d0
  !   !
  !   do ik=1+rank,1,mpiSize  !+- this is the external momentum
  !      !
  !      if(master) write(*,*) 'Pi-irreducible',ik,Lk
       
  !      do jk=1,Lk   !+- this is the momentum over which the integration is performed
  !         !

  !         if(master) write(*,*) 'jk',jk


  !         iik=ik_diff(ik,jk)   !+- this is the grid difference
  !         !
  !         call get_Akw_serial(iik,Sreal(:,1,1,1,1,:),Akw)  !+- here set the array of fermionic greens function
  !         !
  !         do iw=1,Lx          
  !            wcmplx = wrX(iw) + xi*epsX
  !            !
  !            wpX = wX + 0.5*h0*dot_product(kgrid(jk,:),kgrid(jk,:))/am_nm/am_nm !+- exciton energy
  !            !
  !            if(allocated(int_array)) deallocate(int_array)
  !            allocate(int_array(2,Lx)); int_array=0.d0
             
  !            do jw=1,Lx
  !               !
  !               wtmp = wrX(jw)-wpX
  !               !get  wtmp in dimensionless units
  !               !wtmp is in electron Volt: rescale using the thop=1 dimensionless scale
  !               wtmp_rescale = wtmp/thop_meV/1d-3

  !               !
  !               Akw_tmp=0d0
  !               if(wtmp_rescale.gt.wini.and.wtmp_rescale.lt.wfin) then
  !                  !
  !                  !+- here interpolate the spectral functions
  !                  call linear_spline(wr(:),Akw(1,:),wtmp_rescale,Akw_tmp(1))   
  !                  call linear_spline(wr(:),Akw(2,:),wtmp_rescale,Akw_tmp(2))
  !                  !
  !               ! else
  !               !    !+- approx the self-energy with [wini] and [wfin] values
  !               !    if(wtmp_rescale.lt.wini) then                
  !               !       Sigma_tmp=Sreal(:,1,1,1,1,1)
  !               !       !win_tmp = wini
  !               !    else
  !               !       Sigma_tmp=Sreal(:,1,1,1,1,Lreal)
  !               !       !win_tmp = wfin
  !               !    end if
  !               !    call get_Akw_w(iik,wtmp_rescale,Sigma_tmp,Akw_tmp)   
  !               end if
  !               !+- Akw_tmp is in dimensionless units: rescale back to eV^{-1}
  !               !   Akw_tmp[ev^{-1}] = Ake_tmp/thop(meV)/1d-3
  !               Akw_tmp = Akw_tmp/thop_mev/1d-3 
  !               !
  !               int_array(:,jw) = (fermi(wtmp_rescale,beta)-1.d0)/(wcmplx-wrX(jw))*Akw_tmp(:)
  !            end do

  !            !+- add to the integral
  !            pi_irreducible_tmp(1,ik,iw) = pi_irreducible_tmp(1,ik,iw) + & 
  !                 trapz(int_array(1,:),wrX(1),wrX(Lx))/dble(Lk)
  !            pi_irreducible_tmp(2,ik,iw) = pi_irreducible_tmp(2,ik,iw) + & 
  !                 trapz(int_array(2,:),wrX(1),wrX(Lx))/dble(Lk)

  !            !
  !         end do
  !      end do
  !      !
  !      !kkt
  !      ! call get_kkt_serial(im_pi_tmp(1,ik,:),re_pi_tmp(1,ik,:),wrX,'IR')
  !      ! call get_kkt_serial(im_pi_tmp(2,ik,:),re_pi_tmp(2,ik,:),wrX,'IR')       
  !      ! !
  !      ! pi_irreducible_tmp(:,ik,:) = xi*im_pi_tmp(:,ik,:)+re_pi_tmp(:,ik,:)
  !      !
  !   end do
  !   CALL MPI_ALLREDUCE(pi_irreducible_tmp,pi_irreducible,2*Lk*Lx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)    
  ! end subroutine get_pi_irreducible





  ! subroutine get_pi_irreducible_kkt(pi_irreducible)
  !   complex(8),dimension(:,:,:),allocatable,intent(out) :: pi_irreducible
  !   complex(8),dimension(:,:,:),allocatable :: pi_irreducible_tmp
  !   real(8),dimension(:,:,:),allocatable :: re_pi_tmp,im_pi_tmp
  !   real(8) :: wtmp,wpX,wtmp_rescale
  !   real(8),dimension(:,:),allocatable :: Akw
  !   real(8),dimension(2) :: Akw_tmp
  !   integer :: ik,jk,iik
  !   complex(8),dimension(2) :: sigma_tmp
  !   real(8) :: win_tmp
  !   !
  !   allocate(re_pi_tmp(2,Lk,Lx)); re_pi_tmp = 0.d0
  !   allocate(im_pi_tmp(2,Lk,Lx)); im_pi_tmp = 0.d0
  !   !
  !   if(allocated(pi_irreducible)) deallocate(pi_irreducible)
  !   allocate(pi_irreducible(2,Lk,Lx)); pi_irreducible = 0.d0
  !   allocate(pi_irreducible_tmp(2,Lk,Lx)); pi_irreducible_tmp = 0.d0
  !   !
  !   do ik=1+rank,1,mpiSize  !+- this is the external momentum
  !      !
  !      if(master) write(*,*) 'Pi-irreducible',ik,Lk
       
  !      do jk=1,Lk   !+- this is the momentum over which the integration is performed
  !         !
  !         iik=ik_diff(ik,jk)   !+- this is the grid difference
  !         !
  !         call get_Akw_serial(iik,Sreal(:,1,1,1,1,:),Akw)  !+- here set the array of fermionic greens function
  !         !
  !         do iw=1,Lx          
  !            !
  !            wpX = wX + 0.5*h0*dot_product(kgrid(jk,:),kgrid(jk,:))/am_nm/am_nm !+- exciton energy
  !            wtmp = wrX(iw) - wpX
             
  !            !get  wtmp in dimensionless units
  !            !wtmp is in electron Volt: rescale using the thop=1 dimensionless scale
  !            wtmp_rescale = wtmp/thop_meV/1d-3
  !            !
  !            Akw_tmp=0d0
  !            if(wtmp_rescale.gt.wini.and.wtmp_rescale.lt.wfin) then
  !               !
  !               !+- here interpolate the spectral functions
  !               call linear_spline(wr(:),Akw(1,:),wtmp_rescale,Akw_tmp(1))   
  !               call linear_spline(wr(:),Akw(2,:),wtmp_rescale,Akw_tmp(2))
  !               !
  !            else
  !               !+- approx the self-energy with [wini] and [wfin] values
  !               if(wtmp_rescale.lt.wini) then                
  !                  Sigma_tmp=Sreal(:,1,1,1,1,1)
  !                  !win_tmp = wini
  !               else
  !                  Sigma_tmp=Sreal(:,1,1,1,1,Lreal)
  !                  !win_tmp = wfin
  !               end if
  !               call get_Akw_w(iik,wtmp_rescale,Sigma_tmp,Akw_tmp)   
  !            end if
  !            !+- Akw_tmp is in dimensionless units: rescale back to eV^{-1}
  !            !   Akw_tmp[ev^{-1}] = Ake_tmp/thop(meV)/1d-3
  !            Akw_tmp = Akw_tmp/thop_mev/1d-3 
  !            !
  !            !+- add to the integral             
  !            im_pi_tmp(:,ik,iw) = im_pi_tmp(:,ik,iw) - pi*Akw_tmp(:)*(fermi(wtmp,beta/thop_mev/1d-3)-1.d0)/dble(Lk) ! the -1.0 is  bose(-\o_p)             
  !         end do
  !      end do
  !      !
  !      !kkt
  !      call get_kkt_serial(im_pi_tmp(1,ik,:),re_pi_tmp(1,ik,:),wrX,'IR')
  !      call get_kkt_serial(im_pi_tmp(2,ik,:),re_pi_tmp(2,ik,:),wrX,'IR')       
  !      !
  !      pi_irreducible_tmp(:,ik,:) = xi*im_pi_tmp(:,ik,:)+re_pi_tmp(:,ik,:)
  !      !
  !   end do
  !   CALL MPI_ALLREDUCE(pi_irreducible_tmp,pi_irreducible,2*Lk*Lx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)    
  ! end subroutine get_pi_irreducible_kkt




  ! subroutine get_pi_irreducible_kkt_(pi_irreducible)
  !   complex(8),dimension(:,:,:),allocatable,intent(out) :: pi_irreducible
  !   complex(8),dimension(:,:,:),allocatable :: pi_irreducible_tmp
  !   real(8),dimension(:,:,:),allocatable :: re_pi_tmp,im_pi_tmp
  !   real(8) :: wtmp,wpX,wtmp_rescale
  !   real(8),dimension(:,:),allocatable :: Akw
  !   real(8),dimension(2) :: Akw_tmp
  !   integer :: ik,jk,iik,iw,jw
  !   complex(8),dimension(2) :: sigma_tmp
  !   real(8) :: win_tmp

  !   real(8),dimension(:),allocatable :: int_array
  !   !complex(8) :: wcmplx
  !   !
  !   real(8),allocatable,dimension(:,:) :: fermi_bose,fermi_bose_tmp
    
  !   !
  !   allocate(fermi_bose(LK,Lx));     fermi_bose=0d0
  !   allocate(fermi_bose_tmp(LK,Lx)); fermi_bose_tmp=0d0
  !   allocate(int_array(LX)) ; int_array=0d0
  !   !
  !   do ik=1+rank,Lk,mpiSize  !+- this is the external momentum
  !      !
  !      wpX = wX + 0.5*h0*dot_product(kgrid(ik,:),kgrid(ik,:))/am_nm/am_nm !+- exciton energy
  !      !
  !      do iw=1,Lx
  !         !
  !         int_array=0d0
  !         do jw=1,Lx
  !            int_array(jw) = -epsX*(fermi(wrX(jw),beta/thop_mev/1d-3)-1d0)/((wrX(iw)-wrX(jw)-wpX)**2.d0+epsX**2.d0)
  !         end do          
  !         fermi_bose_tmp(ik,iw) = trapz(int_array,wrX(1),wrX(Lx))          
  !      end do
  !   end do
  !   CALL MPI_ALLREDUCE(fermi_bose_tmp,fermi_bose,Lk*Lx,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD,MPIerr)    
    
  !   !
  !   allocate(re_pi_tmp(2,Lk,Lx)); re_pi_tmp = 0.d0
  !   allocate(im_pi_tmp(2,Lk,Lx)); im_pi_tmp = 0.d0
  !   !
  !   if(allocated(pi_irreducible)) deallocate(pi_irreducible)
  !   allocate(pi_irreducible(2,Lk,Lx)); pi_irreducible = 0.d0
  !   allocate(pi_irreducible_tmp(2,Lk,Lx)); pi_irreducible_tmp = 0.d0
  !   !
  !   do ik=1+rank,1,mpiSize  !+- this is the external momentum
  !      !
  !      if(master) write(*,*) 'Pi-irreducible',ik,Lk
       
  !      do jk=1,Lk   !+- this is the momentum over which the integration is performed
  !         !
  !         iik=ik_diff(ik,jk)   !+- this is the grid difference
  !         !
  !         call get_Akw_serial(iik,Sreal(:,1,1,1,1,:),Akw)  !+- here set the array of fermionic greens function
  !         !
  !         do iw=1,Lx          
  !            !
  !            wpX = wX + 0.5*h0*dot_product(kgrid(jk,:),kgrid(jk,:))/am_nm/am_nm !+- exciton energy
  !            wtmp = wrX(iw) - wpX
             
  !            !get  wtmp in dimensionless units
  !            !wtmp is in electron Volt: rescale using the thop=1 dimensionless scale
  !            wtmp_rescale = wtmp/thop_meV/1d-3
  !            !
  !            Akw_tmp=0d0
  !            if(wtmp_rescale.gt.wini.and.wtmp_rescale.lt.wfin) then
  !               !
  !               !+- here interpolate the spectral functions
  !               call linear_spline(wr(:),Akw(1,:),wtmp_rescale,Akw_tmp(1))   
  !               call linear_spline(wr(:),Akw(2,:),wtmp_rescale,Akw_tmp(2))
  !            !    !
  !            else
  !               !+- approx the self-energy with [wini] and [wfin] values
  !               if(wtmp_rescale.lt.wini) then                
  !                  Sigma_tmp=Sreal(:,1,1,1,1,1)
  !                  !win_tmp = wini
  !               else
  !                  Sigma_tmp=Sreal(:,1,1,1,1,Lreal)
  !                  !win_tmp = wfin
  !               end if
  !               call get_Akw_w(iik,wtmp_rescale,Sigma_tmp,Akw_tmp)   
  !            end if
  !            !+- Akw_tmp is in dimensionless units: rescale back to eV^{-1}
  !            !   Akw_tmp[ev^{-1}] = Ake_tmp/thop(meV)/1d-3
  !            Akw_tmp = Akw_tmp/thop_mev/1d-3 
  !            !
  !            !+- add to the integral: the fermi-bose factor has been computed once for all above
  !            !im_pi_tmp(:,ik,iw) = im_pi_tmp(:,ik,iw) - pi*Akw_tmp(:)*(fermi(wtmp,beta)-1.d0)/dble(Lk) ! the -1.0 is  bose(-\o_p)             
  !            im_pi_tmp(:,ik,iw) = im_pi_tmp(:,ik,iw) + Akw_tmp(:)*fermi_bose(jk,iw)/dble(Lk)       
  !         end do
  !      end do
  !      !
  !      !kkt
  !      call get_kkt_serial(im_pi_tmp(1,ik,:),re_pi_tmp(1,ik,:),wrX,'IR')
  !      call get_kkt_serial(im_pi_tmp(2,ik,:),re_pi_tmp(2,ik,:),wrX,'IR')       
  !      !
  !      pi_irreducible_tmp(:,ik,:) = xi*im_pi_tmp(:,ik,:)+re_pi_tmp(:,ik,:)
  !      !
  !   end do
  !   CALL MPI_ALLREDUCE(pi_irreducible_tmp,pi_irreducible,2*Lk*Lx,MPI_DOUBLE_COMPLEX,MPI_SUM,MPI_COMM_WORLD,MPIerr)    
  ! end subroutine get_pi_irreducible_kkt_
