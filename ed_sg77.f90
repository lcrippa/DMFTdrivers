program ed_sg77
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI

  implicit none
  integer                                     :: iloop,Lk,Nso,Lakw
  logical                                     :: converged
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,print_mode
  real(8),allocatable                         :: Bath(:),Bath_(:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:),Weiss_(:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:),Sreal(:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:),Greal(:,:,:,:,:)
  complex(8),allocatable,dimension(:)         :: Gtest
  !hamiltonian input:
  complex(8),allocatable                      :: Hk(:,:,:),sg77Hloc(:,:),sigmaBHZ(:,:),Zmats(:,:),impRho(:,:)
  real(8),allocatable                         :: Wtk(:)
  real(8),allocatable                         :: kxgrid(:),kygrid(:)
  integer,allocatable                         :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: mh,lambda,wmixing,akrange
  character(len=30)                           :: Params
  character(len=16)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,usez,mixG0
  !
  real(8),dimension(2)                        :: Eout
  real(8),allocatable                         :: dens(:)
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma5,GammaN
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  complex(8),dimension(4,4)                   :: GammaR0,GammaRx,GammaRy,GammaRz
  real(8),dimension(:),allocatable            :: lambdasym_vector
  complex(8),dimension(:,:,:,:,:),allocatable :: Hsym_basis
  !MPI Vars:
  integer                                     :: irank,comm,rank,size2,ierr
  logical                                     :: master,getbands,getakw

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)

  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_SG77.in')  
  call parse_input_variable(Params,"Params",finput,default="E0EzEx",&
       comment="Ex; EzEx; E0Ex; ExEy; E0Ez; E0EzEx; E0EzExEy")
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(mh,"MH",finput,default=0.d0)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.d0)
  call parse_input_variable(usez,"USEZ",finput,default=.false.)
  call parse_input_variable(getbands,"GETBANDS",finput,default=.false.)
  call parse_input_variable(getakw,"GETAKW",finput,default=.false.)
  call parse_input_variable(Lakw,"LAKW",finput,default=250)
  call parse_input_variable(akrange,"AKRANGE",finput,default=5d0)
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

  if(Nspin/=2.OR.Norb/=2)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats));Weiss=zero
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats));Weiss_=zero
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats));Smats=zero
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats));Gmats=zero
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal));Sreal=zero
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal));Greal=zero
  allocate(Gtest(Lmats));Gtest=zero
  allocate(dens(Norb));dens=zero
  allocate(SigmaBHZ(Nso,Nso));SigmaBHZ=zero
  allocate(Zmats(Nso,Nso));Zmats=zero
  
  !call set_sigmaBHZ()
  call build_hk(trim(hkfile))
  Smats=zero
  call dmft_gloc_matsubara(Hk,Gmats,Smats)
  call dmft_print_gf_matsubara(Gmats,"Gtest",iprint=1)

  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath)
  
  call ed_solve(comm,bath,Hloc=j2so(sg77Hloc))

  call finalize_MPI()



contains


  !---------------------------------------------------------------------
  !PURPOSE: GET SG77 HAMILTONIAN 
  !---------------------------------------------------------------------
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,j,ik=0
    integer                             :: ix,iy
    real(8)                             :: kx,ky    
    integer                             :: iorb,jorb
    integer                             :: isporb,jsporb
    integer                             :: ispin,jspin
    real(8)                             :: foo
    integer                             :: unit
    complex(8),dimension(Nso,Nso,Lmats) :: Gmats
    complex(8),dimension(Nso,Nso,Lreal) :: Greal
    real(8)                             :: wm(Lmats),wr(Lreal),dw
    real(8),dimension(Nk*Nk*Nk,3)       :: kgrid ![Nk][Ndim]
    integer,dimension(3)                :: Nkvec
    !
    call TB_set_bk(bkx=[pi2,0d0,0d0],bky=[0d0,pi2,0d0],bkz=[0d0,0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) for SG77:"
    Lk=Nk**3
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    !
    Nkvec = [Nk,Nk,Nk]
    call TB_build_kgrid(Nkvec,kgrid)
    do ik=1,Nk*Nk*Nk
      Hk(:,:,ik) = hk_sg77(kgrid(ik,:),Nso)
    enddo
    allocate(sg77Hloc(Nso,Nso))
    sg77Hloc = zero
    sg77Hloc = sum(Hk,dim=3)/Lk
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !SG77 HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_sg77(kvec,N) result(hk)
   integer                   :: N
   real(8),dimension(:)      :: kvec
   complex(8),dimension(N,N) :: hk
   real(8)                   :: kx,ky,kz
   complex(8)                :: hz,h12,hzm,h12m
   integer                   :: ii
   if(N/=Nso)stop "hk_bhz error: N != Nspin*Norb == 4"
   kx=kvec(1)
   ky=kvec(2)
   kz=kvec(3)

   !hz = 1d0*(cos(kx)-cos(ky))+1d0*sin(kx)*sin(ky)
   !hzm = 1d0*(cos(kx)-cos(ky))+1d0*sin(-kx)*sin(-ky)
   !h12 = 0.5*(cos(kx)+exp((0,1)*kz)*cos(ky)) +(1+exp((0,1)*kz))
   !h12m = conjg(0.5*(cos(kx)+exp(-(0,1)*kz)*cos(ky)) +(1+exp(-(0,1)*kz)))
   !hk(1,1)=hz 
   !hk(2,2)=-hz 
   !hk(1,2)=h12
   !hk(2,1) = conjg(h12)
   !hk(3,3)= hzm
   !hk(4,4)=-hzm
   !hk(3,4) = h12m
   !hk(4,3) = conjg(h12m)
   
   Hk=(cos(kx)+cos(ky)+cos(kz))*zeye(N)
  
 end function hk_sg77


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


end program ed_sg77
