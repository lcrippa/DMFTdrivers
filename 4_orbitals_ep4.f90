program ed_bhz
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  integer                                     :: iloop,Lk,Nso
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
  complex(8),allocatable                      :: Hk(:,:,:),bhzHloc(:,:),sigmaBHZ(:,:),Zmats(:,:),impRho(:,:)
  real(8),allocatable                         :: Wtk(:)
  real(8),allocatable                         :: kxgrid(:),kygrid(:)
  integer,allocatable                         :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: mh,lambda,wmixing,akrange,rh
  character(len=30)                           :: Params
  character(len=16)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,usez,mixG0
  real(8),allocatable                                                    :: wt(:),wr(:)
  complex(8),allocatable                                                 :: wm(:)
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
  logical                                     :: master,getbands

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  size2 = get_Size_MPI(comm)
  master = get_Master_MPI(comm)
    
  !Parse additional variables && read Input && read H(k)^4x4
  call parse_cmd_variable(finput,"FINPUT",default='inputED_BHZ.in')  
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

  if(.not.allocated(wm))allocate(wm(Lmats))
  if(.not.allocated(wr))allocate(wr(Lreal))
  wm     = xi*pi/beta*real(2*arange(1,Lmats)-1,8)
  wr     = linspace(wini,wfin,Lreal)


  Nso=Nspin*Norb

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  allocate(SigmaBHZ(Nso,Nso))
  allocate(Zmats(Nso,Nso))


  !Buil the Hamiltonian on a grid or on  path
  !call build_hk(trim(hkfile))
  call build_hk_GXMG()
  STOP
  greal=zero
  sreal=zero
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)
  call g0test()

  
  
  Nb=ed_get_bath_dimension()
  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath,Hloc=j2so(bhzHloc))
  !
    !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath)
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)


     !Get GLOC:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)


     !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,j2so(bhzHloc),cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=print_mode)


     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif

     !Fit the new bath, starting from the old bath + the supplied Weiss
     select case(ed_mode)
     case default
        stop "ed_mode!=Normal/Nonsu2"
     case("normal")
        call ed_chi2_fitgf(comm,Weiss,bath,ispin=1)
        call ed_spin_symmetrize_bath(bath,save=.true.)
        !else
          !call ed_spin_symmetrize_bath(bath,save=.true.)
        !endif
     case("nonsu2")
        call ed_chi2_fitgf(comm,Weiss,bath)
     end select

     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !

     !Check convergence (if required change chemical potential)
     ! Gtest=zero
     ! do ispin=1,Nspin
     !    do iorb=1,Norb
     !       Gtest=Gtest+Weiss(ispin,ispin,iorb,iorb,:)/Norb/Nspin
     !    enddo
     ! enddo
     Gtest=Weiss(1,1,1,1,:)
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo


  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)


contains


  !---------------------------------------------------------------------
  !PURPOSE: GET THE BHZ HAMILTONIAN ALONG THE Gamma-X-M-Gamma path
  !---------------------------------------------------------------------
  subroutine build_hk_GXMG(kpath_)
    integer                            :: i,j
    integer                            :: Npts,Lkpath
    real(8),dimension(:,:),optional        :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    character(len=64)                      :: file
    !This routine build the H(k) along the GXMG path in BZ,
    !Hk(k) is constructed along this path.
    call TB_set_bk(bkx=[pi2,0d0,0d0],bky=[0d0,pi2,0d0],bkz=[0d0,0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) BHZ along the path GXMG:"
    Npts = 8
    Lkpath=(Npts-1)*Nkpath
    allocate(kpath(Npts,3))
    kpath(1,:)=[0,1,0]!X
    kpath(2,:)=[0,0,0]!G
    kpath(3,:)=[1,1,0]!M
    kpath(4,:)=[1,1,1]!R
    kpath(5,:)=[0,0,1]!Z
    kpath(6,:)=[1,0,1]!A
    kpath(7,:)=[0,0,0]!G
    kpath(8,:)=[0,0,1]!Z
    kpath=kpath*pi
    file="Eigenbands.nint"
    !
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lkpath))
    allocate(wtk(Lkpath))
    call TB_build_model(Hk,hk_bhz,Nso,kpath,Nkpath)
    wtk = 1d0/Lkpath
    if(master)  call TB_Solve_model(hk_bhz,Nso,kpath,Nkpath,&
         colors_name=[orange,blue,red,darkgreen],&
         points_name=[character(len=20) :: "X","G","M","R","Z","A","G","Z"],&
         file=reg(file))
  end subroutine build_hk_GXMG
  
  
  subroutine build_hk(file)
    character(len=*),optional           :: file
    integer                             :: i,ii,j,ik=0
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
    real(8),dimension(Nk**3,3)          :: kgrid
    !
    call TB_set_bk(bkx=[pi2,0d0,0d0],bky=[0d0,pi2,0d0],bkz=[0d0,0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) for BHZ:"
    Lk=Nk**3
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    allocate(wtk(Lk))
    call TB_build_kgrid([Nk,Nk,Nk],kgrid)
    do ii=1,Lk
      kgrid(ii,:)=kgrid(ii,:)-[pi,pi,pi]
    enddo
    !
    call TB_build_model(Hk,hk_bhz,Nso,kgrid)
    wtk = 1d0/Lk
    if(master.AND.present(file))then
       call TB_write_hk(Hk,trim(file),&
            Nlat=1,&
            Nspin=1,&
            Norb=Norb,&
            Nkvec=[Nk,Nk,Nk])
    endif
    allocate(bhzHloc(Nso,Nso))
    bhzHloc = zero
    bhzHloc = sum(Hk,dim=3)/Lk
    where(abs(dreal(bhzHloc))<1d-6)bhzHloc=zero
    if(master)  call TB_write_Hloc(bhzHloc)
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !BHZ HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_bhz(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky,kz
    integer                   :: ii
    kx=kvec(1)
    ky=kvec(2)
    kz=kvec(3)
    Hk(1,:) = [0d0,-cos(kx),0d0,-cos(kz)]
    Hk(2,:) = [-cos(kx),0d0,-cos(ky),0d0]
    Hk(3,:) = [0d0,-cos(ky),0d0, cos(kx)]
    Hk(4,:) = [-cos(kz),0d0, cos(kx),0d0]
  end function hk_bhz


    !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma)
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5);print*,L
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 5. Lmats/Lreal"
    do iorb=1,Norb
       unit=free_unit()
       if(L==Lreal)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_realw.ed"
       elseif(L==Lmats)then
          suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(iorb))//"_iw.ed"
       endif
       write(*,*)"read from file=","impSigma"//reg(suffix)
       open(unit,file="impSigma"//reg(suffix),status='old')
       do i=1,L
          read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
          forall(ispin=1:Nspin)sigma(ispin,ispin,iorb,iorb,i)=dcmplx(reS(ispin),imS(ispin))
       enddo
       close(unit)
    enddo
  end subroutine read_sigma


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
  
  
  
   subroutine g0test() 
      integer                                                     :: ilat,jlat,ispin,iorb,jorb,ii,iw
      integer,dimension(:),allocatable                            :: ind1,ind2
      complex(8),dimension(Nspin*Norb,Nspin*Norb)           :: s_lso
      complex(8),dimension(Nspin*Norb,Nspin*Norb,Lreal)           :: g_lso
      logical                                                     :: wprint_
      real(8),dimension(Nk**3,3)          :: kgrid
      !
      !
      call TB_build_kgrid([Nk,Nk,Nk],kgrid)
      Lk=Nk**3
      print*,Lreal
      do ii=1,Lk
        kgrid(ii,:)=kgrid(ii,:)-[pi,pi,pi]
      enddo
      !
      s_lso=zero
      g_lso=zero
      !
      do iw=1,Lreal
        do ii=1,Lk
          s_lso=(wr(iw)+xi*0.1)*eye(Nso) - Hk_bhz(kgrid(ii,:),4)
          call inv(s_lso)
          G_lso(:,:,iw)=G_lso(:,:,iw)+s_lso/Lk
        enddo
      enddo
      do iorb=1,Nso
        call splot("g0_"//reg(txtfy(iorb-1))//""//reg(txtfy(iorb-1))//".ed",wr,g_lso(iorb,iorb,:))
      enddo
      !   
   end subroutine g0test


end program ed_bhz



