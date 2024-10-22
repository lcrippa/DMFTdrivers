program ed_sg77
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
  implicit none
  integer                                     :: iloop,Lk,Nso,Nlso,Lakw,ilat,Nlat
  logical                                     :: converged
  !Bath:
  integer                                     :: Nb,iorb,jorb,ispin,print_mode
  real(8),allocatable                         :: Bath(:,:),Bath_(:,:)
  !The local hybridization function:
  complex(8),allocatable                      :: Weiss(:,:,:,:,:,:),Weiss_(:,:,:,:,:,:)
  complex(8),allocatable                      :: Smats(:,:,:,:,:,:),Sreal(:,:,:,:,:,:)
  complex(8),allocatable                      :: Gmats(:,:,:,:,:,:),Greal(:,:,:,:,:,:)
  !hamiltonian input:
  complex(8),allocatable                      :: Hk(:,:,:),sg77Hloc(:,:),sigma_sg77(:,:)
  complex(8),allocatable,dimension(:,:,:,:,:) :: Hloc
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
  logical                                     :: master,getbands,getakw,lrsym

   
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
  call parse_input_variable(nk,"NK",finput,default=50)
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

  !if(Nspin/=2.OR.Norb/=1)stop "Wrong setup from input file: Nspin=Norb=2 -> 4Spin-Orbitals"
   if(Norb/=1)stop "Wrong setup from input file: Norb=1"
   if(Nspin/=1.AND.Nspin/=2)stop "Wrong setup from input file: Nspin=1 or Nspin=2"

  Nlat=2
  Nso=Nspin*Norb
  Nlso=Nlat*Nspin*Norb
  
  !Allocate Weiss Field:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(dens(Norb))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(Sigma_sg77(Nso,Nso))
  
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaR0=kron_pauli( pauli_sigma_0, pauli_tau_y )

 
  call build_hk(trim(hkfile))
  Hloc = lso2nnn(sg77Hloc,Nlat,Nspin,Norb)

 !Setup solver
 
  if(bath_type=="replica")then
      allocate(lambdasym_vector(3))
      allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
      Hsym_basis(:,:,:,:,1)=j2so(GammaE0)  ;lambdasym_vector(1)=1
      Hsym_basis(:,:,:,:,2)=j2so(Gamma5) ;lambdasym_vector(2)=0.5
      Hsym_basis(:,:,:,:,3)=j2so(gammaR0) ;lambdasym_vector(3)=0.75
     call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
     Nb=ed_get_bath_dimension(Hsym_basis)
  else     
     Nb=ed_get_bath_dimension()
  endif
  
  allocate(Bath(Nlat,Nb))
  allocate(Bath_(Nlat,Nb))
  call ed_init_solver(comm,bath)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     call ed_solve(comm,bath,Hloc,mpi_lanc=.false.)
     
     call ed_get_sigma_matsubara(Smats,Nlat)
     call ed_get_sigma_realaxis(Sreal,Nlat)
     call ed_get_dens(dens)

     call dmft_gloc_matsubara(Hk,Gmats,Smats) !tridiag option off
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)


     call dmft_self_consistency(Gmats,Smats,Weiss,Hloc,cg_scheme)
     call dmft_print_gf_matsubara(Weiss,"Weiss",iprint=print_mode)


     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif


     if(spinsym)then
        call ed_chi2_fitgf(Comm,bath,Weiss,Hloc,ispin=1)
        call ed_spin_symmetrize_bath(bath,save=.true.)
     else
        call ed_chi2_fitgf(Comm,bath,Weiss,Hloc,ispin=1)
        call ed_chi2_fitgf(Comm,bath,Weiss,Hloc,ispin=2)
     endif

     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif
     !
     if(master)converged = check_convergence(Weiss(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     call Bcast_MPI(comm,converged)
     call end_loop
  enddo


  call dmft_gloc_realaxis(Hk,Greal,Sreal)
  call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)

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
    integer,dimension(3)     :: Nkvec
    !
    call TB_set_bk(bkx=[pi2,0d0,0d0],bky=[0d0,pi2,0d0],bkz=[0d0,0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) for SG77:"
    Lk=Nk**3
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso

    if(allocated(Hk))deallocate(Hk)
    if(allocated(sg77Hloc))deallocate(sg77Hloc)
    allocate(Hk(Nlso,Nlso,Lk)) ;Hk=zero
    !
    Nkvec = [Nk,Nk,Nk]

    call TB_build_kgrid(Nkvec,kgrid)
    do ik=1,Nk*Nk*Nk
      Hk(:,:,ik) = hk_sg77(kgrid(ik,:),Nlso)
    enddo

    allocate(sg77Hloc(Nlso,Nlso))
    sg77Hloc = zero
    sg77Hloc = sum(Hk,dim=3)/Lk
    where(abs(dreal(sg77Hloc))<1d-6)sg77Hloc=zero
    if(master)  call TB_write_Hloc(sg77Hloc)
    
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_Sigma_sg77(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigma_sg77 = zero;if(present(sigma))sigma_sg77=sigma
  end subroutine set_Sigma_sg77


  !--------------------------------------------------------------------!
  !PURPOSE: Solve the topological Hamiltonian
  !--------------------------------------------------------------------!
  subroutine solve_hk_topological(sigma)
    integer                                :: i,j
    integer                                :: Npts
    complex(8),dimension(Nso,Nso)          :: sigma(Nso,Nso)
    real(8),dimension(:,:),allocatable     :: kpath
    !
    if(master)then
       !This routine build the H(k) along the GXMG path in BZ, Hk(k) is constructed along this path.
       write(LOGfile,*)"Build H_TOP(k) along path:"
       !
       call set_sigma_sg77()
       !
       Npts = 4
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=kpoint_gamma
       kpath(2,:)=kpoint_x1
       kpath(3,:)=kpoint_m1
       kpath(4,:)=kpoint_gamma
       call set_sigma_sg77(sigma)
       call TB_solve_model(hk_sg77,Nso,kpath,Nkpath,&
            colors_name=[red,blue,red,blue],&
            points_name=[character(len=20) :: "{\Symbol G}","X","M","{\Symbol G}"],&
            file="Eig_Htop.ed")
    endif
  end subroutine solve_hk_topological

  !--------------------------------------------------------------------!
  !SG77 HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_sg77(kvec,N) result(hk)
   !! Hamiltonian for SG77.
   integer                   :: N
   real(8),dimension(:)      :: kvec
   complex(8),dimension(N,N) :: hk
   real(8)                   :: kx,ky,kz
   complex(8)                :: hz,h12,hzm,h12m
   integer                   :: ii
   !if(N/=Nso)stop "hk__sg77 error: N != Nspin*Norb == 4"
   kx=kvec(1)
   ky=kvec(2)
   kz=kvec(3)
   if (Nspin == 1) then
      hz = 1d0*(cos(kx)-cos(ky))+1d0*sin(kx)*sin(ky)
      hzm = 1d0*(cos(kx)-cos(ky))+1d0*sin(-kx)*sin(-ky)
      h12 = 0.5*(cos(kx)+exp((0,1)*kz)*cos(ky)) +(1+exp((0,1)*kz))
      h12m = conjg(0.5*(cos(kx)+exp(-(0,1)*kz)*cos(ky)) +(1+exp(-(0,1)*kz)))
      hk(1,1)=hz 
      hk(2,2)=-hz 
      hk(1,2)=h12
      hk(2,1) = conjg(h12)
   endif
   if (Nspin == 2) then
      hz = 1d0*(cos(kx)-cos(ky))+1d0*sin(kx)*sin(ky)
      hzm = 1d0*(cos(kx)-cos(ky))+1d0*sin(-kx)*sin(-ky)
      h12 = 0.5*(cos(kx)+exp((0,1)*kz)*cos(ky)) +(1+exp((0,1)*kz))
      h12m = conjg(0.5*(cos(kx)+exp(-(0,1)*kz)*cos(ky)) +(1+exp(-(0,1)*kz)))
      hk(1,1)=hz 
      hk(3,3)=-hz 
      hk(1,3)=h12
      hk(3,1) = conjg(h12)
      hk(2,2)= hzm
      hk(4,4)=-hzm
      hk(2,4) = h12m
      hk(4,2) = conjg(h12m)
   endif

   !hk(1,1)=hz 
   !hk(2,2)=-hz 
   !hk(1,2)=h12
   !hk(2,1) = conjg(h12)
   !hk(3,3)= hzm
   !hk(4,4)=-hzm
   !hk(3,4) = h12m
   !hk(4,3) = conjg(h12m)
  
 end function hk_sg77












  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!


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



 function lso2nnn(Hlso,Nlat,Nspin,Norb) result(Hnnn)
   !!Convert from 2 dim (Nlat*Nspin*Norb,Nlat*Nspin*Norb) to 5 dim (Nlat,Nspin,Nspin,Norb,Norb) Hamiltonian.
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
                 ! print * , is, js, ilat, ispin, jspin, iorb, jorb
                  Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
               enddo
            enddo
         enddo
      enddo
   enddo
 end function lso2nnn

end program ed_sg77
