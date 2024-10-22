program ed_pyro
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
  complex(8),allocatable                      :: Hk(:,:,:),pyroHloc(:,:),sigmapyro(:,:),Zmats(:,:),impRho(:,:)
  real(8),allocatable                         :: Wtk(:)
  real(8),allocatable                         :: kxgrid(:),kygrid(:)
  real(8),dimension(3,3)                      :: cmat,rotation_matrix
  real(8),dimension(3)                        :: b1,b2,b3
  integer,allocatable                         :: ik2ix(:),ik2iy(:)
  !variables for the model:
  integer                                     :: Nk,Nkpath
  real(8)                                     :: mh,lambda,wmixing,akrange,ts,rh
  character(len=30)                           :: Params
  character(len=16)                           :: finput
  character(len=32)                           :: hkfile
  logical                                     :: spinsym,usez,mixG0
  !
  real(8),dimension(2)                        :: Eout
  real(8),allocatable                         :: dens(:)
  complex(8),dimension(4,4)                   :: Gamma1,Gamma2,Gamma3,Gamma5,GammaN
  complex(8),dimension(4,4)                   :: GammaE0,GammaEx,GammaEy,GammaEz
  complex(8),dimension(4,4)                   :: GammaR0,GammaRx,GammaRy,GammaRz,GammaRtest
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
  call parse_cmd_variable(finput,"FINPUT",default='inputED_pyro.in')  
  call parse_input_variable(hkfile,"HKFILE",finput,default="hkfile.in")
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ts,"TS",finput,default=0.5d0)
  call parse_input_variable(nkpath,"NKPATH",finput,default=500)
  call parse_input_variable(wmixing,"WMIXING",finput,default=0.75d0)
  call parse_input_variable(spinsym,"SPINSYM",finput,default=.true.)
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
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

  !if(Nspin/=2.OR.Norb/=4)stop "Wrong setup from input file"
  Nso=Nspin*Norb
  usez=.false.

  !Allocate Weiss Field:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  allocate(Sigmapyro(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  gamma1=kron_pauli( pauli_sigma_z, pauli_tau_x)
  gamma2=kron_pauli( pauli_sigma_0,-pauli_tau_y)
  gamma3=kron_pauli( pauli_sigma_x, pauli_tau_x)
  gamma5=kron_pauli( pauli_sigma_0, pauli_tau_z)
  gammaN=kron_pauli( pauli_sigma_0, pauli_tau_0)
  !
  gammaE0=kron_pauli( pauli_sigma_0, pauli_tau_x )
  gammaEx=kron_pauli( pauli_sigma_x, pauli_tau_x )
  gammaEy=kron_pauli( pauli_sigma_y, pauli_tau_x )
  gammaEz=kron_pauli( pauli_sigma_z, pauli_tau_x )
  !
  gammaR0=kron_pauli( pauli_sigma_0, pauli_tau_y )
  gammaRx=kron_pauli( pauli_sigma_x, pauli_tau_y )
  gammaRy=kron_pauli( pauli_sigma_y, pauli_tau_y )
  gammaRz=kron_pauli( pauli_sigma_z, pauli_tau_y )
  !
  gammaRtest(1,:)=[zero,  zero,   one,  zero]
  gammaRtest(2,:)=[zero,  zero,  zero,  -one]
  gammaRtest(3,:)=[ one,  zero,  zero,  zero]
  gammaRtest(4,:)=[-one,    xi,  zero,  zero]

  !Buil the Hamiltonian on a grid or on  path
  call set_sigmapyro()
  call build_hk(trim(hkfile))

  print_mode=1
  if(ed_mode=="nonsu2")print_mode=4

  if(getbands)then
     !call read_array("Smats",Smats)
     call solve_hk_topological(so2j(Smats(:,:,:,:,1)))
     call finalize_MPI()
     stop
  endif

  !Setup solver
  if(bath_type=="replica")then
     select case(trim(Params))
     case default
        stop "Params not in [Mh; Ex; EzEx; E0Ex; ExEy; E0Ez; E0EzEx; E0EzExEy]"
     case("Mh")
        allocate(lambdasym_vector(2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaRtest) ;lambdasym_vector(2)=rh
        
     case("Ex")
        allocate(lambdasym_vector(2))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,2))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEx) ;lambdasym_vector(2)=-sb_field

     case("EzEx")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEz) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx) ;lambdasym_vector(3)=-sb_field

     case("E0Ex")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEx) ;lambdasym_vector(3)=-sb_field

     case("ExEy")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaEx) ;lambdasym_vector(2)=-sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEy) ;lambdasym_vector(3)=sb_field

     case("E0Ez")
        allocate(lambdasym_vector(3))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,3))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field

     case("E0EzEx")
        allocate(lambdasym_vector(4))
        allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,4))
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(4)=-sb_field

     case("E0EzExEy")
        Hsym_basis(:,:,:,:,1)=j2so(Gamma5)  ;lambdasym_vector(1)=Mh
        Hsym_basis(:,:,:,:,2)=j2so(GammaE0) ;lambdasym_vector(2)=sb_field
        Hsym_basis(:,:,:,:,3)=j2so(GammaEz) ;lambdasym_vector(3)=sb_field
        Hsym_basis(:,:,:,:,4)=j2so(GammaEx) ;lambdasym_vector(4)=-sb_field
        Hsym_basis(:,:,:,:,5)=j2so(GammaEy) ;lambdasym_vector(5)=sb_field

     end select
     call ed_set_Hreplica(Hsym_basis,lambdasym_vector)
     Nb=ed_get_bath_dimension(Hsym_basis)
  else     
     Nb=ed_get_bath_dimension()
  endif
     allocate(Bath(Nb))
     allocate(Bath_(Nb))
     call ed_init_solver(comm,bath)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,Hloc=j2so(pyroHloc))
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)


     !Get GLOC:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)


     !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,j2so(pyroHloc),cg_scheme)
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
        !if(.not.spinsym)then
           call ed_chi2_fitgf(comm,Weiss,bath,ispin=2)
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

  call dmft_kinetic_energy(Hk,Smats)

  call solve_hk_topological(so2j(Smats(:,:,:,:,1)))

  call save_array("Smats",Smats)
  call save_array("Sreal",Sreal)

  call finalize_MPI()



contains


  !---------------------------------------------------------------------
  !PURPOSE: GET pyro HAMILTONIAN (from the NonInteracting code)
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
    !
    b1=[-1d0,1d0,1d0]
    b2=[1d0,-1d0,1d0]
    b3=[1d0,1d0,-1d0]    
    !
    rotation_matrix=eye(3)
    !
    !rotation_matrix(1,:)=[Sqrt(2d0/3d0),1/Sqrt(6d0), 1/Sqrt(6d0)]
    !rotation_matrix(2,:)=[-(1d0/Sqrt(6d0)),(3d0 + Sqrt(6d0))/6d0,(-3d0 + Sqrt(6d0))/6d0]
    !rotation_matrix(3,:)=[-(1d0/Sqrt(6d0)),(-3d0 + Sqrt(6d0))/6d0,(3d0 + Sqrt(6d0))/6d0]
    !
    call matvec(rotation_matrix,b1)
    call matvec(rotation_matrix,b2)
    call matvec(rotation_matrix,b3)
    !
    call TB_set_bk(bkx=pi2*b1,bky=pi2*b2,bkz=pi2*b3)
    !
    cmat(1,:)=[dot_product(b1,b1),dot_product(b1,b2),dot_product(b1,b3)]
    cmat(2,:)=[dot_product(b2,b1),dot_product(b2,b2),dot_product(b2,b3)]
    cmat(3,:)=[dot_product(b3,b1),dot_product(b3,b2),dot_product(b3,b3)]
    call inv(cmat)
    !
    if(master)write(LOGfile,*)"Build H(k) for pyro:"
    Lk=Nk**3
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    allocate(wtk(Lk))
    !
    call TB_build_model(Hk,hk_pyro,Nso,[Nk,Nk,Nk])
    wtk = 1d0/Lk
    if(master.AND.present(file))then
       call TB_write_hk(Hk,trim(file),&
            Nlat=1,&
            Nspin=1,&
            Norb=Norb,&
            Nkvec=[Nk,Nk,Nk])
    endif
    allocate(pyroHloc(Nso,Nso))
    pyroHloc = zero
    pyroHloc = sum(Hk,dim=3)/Lk
    if(master)  call TB_write_Hloc(pyroHloc)
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_Sigmapyro(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmapyro = zero;if(present(sigma))sigmapyro=sigma
  end subroutine set_Sigmapyro


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
       write(LOGfile,*)"Build H_TOP(k) pyro along path:"
       !
       call set_sigmapyro()
       !
       Npts = 7
       Lk=(Npts-1)*Nkpath
       allocate(kpath(Npts,3))
       kpath(1,:)=[0.00, 0.00, 0.00]
       kpath(2,:)=[1.00, 0.00, 0.00]
       kpath(3,:)=[0.00, 0.00, 0.00]
       kpath(4,:)=[0.00, 1.00, 0.00]
       kpath(5,:)=[0.00, 0.00, 0.00]
       kpath(6,:)=[0.00, 0.00, 1.00]
       kpath(7,:)=[0.00, 0.00, 0.00]
       kpath=kpath*pi2
       !
       call matvec(rotation_matrix,kpath(1,:))
       call matvec(rotation_matrix,kpath(2,:))
       call matvec(rotation_matrix,kpath(3,:))
       call matvec(rotation_matrix,kpath(4,:))
       call matvec(rotation_matrix,kpath(5,:))
       call matvec(rotation_matrix,kpath(6,:))
       call matvec(rotation_matrix,kpath(7,:))
       !
       call set_sigmapyro(sigma)
      sigmapyro=zero
      sigmapyro(1,1)=-1d0*xi
      sigmapyro(2,2)=-1.2d0*xi
      sigmapyro(3,3)=-1.4d0*xi
      sigmapyro(4,4)=-1.6d0*xi
      sigmapyro(5,5)=sigmapyro(1,1)
      sigmapyro(6,6)=sigmapyro(2,2)
      sigmapyro(7,7)=sigmapyro(3,3)
      sigmapyro(8,8)=sigmapyro(4,4)
       !
       !
       call solve_nh_model(hk_pyro,Nso,kpath,Nkpath,&
            colors_name=[red,blue,red,blue],&
            points_name=[character(len=20) :: "G","X","G","Y","G","Z","G"],&
            file="Eig_Htop.ed")
    endif
  end subroutine solve_hk_topological




  !--------------------------------------------------------------------!
  !pyro HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_pyro(kvec,N) result(hk)
    integer                         :: N
    real(8),dimension(:)            :: kvec
    real(8),dimension(3)            :: avec
    complex(8),dimension(N,N)       :: hk
    complex(8),dimension(Norb,Norb) :: hk_block
    real(8)                         :: ek,theta1,theta2,theta3
    integer                         :: ii,jj
    if(N/=Nso)stop "hk_pyro error"
    !
    avec=[dot_product(kvec,b1),dot_product(kvec,b2),dot_product(kvec,b3)]
    !
    theta1=dot_product(cmat(1,:),avec)
    theta2=dot_product(cmat(2,:),avec)
    theta3=dot_product(cmat(3,:),avec)
    !
    Hk = zero
    Hk_block=zero
    !
    Hk_block(1,2) = -ts*(1.d0 + exp(xi*(theta2-theta1)))
    Hk_block(1,3) = -ts*(1.d0 + exp(xi*(theta2-theta3)))
    Hk_block(1,4) = -ts*(1.d0 + exp(xi*theta2))
    Hk_block(2,3) = -ts*(1.d0 + exp(xi*(theta1-theta3)))
    Hk_block(2,4) = -ts*(1.d0 + exp(xi*theta1))
    Hk_block(3,4) = -ts*(1.d0 + exp(xi*theta3))
    !
    do ii = 1,4
      do jj = 1,ii
        Hk_block(ii,jj) = conjg(Hk_block(jj,ii))
      enddo
    enddo
    !
    !
    if(Nspin .eq. 1) then
      Hk=Hk_block
    else
      Hk(1:4,1:4) = Hk_block
      Hk(5:8,5:8) = Hk_block
    endif
    !
    Hk = Hk + Sigmapyro
    if (usez) then
       Zmats=zero
       do ii=1,Nso
          Zmats(ii,ii)  = 1.d0/abs( 1.d0 +  abs(dimag(sigmapyro(ii,ii))/(pi/beta)) )
       end do
       Hk = matmul(Zmats,Hk)
    endif
  end function hk_pyro



  !--------------------------------------------------------------------!
  !TRANSFORMATION BETWEEN DIFFERENT BASIS AND OTHER ROUTINES
  !--------------------------------------------------------------------!
  subroutine read_sigma(sigma) ! TODO: fix this!
    complex(8)        :: sigma(:,:,:,:,:)
    integer           :: iorb,ispin,i,L,unit
    real(8)           :: reS(Nspin),imS(Nspin),ww
    character(len=20) :: suffix
    if(size(sigma,1)/=Nspin)stop "read_sigma: error in dim 1. Nspin"
    if(size(sigma,3)/=Norb)stop "read_sigma: error in dim 3. Norb"
    L=size(sigma,5);print*,L
    if(L/=Lmats.AND.L/=Lreal)stop "read_sigma: error in dim 5. Lmats/Lreal"
    do iorb=1,Norb
      do jorb=1,Norb
         unit=free_unit()
         if(L==Lreal)then
            suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_realw.ed"
         elseif(L==Lmats)then
            suffix="_l"//reg(txtfy(iorb))//"_m"//reg(txtfy(jorb))//"_iw.ed"
         endif
         write(*,*)"read from file=","impSigma"//reg(suffix)
         open(unit,file="impSigma"//reg(suffix),status='old')
         do i=1,L
            read(unit,"(F26.15,6(F26.15))")ww,(imS(ispin),reS(ispin),ispin=1,Nspin)
            forall(ispin=1:Nspin)sigma(ispin,ispin,iorb,jorb,i)=dcmplx(reS(ispin),imS(ispin))
         enddo
         close(unit)
      enddo
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
  
  subroutine matvec(mat,vec)
    real(8),dimension(3,3)    :: mat
    real(8),dimension(3)      :: vec,vec_
    integer                   :: i
    !
    vec_=zero
    do i=1,3
      vec_(i)=dot_product(mat(i,:),vec)
    enddo
    vec=vec_
    !
  end subroutine matvec
   
  
  
  subroutine solve_nh_model(hk_model,Nlso,kpath,Nkpath,colors_name,points_name,file,iproject)
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    integer                                   :: Nlso,io
    real(8),dimension(:,:)                    :: kpath
    integer                                   :: Nkpath
    type(rgb_color),dimension(Nlso)           :: colors_name
    character(len=*),dimension(size(kpath,1)) :: points_name
    character(len=*),optional                 :: file
    logical,optional                          :: iproject
    character(len=256)                        :: file_,file_real,file_imag
    logical                                   :: iproject_
    character(len=256)                        :: xtics
    integer                                   :: Npts,Ndim,Nktot
    integer                                   :: ipts,ik,ic,unit,iorb
    real(8),dimension(size(kpath,2))          :: kstart,kstop,kpoint,kdiff,bk_x,bk_y,bk_z
    real(8)                                   :: coeff(Nlso),klen,ktics(size(Kpath,1))
    complex(8)                                :: h(Nlso,Nlso),evec(Nlso,Nlso),eval(Nlso)
    type(rgb_color)                           :: corb(Nlso),c(Nlso)
    character(len=10)                         :: chpoint
    character(len=32)                         :: fmt
    real(8),allocatable                       :: kseg(:)
    complex(8),allocatable                    :: Ekval(:,:)
    real(8),allocatable                       :: Ekval_sorted(:,:)
    integer,allocatable                       :: Ekcol(:,:)
    !
    master=.true.
    !
    file_    = "Eigenbands.tb";if(present(file))file_=file
    file_real=reg(file_)//".real"
    file_imag=reg(file_)//".imag"
    iproject_= .false.
    !
    Npts = size(kpath,1)
    Ndim = size(kpath,2)
    Nktot= (Npts-1)*Nkpath
    do iorb=1,Nlso
       corb(iorb) = colors_name(iorb)
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
    !
    if(master)then
       write(*,*)"Solving model along the path:"
       write(fmt,"(A3,I0,A)")"(A,",size(kpath,2),"F7.4,A1)"
       do ipts=1,Npts
          write(*,fmt)"Point"//str(ipts)//": [",(kpath(ipts,ic),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic = 0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nlso))
    allocate(ekval_sorted(Nktot,Nlso))
    allocate(ekcol(Nktot,Nlso))
    klen=0d0  
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/dble(Nkpath)
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hk_model(kpoint,Nlso)
          call eig(h,Eval,Evec)
          do iorb=1,Nlso
             coeff(:)=h(:,iorb)*conjg(h(:,iorb))
             c(iorb) = coeff.dot.corb
             Ekval(ic,iorb) = Eval(iorb)
             Ekcol(ic,iorb) = rgb(c(iorb))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = kseg(ic-1)
    !
    if(master)then
       open(free_unit(unit),file=str(file_real))
       do ic=1,Nktot
        Ekval_sorted(ic,:)=lazy_sort(REAL(Ekval(ic,:)))
       enddo
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval_sorted(ic,io),Ekcol(ic,io)
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
       write(unit,*)"plot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_real)//".gp")
    endif
    !
    if(master)then
       open(free_unit(unit),file=str(file_imag))
       do ic=1,Nktot
        Ekval_sorted(ic,:)=lazy_sort(IMAG(Ekval(ic,:)))
       enddo
       do io=1,Nlso
          do ic=1,Nktot
             write(unit,*)kseg(ic),Ekval_sorted(ic,io),Ekcol(ic,io)
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
       write(unit,*)"plot '"//reg(file_imag)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_imag)//".gp")
    endif
  end subroutine solve_nh_model
  
  subroutine bands_3d(hk_model,Nlso,Nk_)
    interface 
       function hk_model(kpoint,N)
         real(8),dimension(:)      :: kpoint
         integer                   :: N
         complex(8),dimension(N,N) :: hk_model
       end function hk_model
    end interface
    integer                                   :: Nlso,io,Nk_
    real(8),dimension(Nk_**2,2)               :: kpath
    type(rgb_color),dimension(Nlso)           :: colors_name
    character(len=256)                        :: file_,file_real,file_imag
    integer                                   :: Ndim,Nktot
    integer                                   :: ik,jk,unit,iorb,ipoint
    real(8),dimension(size(kpath,2))          :: bk_x,bk_y,kpoint
    complex(8)                                :: h(Nlso,Nlso),evec(Nlso,Nlso),eval(Nlso)
    character(len=32)                         :: fmt
    complex(8),allocatable                    :: Ekval(:,:)
    real(8),allocatable                       :: Ekval_sorted(:,:)
    !
    master=.true.
    !
    file_    = "bands_3d.ed"
    file_real=reg(file_)//".real"
    file_imag=reg(file_)//".imag"
    !
    Ndim = size(kpath,2)
    Nktot= Nk_**Ndim
    !
    bk_x=[pi2,0d0]
    bk_y=[0d0,pi2]
    !  
    !
    if(master)then
       write(*,*)"3d bands model:"
    endif
    !
    allocate(ekval(Nktot,Nlso))
    allocate(ekval_sorted(Nktot,Nlso))
    !
    call TB_build_kgrid([Nk_,Nk_],kpath)
    
    do ik=1,Nktot
      kpath(ik,:)=kpath(ik,:)-[pi,pi]
    enddo
    !
    do ik=1,Nktot
       kpoint = kpath(ik,:)
       h = hk_model(kpoint,Nlso)
       call eig(h,Eval,Evec)
       do iorb=1,Nlso
          Ekval(ik,iorb) = Eval(iorb)
       enddo
    enddo
    !
    if(master)then
       open(free_unit(unit),file=str(file_real))
       do ik=1,Nktot
        Ekval_sorted(ik,:)=lazy_sort(REAL(Ekval(ik,:)))
       enddo
       do io=1,2
          do ik=1,Nk_
             do jk=1,Nk_
               ipoint=ik*(Nk_-1)+jk
               write(unit,*)kpath(ipoint,1),kpath(ipoint,2),Ekval_sorted(ipoint,io)
             enddo
             write(unit,*)""
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       !
       open(unit,file=reg(file_real)//".gp")
       write(unit,*)"unset key"
       !
       write(unit,*)"splot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_real)//".gp")
    endif
    !
    if(master)then
       open(free_unit(unit),file=str(file_imag))
       do ik=1,Nktot
        Ekval_sorted(ik,:)=lazy_sort(IMAG(Ekval(ik,:)))
       enddo
       do io=1,Nlso
          do ik=1,Nk_
             do jk=1,Nk_
               ipoint=ik*(Nk_-1)+jk
               write(unit,*)kpath(ipoint,1),kpath(ipoint,2),Ekval_sorted(ipoint,io)
             enddo
             write(unit,*)""
          enddo
          write(unit,*)""
       enddo
       close(unit)
       !
       !
       !
       open(unit,file=reg(file_imag)//".gp")
       write(unit,*)"unset key"
       !
       write(unit,*)"plot '"//reg(file_imag)//"' every :::0 u 1:2:3 w l lw 3 lc rgb 'black'"
       !
       close(unit)
       !
       call system("chmod +x "//reg(file_imag)//".gp")
    endif
  end subroutine bands_3d
  
  
  


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

  



end program ed_pyro



