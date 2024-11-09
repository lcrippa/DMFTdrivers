program realistic_soc
  USE SCIFOR
  USE DMFT_TOOLS
  USE EDIPACK2
  USE MPI
  implicit none
  
  
  integer                                       :: Nb,Nso,Nx,NsymMats,Lk,iloop,iorb,ispin
  logical                                       :: converged
  real(8)                                       :: t0,t1,t2,t3,t4,ecf
  real(8)                                       :: wmixing,lambda

  real(8),dimension(:),allocatable              :: wm,wr
  
  real(8),dimension(:),allocatable              :: dens

  integer                                       :: comm,rank,unit
  logical                                       :: master
  logical                                       :: mixG0
  
  character(len=32)                             :: finput
  
  
  real(8),allocatable                           :: Bath(:),Bath_(:)
  real(8),dimension(:,:),allocatable            :: lambdasym_vectors
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_basis
  complex(8),dimension(:,:,:,:,:),allocatable   :: Hsym_t2g
  
  complex(8),dimension(:,:,:,:),allocatable     :: spinorbit_matrix_t2g_nn
  complex(8),dimension(:,:),allocatable         :: spinorbit_matrix_t2g, u_jbasis
  
  
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats,Greal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats,Sreal
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:)           :: Gtest
 
  complex(8),allocatable,dimension(:,:,:,:)     :: S0
  real(8),allocatable,dimension(:,:)            :: Zmats
 
  complex(8),allocatable                        :: Hloc(:,:,:,:)
  complex(8),allocatable                        :: Hloc_so(:,:)
  complex(8),allocatable                        :: Hk(:,:,:)
  
  complex(8),allocatable                        :: tmpmat(:,:)
  
  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm) 

  
  call parse_cmd_variable(finput,"FINPUT",default='input.in')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(Nx,"Nx",finput,default=100,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(lambda,"lambda",finput,default=0d0,comment="soc parameter")

  call ed_read_input(trim(finput))
  Nso=Nspin*Norb
  
  t0 = 0.42d0
  t1 = 0.17d0
  t2 = 0.30d0
  t3 = 0.03d0
  t4 = 0.04d0
  ecf = 0.d0
 
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(nbath,"nbath")
  call add_ctrl_var(ed_hw_bath,"ed_hw_bath")
  
  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)
  
  !The spin-orbit matrix, initialized by hand

  allocate(spinorbit_matrix_t2g(6,6));spinorbit_matrix_t2g=zero
  allocate(u_jbasis(6,6));u_jbasis=zero
  allocate(spinorbit_matrix_t2g_nn(Nspin,Nspin,3,3));spinorbit_matrix_t2g_nn=zero
  !
  spinorbit_matrix_t2g(1,:) = (-0.5)*[  zero,   -xi,   zero,   zero,  zero,   one]
  spinorbit_matrix_t2g(2,:) = (-0.5)*[    xi,  zero,   zero,   zero,  zero,   -xi]
  spinorbit_matrix_t2g(3,:) = (-0.5)*[  zero,  zero,   zero,   -one,    xi,  zero]
  spinorbit_matrix_t2g(4,:) = (-0.5)*[  zero,  zero,   -one,   zero,    xi,  zero]
  spinorbit_matrix_t2g(5,:) = (-0.5)*[  zero,  zero,    -xi,    -xi,  zero,  zero]  
  spinorbit_matrix_t2g(6,:) = (-0.5)*[   one,   xi,    zero,   zero,  zero,  zero] 
  spinorbit_matrix_t2g_nn = so2nn(spinorbit_matrix_t2g,2,3)
  !
  u_jbasis(1,:) = (1.0/sqrt(6d0))*[ -one*sqrt(3d0),     one,          zero,  zero,           zero, -one*sqrt(2d0)]
  u_jbasis(2,:) = (1.0/sqrt(6d0))*[  -xi*sqrt(3d0),     -xi,          zero,  zero,           zero,   xi*sqrt(2d0)]
  u_jbasis(3,:) = (1.0/sqrt(6d0))*[           zero,    zero,          zero, 2*one, -one*sqrt(2d0),           zero]
  u_jbasis(4,:) = (1.0/sqrt(6d0))*[           zero,    zero, one*sqrt(3d0),  -one, -one*sqrt(2d0),           zero]
  u_jbasis(5,:) = (1.0/sqrt(6d0))*[           zero,    zero, -xi*sqrt(3d0),   -xi,  -xi*sqrt(2d0),           zero]  
  u_jbasis(6,:) = (1.0/sqrt(6d0))*[           zero,   2*one,          zero,  zero,           zero,  one*sqrt(2d0)] 


  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats),Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nspin,Nspin,Norb,Norb))
  allocate(S0(Nspin,Nspin,Norb,Norb))
  allocate(Hloc_so(Nso,Nso))
  allocate(Zmats(Nso,Nso))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  
  Weiss=zero
  Weiss_=zero
  Gmats=zero
  Greal=zero
  Smats=zero
  Sreal=zero
  Hloc=zero
  Gtest=zero
  
  !Build Hk and Hloc for impurity and lattice
  call set_hamiltonian() 
  
  !Setup Bath
  call set_bath()
  
  !Setup solver
  call ed_init_solver(bath)
  
  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
    iloop=iloop+1
    call start_loop(iloop,nloop,"DMFT-loop")
    
    !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
    call ed_solve(bath)

    !Set self-energy matrix and the ones used to calculate Gloc
    call ed_get_sigma(Smats,axis='mats')
    call ed_get_sigma(Sreal,axis='real')
    s0=Smats(:,:,:,:,1)
    
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
    
    call ed_get_dens(dens)
    if(nread/=0d0)call ed_search_chemical_potential(xmu,sum(dens),converged)
    
    call end_loop 
  enddo

  call dmft_kinetic_energy(Hk,Smats)

  call finalize_MPI()
  
contains
    !+---------------------------------------------------------------------------+
    !PURPOSE : model Hamiltonian
    !+---------------------------------------------------------------------------+
    function hk_model(kpoint,N) result(h_k)
      real(8),dimension(:) :: kpoint
      integer              :: N,ih,i,j
      complex(8)           :: h_k(N,N)
      complex(8)           :: tmpmat_nn(2,2,3,3)
      complex(8)           :: hk_t2g(6,6)
      H_k = zero
      hk_t2g = zero
      !
      !the base is ordered as [[yz, xz,xy]_up, [yz, xz,xy]_dw]
      !on-site
      Hk_t2g(1,1) = Hk_t2g(1,1) + ecf !on-site
      Hk_t2g(2,2) = Hk_t2g(2,2) + ecf !on-site
      
      !plus-minus x
      !yz,yz (direction x)
      Hk_t2g(1,1) = Hk_t2g(1,1) - t3*exp(-xi*dot_product(kpoint,[1d0,0d0]))  !right
      Hk_t2g(1,1) = Hk_t2g(1,1) - t3*exp(-xi*dot_product(kpoint,[-1d0,0d0])) !left
      !xz,xz (direction y)
      Hk_t2g(2,2) = Hk_t2g(2,2) - t2*exp(-xi*dot_product(kpoint,[1d0,0d0]))  !right
      Hk_t2g(2,2) = Hk_t2g(2,2) - t2*exp(-xi*dot_product(kpoint,[-1d0,0d0])) !left
      !xy,xy (direction z)
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[1d0,0d0]))  !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[-1d0,0d0])) !left
      !xz,yz
      Hk_t2g(1,2) = Hk_t2g(1,2) + 0d0 !right and left
      Hk_t2g(2,1) = Hk_t2g(2,1) + 0d0 !right and left
      
      !plus-minus y
      !yz,yz (direction x)
      Hk_t2g(1,1) = Hk_t2g(1,1) - t2*exp(-xi*dot_product(kpoint,[0d0,1d0]))  !right
      Hk_t2g(1,1) = Hk_t2g(1,1) - t2*exp(-xi*dot_product(kpoint,[0d0,-1d0])) !left
      !xz,xz (direction y)
      Hk_t2g(2,2) = Hk_t2g(2,2) - t3*exp(-xi*dot_product(kpoint,[0d0,1d0]))  !right
      Hk_t2g(2,2) = Hk_t2g(2,2) - t3*exp(-xi*dot_product(kpoint,[0d0,-1d0])) !left
      !xy,xy (direction z)
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[0d0,1d0]))  !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t0*exp(-xi*dot_product(kpoint,[0d0,-1d0])) !left
      !xz,yz
      Hk_t2g(1,2) = Hk_t2g(1,2) + 0d0
      Hk_t2g(2,1) = Hk_t2g(2,1) + 0d0
      
      !plusx plusy
      !yz,yz (direction x)
      Hk_t2g(1,1) = Hk_t2g(1,1) + 0d0
      Hk_t2g(1,1) = Hk_t2g(1,1) + 0d0
      !xz,xz (direction y)
      Hk_t2g(2,2) = Hk_t2g(2,2) + 0d0
      Hk_t2g(2,2) = Hk_t2g(2,2) + 0d0
      !xy,xy (direction z)
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[1d0,1d0]))   !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[-1d0,-1d0])) !left
      !xz,yz
      Hk_t2g(1,2) = Hk_t2g(1,2) - t4*exp(-xi*dot_product(kpoint,[1d0,1d0]))   !right
      Hk_t2g(1,2) = Hk_t2g(1,2) - t4*exp(-xi*dot_product(kpoint,[-1d0,-1d0])) !left
      Hk_t2g(2,1) = Hk_t2g(2,1) - t4*exp(-xi*dot_product(kpoint,[1d0,1d0]))   !right
      Hk_t2g(2,1) = Hk_t2g(2,1) - t4*exp(-xi*dot_product(kpoint,[-1d0,-1d0])) !left
      
      !plusx minusy
      !yz,yz (direction x)
      Hk_t2g(1,1) = Hk_t2g(1,1) + 0d0
      Hk_t2g(1,1) = Hk_t2g(1,1) + 0d0
      !xz,xz (direction y)
      Hk_t2g(2,2) = Hk_t2g(2,2) + 0d0
      Hk_t2g(2,2) = Hk_t2g(2,2) + 0d0
      !xy,xy (direction z)
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[1d0,-1d0])) !right
      Hk_t2g(3,3) = Hk_t2g(3,3) - t1*exp(-xi*dot_product(kpoint,[-1d0,1d0])) !left
      !xz,yz
      Hk_t2g(1,2) = Hk_t2g(1,2) + t4*exp(-xi*dot_product(kpoint,[1d0,-1d0])) !right
      Hk_t2g(1,2) = Hk_t2g(1,2) + t4*exp(-xi*dot_product(kpoint,[-1d0,1d0])) !left
      Hk_t2g(2,1) = Hk_t2g(2,1) + t4*exp(-xi*dot_product(kpoint,[1d0,-1d0])) !right
      Hk_t2g(2,1) = Hk_t2g(2,1) + t4*exp(-xi*dot_product(kpoint,[-1d0,1d0])) !left
      !
      !
      Hk_t2g(4:6,4:6) = Hk_t2g(1:3,1:3)
      !here add SOC!
      Hk_t2g = Hk_t2g + lambda * spinorbit_matrix_t2g

      tmpmat_nn = so2nn(Hk_t2g,2,3)
      H_k = nn2so(tmpmat_nn(:,:,1:Norb,1:Norb),Nspin,Norb)

    end function hk_model

    !+---------------------------------------------------------------------------+
    !PURPOSE : set Hamiltonian
    !+---------------------------------------------------------------------------+
    
    subroutine set_hamiltonian()
      integer :: i

      call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
      
      Lk = Nx*Nx
      
      allocate(Hk(Nso,Nso,Lk))
      Hk = zero
      
      call TB_build_model(Hk(:,:,:),hk_model,Nso,[Nx,Nx])
            
      Hloc_so   = zero
      Hloc_so = sum(Hk,dim=3)/Lk
      where(abs(Hloc_so)<1.d-10) Hloc_so=0d0
      
      Hloc = so2nn(Hloc_so,nspin,norb)
      call ed_set_hloc(Hloc) 
      
    end subroutine set_hamiltonian    
    
    
    
    !+---------------------------------------------------------------------------+
    !PURPOSE : reshape routines
    !+---------------------------------------------------------------------------+
    function so2nn(Hlso,n_spin,n_orb) result(Hnnn)
      complex(8),dimension(n_spin*N_orb,n_spin*n_orb) :: Hlso
      complex(8),dimension(n_spin,n_spin,n_orb,n_orb) :: Hnnn
      integer                                         :: iorb,jorb,n_spin,n_orb,Nso
      integer                                         :: ispin,jspin
      integer                                         :: is,js
      Hnnn=zero
      
      
      do ispin=1,n_spin
        do jspin=1,n_spin
          do iorb=1,n_orb
            do jorb=1,n_orb
              is = iorb + (ispin-1)*n_orb
              js = jorb + (jspin-1)*n_orb
              Hnnn(ispin,jspin,iorb,jorb) = Hlso(is,js)
            enddo
          enddo
        enddo
      enddo
    end function so2nn

    function nn2so(Hnnn,n_spin,n_orb) result(Hlso)
      complex(8),dimension(n_spin*n_orb,n_spin*n_orb) :: Hlso
      complex(8),dimension(n_spin,n_spin,n_orb,n_orb) :: Hnnn
      integer                                         :: iorb,jorb,n_spin,n_orb
      integer                                         :: ispin,jspin
      integer                                         :: is,js
      
      
      Hlso=zero
      do ispin=1,n_spin
        do jspin=1,n_spin
          do iorb=1,n_orb
            do jorb=1,n_orb
              is = iorb + (ispin-1)*n_orb
              js = jorb + (jspin-1)*n_orb
              Hlso(is,js) = Hnnn(ispin,jspin,iorb,jorb)
            enddo
          enddo
        enddo
      enddo
    end function nn2so

  !+---------------------------------------------------------------------------+
  !set bath
  !+---------------------------------------------------------------------------+
    
  subroutine set_bath()

    if(bath_type =="replica" .or. bath_type =="general")then
      call parametrization_symmetric()

      if(bath_type =="replica")then
        call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
      else
        call ed_set_Hgeneral(Hsym_basis,lambdasym_vectors)
      endif

      Nb=ed_get_bath_dimension(NSymMats)
    else
      Nb=ed_get_bath_dimension()
    endif

    allocate(bath(Nb))
    allocate(bath_(Nb))

    bath  = zero
    bath_ = zero

  end subroutine set_bath
    
  subroutine parametrization_symmetric()
    integer              :: irepl,onsite
        
    allocate(tmpmat(Nso,Nso))
    NsymMats=2

    allocate(lambdasym_vectors(Nbath,NSymMats)); lambdasym_vectors=zero
    allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,NSymMats)); Hsym_basis=zero
    allocate(Hsym_t2g(2,2,3,3,NSymMats)); Hsym_t2g=zero

    !symmetry 1: diagonal
    Hsym_t2g(:,:,:,:,1) = so2nn(zeye(3*Nspin),nspin,3)

    !symmetry 2: soc matrix
    Hsym_t2g(:,:,:,:,2) = spinorbit_matrix_t2g_nn * 2.0

    !cut out the block
    Hsym_basis = Hsym_t2g(:,:,1:Norb,1:Norb,:)

    
    write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)

    do irepl=1,Nbath
      onsite = irepl -1 - (Nbath-1)/2d0          ![-(Nbath-1)/2:(Nbath-1)/2]
      if(Nbath>1) onsite = onsite * 2*ed_hw_bath/(Nbath-1)   !P-H symmetric band, -ed_hw_bath:ed_hw_bath
      lambdasym_vectors(irepl,1) = onsite        !Multiplies the suitable identity
    enddo

    !random values around lambda for off-diagonal components
    call random_init(.true., .true.)
    call RANDOM_NUMBER(lambdasym_vectors(:,2:NSymMats))
    lambdasym_vectors(:,2:NSymMats) = (lambdasym_vectors(:,2:NSymMats) - 0.5d0)*0.1 + lambda

    if(mod(Nbath,2)==0 .and. Nbath>2)then
      lambdasym_vectors(Nbath/2,1) = -1d-1    !Much needed small energies around
      lambdasym_vectors(Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
    endif    

  end subroutine parametrization_symmetric
  
  !+---------------------------------------------------------------------------+
  !get glocal
  !+---------------------------------------------------------------------------+
  subroutine get_glocal()
    !    
    call dmft_get_gloc(Hk,Gmats,Smats,axis='mats')
    call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=6)
    !
    call dmft_get_gloc(Hk,Greal,Sreal,axis='real')
    call dmft_write_gf(Greal,"Gloc",axis='real',iprint=6)    
    !
  end subroutine get_glocal
  
  !+---------------------------------------------------------------------------+
  !get weiss field
  !+---------------------------------------------------------------------------+
  
  subroutine get_weissfield()
    if(cg_scheme=='delta')then
      call dmft_self_consistency(Gmats,Smats,Weiss,Hloc)
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
    
    call ed_chi2_fitgf(Weiss,bath)
    
    select case(bath_type)
      case("normal")
        call ed_ph_symmetrize_bath(bath,save=.true.)
        !call ed_spin_symmetrize_bath(bath,save=.true.)
        !call ed_orb_symmetrize_bath(bath,save=.true.)
      case("hybrid")
        call ed_ph_symmetrize_bath(bath,save=.true.)
        !call ed_spin_symmetrize_bath(bath,save=.true.)
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


  subroutine get_zmats()
    complex(8),dimension(Nso,Nso) :: s0so,Uinv,tmpmat,s0nh
    !
    !
    Zmats=zero
    s0so=nn2so(s0,Nspin,Norb)
    !Uinv=u_jbasis
    !call inv(Uinv)
    
    !if(ZJBASIS)then
    !   s0so=matmul(s0so,u_jbasis)
    !   s0so=matmul(Uinv,s0so)
    !endif
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
  !print matrices 
  !+---------------------------------------------------------------------------+
    
  subroutine print_a_matrix(mat)
    real(8),dimension(:,:)               :: mat
    integer                              :: is,js,NN,unit
    character(len=32)                    :: fmt
    !
    unit=LOGfile
    !
    NN = size(mat,1)
    do is=1,NN
       write(unit,"(12(A1,F8.4,A1,2x))")&
            ('(',mat(is,js),')',js =1,Nso)
    enddo
    write(unit,*)""
    !
  end subroutine print_a_matrix
    
  subroutine print_c_matrix(mat)
    complex(8),dimension(:,:)            :: mat
    integer                              :: is,js,NN,unit
    character(len=32)                    :: fmt
    !
    unit=LOGfile
    !
    NN = size(mat,1)
    do is=1,NN
       write(unit,"(20(A1,F8.4,A1,F8.4,A1,2x))")&
            ('(',real(mat(is,js)),',',imag(mat(is,js)),')',js =1,NN)
    enddo
    write(unit,*)""
    !
  end subroutine print_c_matrix




end program realistic_soc
