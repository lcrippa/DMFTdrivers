!   Solve the Hubbard model with AFM 2 atoms in the basis 
program ed_hm_square_afm2
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: ip,iloop,ilat,ineq,Lk,Nso,Nlso,ispin,iorb
  logical                                       :: converged
  integer                                       :: Nineq
  !Bath:
  integer                                       :: Nb
  real(8),allocatable                           :: Bath_ineq(:,:),Bath_prev(:,:)

  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Weiss_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Smats,Smats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Sreal,Sreal_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Gmats,Gmats_ineq
  complex(8),allocatable,dimension(:,:,:,:,:,:) :: Greal,Greal_ineq
  !Hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hk ![Nineq*Nspin*Norb,Nineq*Nspin*Norb,Nk]
  complex(8),allocatable,dimension(:,:)         :: modelHloc ![Nineq*Nspin*Norb,Nineq*Nspin*Norb]
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc_ineq
  !variables for the model:
  character(len=16)                             :: finput
  real(8)                                       :: ts,wmixing
  integer                                       :: Nktot,Nk,Nkpath
  logical                                       :: spinsym,neelsym

  integer                                       :: comm,rank
  logical                                       :: master


  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  call parse_cmd_variable(finput   , "FINPUT" , default='inputED.conf')
  call parse_input_variable(ts     , "TS"     , finput, default=1.d0)
  call parse_input_variable(Nk    , "Nk"    , finput, default=25)
  call parse_input_variable(nkpath , "NKPATH" , finput, default=500)
  call parse_input_variable(wmixing, "WMIXING", finput, default=0.75d0)
  call parse_input_variable(NINEQ, "NINEQ", finput, default=1)

  call ed_read_input(trim(finput),comm)


  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  Nso=Nspin*Norb
  Nlso=Nineq*Nso
  if(Norb/=1)stop  "Norb != 1"
  if(Nspin/=1)stop  "Nspin != 1"
  !
  allocate(Weiss_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gmats_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Greal_ineq(Nineq,Nspin,Nspin,Norb,Norb,Lreal))
  !
  allocate(Hloc_ineq(Nineq,Nspin,Nspin,Norb,Norb))


  call build_hk()
  Hloc_ineq = lso2nnn_reshape(modelHloc,Nineq,Nspin,Norb)

  !Setup solver
  Nb=ed_get_bath_dimension()
  allocate(Bath_ineq(Nineq,Nb))
  allocate(Bath_prev(Nineq,Nb))
  call ed_init_solver(comm,Bath_ineq)

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")
     !
     !solve the impurity problem:
     call ed_solve(Comm,Bath_ineq,Hloc_ineq)
     !
     !retrieve inequivalent self-energies:
     do ilat=1,Nineq
      call ed_get_sigma_matsubara(Smats_ineq(ilat,:,:,:,:,:),Nineq)
     enddo
     !

     !
     ! compute the local gf:
     call dmft_gloc_matsubara(Hk,Gmats_ineq,Smats_ineq)
     call dmft_print_gf_matsubara(Gmats_ineq,"Gloc",iprint=4)
     !
     ! compute the Weiss field (only the Nineq ones)
     call dmft_self_consistency(Gmats_ineq,Smats_ineq,Weiss_ineq,Hloc_ineq,cg_scheme)
     !
     !
     ! fit baths and mix result with old baths
     call ed_chi2_fitgf(Comm,Bath_ineq,Weiss_ineq,Hloc_ineq,ispin=1)
     !
     !
     ! Mixing:
     if(iloop>1)Bath_ineq = wmixing*Bath_ineq + (1.d0-wmixing)*Bath_prev
     Bath_prev=Bath_ineq

     ! Convergence
     converged = check_convergence(Weiss_ineq(:,1,1,1,1,:),dmft_error,nsuccess,nloop)
     !
     call end_loop
  enddo


  call dmft_gloc_realaxis(Hk,Greal_ineq,Sreal_ineq)
  call dmft_print_gf_realaxis(Greal_ineq,"Gloc",iprint=4)

  call finalize_MPI()


contains



  !--------------------------------------------------------------------!
  !PURPOSE: BUILD THE H(k) FOR THE BHZ-AFM MODEL.
  !--------------------------------------------------------------------!
  subroutine build_hk()
    integer                                 :: Npts
    integer                                 :: i,j,k,ik,iorb,jorb
    integer                                 :: ix,iy,iz
    real(8)                                 :: kx,ky,kz
    real(8),dimension(:),allocatable        :: kxgrid,kygrid
    real(8),dimension(:,:),allocatable      :: kpath
    real(8),dimension(2)                    :: bk1,bk2,kvec
    real(8)                                 :: n(Nlso)
    complex(8)                              :: w
    complex(8)                              :: Gmats(Nlso,Nlso,Lmats),Greal(Nlso,Nlso,Lreal)
    complex(8)                              :: Smats(Nlso,Nlso,Lmats),Sreal(Nlso,Nlso,Lreal)
    !
    Nktot=Nk*Nk
    write(LOGfile,*)"Using Nk_total="//txtfy(Nktot)
    !
    !
    !>Reciprocal lattice basis vector  
    bk1= [ pi2,  0d0 ]/Nineq
    bk2= [ 0d0,  pi2 ]
    call TB_set_bk(bk1,bk2)
    !
    !
    !>Get TB Hamiltonian matrix
    allocate(Hk(Nlso,Nlso,Nktot))
    allocate(modelHloc(Nlso,Nlso))
    call TB_build_model(Hk,hk_model,Nlso,[Nk,Nk])
    modelHloc = sum(Hk(:,:,:),dim=3)/Nktot
    where(abs((modelHloc))<1.d-9)modelHloc=0.d0
    !
  end subroutine build_hk




  function hk_model(kpoint,N) result(hk)
    real(8),dimension(:)          :: kpoint
    integer                       :: N,i
    real(8)                       :: kx,ky
    complex(8),dimension(N,N)     :: hk
    complex(8),dimension(N,N)     :: h0,tk
    !
    if(N/=Nlso)stop "hk_model error: N != Nlso" 
    kx = kpoint(1)
    ky= kpoint(2)
    !x hopping
    hk=zero
    if(N>1)then
      do i=1,N-1
        hk(i,i+1) = -ts*one
        hk(i+1,i) = -ts*one
      enddo
    endif
    hk(1,N)=hk(1,N)-ts*exp(xi*N*kx)
    hk(N,1)=hk(N,1)-ts*exp(-xi*N*kx)
    !y hopping
    do i=1,N
      hk(i,i) = hk(i,i)-2*ts*cos(ky)
    enddo
  end function hk_model
  

  function lso2nnn_reshape(Hlso,Nineq,Nspin,Norb) result(Hnnn)
    integer                                                :: Nineq,Nspin,Norb
    complex(8),dimension(Nineq*Nspin*Norb,Nineq*Nspin*Norb) :: Hlso
    complex(8),dimension(Nineq,Nspin,Nspin,Norb,Norb)      :: Hnnn
    integer                                                :: iorb,ispin,ilat,is
    integer                                                :: jorb,jspin,js
    Hnnn=zero
    do ilat=1,Nineq
       do ispin=1,Nspin
          do jspin=1,Nspin
             do iorb=1,Norb
                do jorb=1,Norb
                   is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   js = jorb + (jspin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                   Hnnn(ilat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                enddo
             enddo
          enddo
       enddo
    enddo
  end function lso2nnn_reshape
  
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


end program ed_hm_square_afm2



