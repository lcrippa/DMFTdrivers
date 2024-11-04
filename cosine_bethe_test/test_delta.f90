program soc_test_model
  USE SCIFOR
  USE DMFT_TOOLS
  USE EDIPACK2 
  USE MPI
  implicit none

  integer                                       :: Nb,Nso,le,iorb,ifreq
  real(8),dimension(:),allocatable              :: Wbethe,Dbethe
  real(8),dimension(:,:),allocatable            :: Dbands
  real(8),dimension(:,:),allocatable            :: Ebands
  real(8),dimension(:),allocatable              :: de,dens
  real(8),dimension(:),allocatable              :: Wband,wm,wr
  !Bath:
  real(8),allocatable                           :: Bath(:)
  !The local hybridization function:
  complex(8),dimension(:,:),allocatable         :: Hloc
  real(8),dimension(:),allocatable              :: H0 !elements on the diagonal of Hloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:)       :: Gmats_so
  complex(8),allocatable,dimension(:,:,:)       :: Smats_so
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Weiss,Weiss_
  !
  character(len=32)                             :: finput
  !
  integer                                       :: comm,rank,unit
  logical                                       :: master
  logical                                       :: betheSC

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  
  call parse_cmd_variable(finput,"FINPUT",default='input.in')
  call parse_input_variable(betheSC,"BETHESC",finput,default=.false.)
  call parse_input_variable(Le,"LE",finput,default=500)  
  call ed_read_input(trim(finput))

  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(nbath,"nbath")
  call add_ctrl_var(ed_hw_bath,"ed_hw_bath")


  
  Nso=Nspin*Norb


  allocate(wm(Lmats),wr(Lreal))
  wm = pi/beta*(2*arange(1,Lmats)-1)
  wr = linspace(wini,wfin,Lreal)
  
  allocate(Dbethe(6))
  allocate(Wbethe(6))
  Wbethe=[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]
  Dbethe=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


  !Allocate Fields:
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats_so(Nso,Nso,Lmats))
  allocate(Smats_so(Nso,Nso,Lmats))
  allocate(Hloc(Nso,Nso)) 


  Weiss=zero
  Gmats=zero
  Gmats_so=zero
  Smats=zero
  Smats_so=zero
  Hloc=zero


  !Build Hk and Hloc for impurity and lattice
  call set_dos()


  !Setup Bath
  Nb=ed_get_bath_dimension()
  allocate(bath(Nb))

  !Setup solver
  call ed_init_solver(bath)
  !
  call ed_solve(bath) 

  !Set self-energy matrix and the ones used to calculate Gloc
  call ed_get_sigma(Smats,axis='mats')
  do ifreq=1,Lmats
    Smats_so(:,:,ifreq) = nn2so(Smats(:,:,:,:,ifreq))
  enddo
    
  
  !Compute the local gfs on the imaginary axis:
  call dmft_get_gloc(Ebands,Dbands,H0,Gmats,Smats,axis='mats',diagonal=.true.)
  !call dmft_get_gloc(Ebands,Dbands,H0,Gmats_so,Smats_so,axis='mats',diagonal=.true.)
  !do ifreq=1,Lmats
  !  Gmats(:,:,:,:,ifreq) = so2nn(Gmats_so(:,:,ifreq))
  !enddo
  
  call dmft_write_gf(Gmats,"Gloc",axis='mats',iprint=6)
  call ed_chi2_fitgf(Weiss,bath)

  !Get the Weiss field/Delta function to be fitted
  call get_weissfield()


  call finalize_MPI()

  contains

  
!+---------------------------------------------------------------------------+
!set hamiltonian
!+---------------------------------------------------------------------------+

subroutine set_dos()
  !         
  allocate(Ebands(Nso,Le));Ebands=zero
  allocate(Dbands(Nso,Le));Dbands=zero
  allocate(Wband(Nso));Wband=zero
  allocate(H0(Nso));H0=zero
  allocate(de(Nso))
  Wband = Wbethe(:Nso)  !half bandwidth
  H0    = Dbethe(:Nso)  !crystal field (centri delle dos)
  !
  do iorb=1,Nso
    Ebands(iorb,:) = linspace(-Wband(iorb),Wband(iorb),Le,mesh=de(iorb)) !Ebands is the range of frequencies (energies)
    Dbands(iorb,:) = dens_bethe(Ebands(iorb,:),Wband(iorb))*de(iorb) !Dbands is the dos of the given Ebands
  enddo
  !Set Hloc of the lattice
  Hloc = diag(H0)
         
  call ed_set_Hloc(Hloc) 
  
end subroutine set_dos
        
!+---------------------------------------------------------------------------+
!get weiss field
!+---------------------------------------------------------------------------+

subroutine get_weissfield()
  !
  if(cg_scheme=='delta')then
    if(betheSC)then
      call dmft_get_delta_normal_bethe(Gmats,Weiss,so2nn(Hloc),Wband)
    else
      call dmft_self_consistency(Gmats,Smats,Weiss,so2nn(Hloc))
    endif
    call dmft_write_gf(Weiss,"Delta",axis='mats',iprint=6)
  else
    if(betheSC)then
      call dmft_get_weiss_normal_bethe(Gmats,Weiss,so2nn(Hloc),Wband)
    else
      call dmft_self_consistency(Gmats,Smats,Weiss)
    endif
    call dmft_write_gf(Weiss,"Weiss",axis='mats',iprint=6)
  endif          
end subroutine get_weissfield



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
    !bethe specific
    !+---------------------------------------------------------------------------+


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





