program ed_ti_slab
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none

  integer                                       :: iloop
  integer                                       :: Nlat,Nlso,Nsites,Nilso
  integer                                       :: ilat,iorb,ispin
  logical                                       :: converged,PBC
  
  !Bath:
  integer                                       :: Nb
  real(8)                                       :: ts
  real(8),allocatable,dimension(:,:)            :: Bath
   
  !hamiltonian input:
  complex(8),allocatable,dimension(:,:,:)       :: Hkr
  complex(8),allocatable,dimension(:,:)         :: tiHloc
  complex(8),allocatable,dimension(:,:,:,:,:)   :: Hloc
  real(8),dimension(:,:),allocatable            :: kpath
  integer                                       :: Nk,Lk,Ly,Nkpath,Nx,Nz
  real(8)                                       :: dens

  !gamma matrices and fields:
  integer                                       :: BAND_CHECK
  complex(8),dimension(4,4)                     :: emat,soxmat,soymat,sozmat
  real(8)                                       :: e0,mh,lambda,dummymag,dummymag_rescale,SCALE_DISPERSION,MASS_OFFSET
  complex(8),dimension(4,4)                     :: NHmat,GapOpeningMat,magmat
  real(8)                                       :: GapOpeningField,NHfield_up,NHfield_dw
  

  !misc
  logical                                       :: tridiag,lrsym,do_realpart,do_eigenbands
  character(len=60)                             :: finput
  complex(8),dimension(:,:,:),allocatable       :: toconverge

  !mpi
  integer                                       :: comm,rank
  logical                                       :: master,DIAG_SIGMA

  call init_MPI()
  comm = MPI_COMM_WORLD
  call StartMsg_MPI(comm)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)
  ! 
  !
  call parse_cmd_variable(finput,"FINPUT",default='inputED_ti_SLAB.conf')
  call parse_cmd_variable(ts,"TS",default=0.5d0)
  call parse_input_variable(nk,"NK",finput,default=100)
  call parse_input_variable(Ly,"Ly",finput,default=20)
  call parse_input_variable(Nkpath,"NKPATH",finput,default=501)
  call parse_input_variable(tridiag,"TRIDIAG",finput,default=.true.)
  call parse_input_variable(mh,"MH",finput,default=1d0)
  call parse_input_variable(lambda,"LAMBDA",finput,default=0.3d0)
  call parse_input_variable(GapOpeningField,"GAPOPENINGFIELD",finput,default=0.0d0)
  call parse_input_variable(NHfield_up,"NHFIELD_UP",finput,default=0.0d0)
  call parse_input_variable(NHfield_dw,"NHFIELD_DW",finput,default=0.0d0)
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
  call parse_input_variable(Nz,"Nz",finput,default=2,comment="Number of cluster sites in z direction")
  call parse_input_variable(e0,"e0",finput,default=1d0)
  call parse_input_variable(PBC,"PBC",finput,default=.false.)
  call parse_input_variable(MASS_OFFSET,"MASS_OFFSET",finput,default=-3.d0)
  call parse_input_variable(SCALE_DISPERSION,"SCALE_DISPERSION",finput,default=0.1d0)
  !
  call ed_read_input(trim(finput),comm)

  !SETUP THE GAMMA MATRICES:
  emat = kron_pauli( pauli_sigma_0, pauli_tau_z)
  !
  soxmat=kron_pauli( pauli_sigma_z, pauli_tau_x)
  soymat=kron_pauli( pauli_sigma_0, pauli_tau_y)
  sozmat=kron_pauli( pauli_sigma_x, pauli_tau_x)
  magmat=kron_pauli( pauli_sigma_z, pauli_tau_0)
  !
  GapOpeningMat=kron_pauli( pauli_sigma_y, pauli_tau_X)
  NHmat =zero
  NHmat(1:2,1:2) = xi * NHFIELD_UP * pauli_tau_0
  NHmat(3:4,3:4) = xi * NHFIELD_DW * pauli_tau_0
  !

  !Add DMFT CTRL Variables:
  call add_ctrl_var(Norb,"norb")
  call add_ctrl_var(Nspin,"nspin")
  call add_ctrl_var(beta,"beta")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,'wini')
  call add_ctrl_var(wfin,'wfin')
  call add_ctrl_var(eps,"eps")

  !set the global number of lattice sites equal to the number of layers along the y-axis
  Nsites = Ly
  Nlat=Nx*Nz
  Nlso=Nlat*Nspin*Norb
  Nilso=Nsites*Nlso
  !
  !Buil the Hamiltonian on a grid or on  path
  !
  call build_eigenbands()
  
  
  call finalize_MPI()

contains


  !+---------------------------------------------------------------------------+!
  !PURPOSE: solve H_ti(k_x,R_y,kz) along the 1d -pi:pi path in the BZ.
  !+---------------------------------------------------------------------------+!
  subroutine build_eigenbands(kpath_)
    real(8),dimension(:,:),optional            :: kpath_
    real(8),dimension(:,:),allocatable         :: kpath
    type(rgb_color),dimension(:,:),allocatable :: colors
    integer                                    :: Npts,iii
    real(8)                                    :: offset
    character(len=64)                          :: file
    !
    !PRINT H(kx,Ry) ALONG A -pi:pi PATH
    if(master)write(LOGfile,*)"Solve H(kx,y,kz) along [-Z:Z]:"
    Npts=5
    allocate(Kpath(Npts,3))
    !
    kpath(1,:)=[0d0,0d0,0d0]*pi
    kpath(2,:)=[1d0/Nx,0d0,0d0]*pi
    kpath(3,:)=[1d0/Nx,0d0,1d0/Nz]*pi
    kpath(4,:)=[0d0,0d0,1d0/Nz]*pi
    kpath(5,:)=[0d0,0d0,0d0]*pi
    !
    file="Eigenbands.nint"
    allocate(colors(Ly,Nlso))
    colors=black
    !do iii=1,Ly
    !  colors(iii,:) = [red1,red1,blue1,blue1]
    !enddo
    call solve_nh_model(ti_edge_model,Nsites,Nlso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "G","X","R","Z","G"],&
         file="Eigenbands.nint",&
         pbc=PBC)
  end subroutine build_eigenbands

   !+------------------------------------------------------------------+
   !PURPOSE  : the BHZ-edge model hamiltonian
   !+------------------------------------------------------------------+
   
 function ti_edge_model(kpoint,Nsites,N,pbc) result(Hrk)
    real(8),dimension(:)                    :: kpoint
    real(8)                                 :: kx
    integer                                 :: Nsites,N
    complex(8),dimension(N,N)               :: Hmat,Tmat,TmatH
    complex(8),dimension(Nsites*N,Nsites*N) :: Hrk
    integer                                 :: i,Idmin,Idmax,Itmin,Itmax
    logical                                 :: pbc
    Hrk=zero
    Hmat=hk_model([kpoint(1),kpoint(3)],N)
    Tmat=t0_rk_bhz(N)
    TmatH=conjg(transpose(Tmat))
    do i=1,Nsites
       Idmin=1+(i-1)*N
       Idmax=      i*N
       Hrk(Idmin:Idmax,Idmin:Idmax)=Hmat 
    enddo
    do i=1,Nsites-1
       Idmin=1 + (i-1)*N
       Idmax=        i*N
       Itmin=1 +     i*N
       Itmax=    (i+1)*N
       Hrk(Idmin:Idmax,Itmin:Itmax)=Tmat
       Hrk(Itmin:Itmax,Idmin:Idmax)=TmatH
    enddo
    if(pbc)then
       Itmin=1+(Nsites-1)*N
       Itmax=0+Nsites*N
       Hrk(1:N,Itmin:Itmax)=TmatH
       Hrk(Itmin:Itmax,1:N)=Tmat
    endif
  end function ti_edge_model


  function t0_rk_bhz(N) result(H)
    integer                    :: N,i,Idmin,Idmax
    complex(8),dimension(N,N)  :: H
    if (N/=Nlso) stop "t0_rk_bhz: wrong dimension, the block has to be Nlso"
    H=zero
    do i=1,N/2
       Idmin=2*i-1
       Idmax=2*i
       H(Idmin:Idmax,Idmin:Idmax)=t_y(ts,lambda)
    enddo
  end function t0_rk_bhz



   !+------------------------------------------------------------------+
   !PURPOSE  : Hcluster for the 2d BHZ model
   !+------------------------------------------------------------------+


  function Hloc_model(N,Mh_,ts_,lambda_) result (H0)
    integer                                               :: N,ilat,jlat,ispin,iorb,jorb,ind1,ind2
    real(8)                                               :: Mh_,ts_,lambda_
    complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: hopping_matrix
    complex(8),dimension(N,N)                             :: H0
    !
    hopping_matrix=zero
    ! Hopping along x
    do ispin=1,Nspin
       do ilat=1,Nx
          do jlat=1,Nz
             ind1=indices2N([ilat,jlat])
             hopping_matrix(ind1,ind1,ispin,ispin,:,:)= t_m(mh_)
             if(ilat<Nx)then
                ind2=indices2N([ilat+1,jlat])
                hopping_matrix(ind2,ind1,ispin,ispin,:,:)= t_x(ts_,lambda_,ispin)
             endif
             if(ilat>1)then
                ind2=indices2N([ilat-1,jlat])
                hopping_matrix(ind2,ind1,ispin,ispin,:,:)= dconjg(transpose(t_x(ts_,lambda_,ispin)))

             endif
          enddo
       enddo
   enddo
   !Hopping along z
  do ilat=1,Nx
    do jlat=1,Nz
         ind1=indices2N([ilat,jlat])
         hopping_matrix(ind1,ind1,1,1,:,:)= hopping_matrix(ind1,ind1,1,1,:,:) + xi*NHFIELD_UP*pauli_tau_0
         hopping_matrix(ind1,ind1,2,2,:,:)= hopping_matrix(ind1,ind1,2,2,:,:) + xi*NHFIELD_DW*pauli_tau_0
         if(jlat<Nz)then
            ind2=indices2N([ilat,jlat+1])
            hopping_matrix(ind2,ind1,1,1,:,:)= hopping_matrix(ind2,ind1,1,1,:,:)-ts*pauli_sigma_z
            hopping_matrix(ind2,ind1,2,2,:,:)= hopping_matrix(ind2,ind1,2,2,:,:)-ts*pauli_sigma_z
            hopping_matrix(ind2,ind1,1,2,:,:)= hopping_matrix(ind2,ind1,1,2,:,:)+lambda_z(ts_,lambda_)
            hopping_matrix(ind2,ind1,2,1,:,:)= hopping_matrix(ind2,ind1,2,1,:,:)+lambda_z(ts_,lambda_)
         endif
         if(jlat>1)then
            ind2=indices2N([ilat,jlat-1])
            hopping_matrix(ind2,ind1,1,1,:,:)= hopping_matrix(ind2,ind1,1,1,:,:)-ts*pauli_sigma_z
            hopping_matrix(ind2,ind1,2,2,:,:)= hopping_matrix(ind2,ind1,2,2,:,:)-ts*pauli_sigma_z
            hopping_matrix(ind2,ind1,1,2,:,:)= hopping_matrix(ind2,ind1,1,2,:,:)+dconjg(transpose(lambda_z(ts_,lambda_)))
            hopping_matrix(ind2,ind1,2,1,:,:)= hopping_matrix(ind2,ind1,2,1,:,:)+dconjg(transpose(lambda_z(ts_,lambda_)))
         endif
      enddo
   enddo
    !
    H0=nnn2lso(hopping_matrix)
    !
  end function hloc_model



   function hk_model(kpoint,N) result(Hk)
      integer                                                      :: N,ilat,jlat,ispin,iorb,jorb,i,j,ind1,ind2
      real(8),dimension(:)                                         :: kpoint
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)        :: hopping_matrix
      complex(8),dimension(N,N)                                    :: hk
      !
      hopping_matrix=zero
      !
      do ispin=1,Nspin
         do ilat=1,Nz
            ind1=indices2N([1,ilat])
            ind2=indices2N([Nx,ilat])
            hopping_matrix(ind2,ind1,ispin,ispin,:,:)=hopping_matrix(ind2,ind1,ispin,ispin,:,:) + dconjg(transpose(t_x(ts,lambda,ispin)))*exp(xi*kpoint(1)*Nx)
            hopping_matrix(ind1,ind2,ispin,ispin,:,:)=hopping_matrix(ind1,ind2,ispin,ispin,:,:) + t_x(ts,lambda,ispin)*exp(-xi*kpoint(1)*Nx)
         enddo
      enddo
      !
      do ilat =1,Nx
        ind1=indices2N([ilat,1])
        ind2=indices2N([ilat,Nz])
        hopping_matrix(ind2,ind1,1,1,:,:)=hopping_matrix(ind2,ind1,1,1,:,:) -ts*pauli_sigma_z*exp(xi*kpoint(2)*Nz)
        hopping_matrix(ind2,ind1,2,2,:,:)=hopping_matrix(ind2,ind1,2,2,:,:) -ts*pauli_sigma_z*exp(xi*kpoint(2)*Nz)
        hopping_matrix(ind1,ind2,1,1,:,:)=hopping_matrix(ind1,ind2,1,1,:,:) -ts*pauli_sigma_z*exp(-xi*kpoint(2)*Nz)
        hopping_matrix(ind1,ind2,2,2,:,:)=hopping_matrix(ind1,ind2,2,2,:,:) -ts*pauli_sigma_z*exp(-xi*kpoint(2)*Nz)
        !
        hopping_matrix(ind2,ind1,1,2,:,:)=hopping_matrix(ind2,ind1,1,2,:,:) +dconjg(transpose(lambda_z(ts,lambda)))*exp(-xi*kpoint(2)*Nz)
        hopping_matrix(ind2,ind1,2,1,:,:)=hopping_matrix(ind2,ind1,2,1,:,:) +dconjg(transpose(lambda_z(ts,lambda)))*exp(-xi*kpoint(2)*Nz)
        hopping_matrix(ind1,ind2,1,2,:,:)=hopping_matrix(ind1,ind2,1,2,:,:) +lambda_z(ts,lambda)*exp(xi*kpoint(2)*Nz)
        hopping_matrix(ind1,ind2,2,1,:,:)=hopping_matrix(ind1,ind2,2,1,:,:) +lambda_z(ts,lambda)*exp(xi*kpoint(2)*Nz)
      enddo
      !
      Hk=nnn2lso(hopping_matrix)+hloc_model(N,Mh,ts,lambda)
      !
   end function hk_model


   !AUXILLIARY HOPPING MATRIX CONSTRUCTORS

   function t_m(mass) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: mass
      !
      tmpmat=zero
      tmpmat=mass*pauli_sigma_z
      !
   end function t_m

   function t_x(hop1,hop2,spinsign) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2,sz
      integer                         :: spinsign
      !
      tmpmat=zero
      sz=(-1.d0)**(spinsign+1)
      tmpmat=-hop1*pauli_sigma_z+0.5d0*sz*xi*hop2*pauli_sigma_x
      !
   end function t_x
   
   function lambda_z(hop1,hop2) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2,sz
      !
      tmpmat=zero
      tmpmat=0.5d0*xi*hop2*pauli_sigma_x
      !
   end function lambda_z 


   function t_y(hop1,hop2) result(tmpmat)
      complex(8),dimension(Norb,Norb) :: tmpmat
      real(8)                         :: hop1,hop2
      !
      tmpmat=zero
      tmpmat=-hop1*pauli_sigma_z
      tmpmat(1,2)=-hop2*0.5d0
      tmpmat(2,1)=hop2*0.5d0
      !
   end function t_y



   !+------------------------------------------------------------------+
   !PURPOSE  : Auxilliary reshape functions
   !+------------------------------------------------------------------+
    function indices2N(indices) result(N)
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      !
      N=Nx*(indices(2)-1)+indices(1)
    end function indices2N

    function N2indices(N) result(indices)
      integer,dimension(2)         :: indices
      integer                      :: N,i
      !
      indices(1)=mod(N,Nx)
      if(indices(1)==0)then
         indices(1)=Nx
         indices(2)=(N-Nx)/Nx+1
      else
         indices(2)=N/Nx+1
      endif
    end function N2indices

   
   

   function lso2nnn(Hlso) result(Hnnn)
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
      integer                                               :: ilat,jlat
      integer                                               :: iorb,jorb
      integer                                               :: ispin,jspin
      integer                                               :: is,js
      Hnnn=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hnnn(ilat,jlat,ispin,jspin,iorb,jorb) = Hlso(is,js)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function lso2nnn


   function nnn2lso(Hnnn) result(Hlso)
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Hnnn
      complex(8),dimension(Nlat*Nspin*Norb,Nlat*Nspin*Norb) :: Hlso
      integer                                               :: ilat,jlat
      integer                                               :: iorb,jorb
      integer                                               :: ispin,jspin
      integer                                               :: is,js
      Hlso=zero
      do ilat=1,Nlat
         do jlat=1,Nlat
            do ispin=1,Nspin
               do jspin=1,Nspin
                  do iorb=1,Norb
                     do jorb=1,Norb
                        is = iorb + (ilat-1)*Norb + (ispin-1)*Norb*Nlat
                        js = jorb + (jlat-1)*Norb + (jspin-1)*Norb*Nlat
                        Hlso(is,js) = Hnnn(ilat,jlat,ispin,jspin,iorb,jorb)
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso
     
  subroutine solve_nh_model(hkr_model,Nsites,Nlso,kpath,Nkpath,colors_name,points_name,file,pbc,iproject)
    interface 
       function hkr_model(kpoint,Nsites,Nlso,pbc)
         real(8),dimension(:)                    :: kpoint
         integer                                 :: Nsites,Nlso
         logical                                 :: pbc
         complex(8),dimension(Nlso,Nlso)         :: hkr_model
       end function hkr_model
    end interface
    integer                                       :: Nsites,Nilso,Nlso
    real(8),dimension(:,:)                        :: kpath
    integer                                       :: Nkpath,Nktot
    type(rgb_color),dimension(Nsites,Nlso)        :: colors_name
    character(len=*),dimension(size(kpath,1))     :: points_name
    character(len=*),optional                     :: file
    logical,optional                              :: pbc,iproject
    character(len=256)                            :: file_general,file_real,file_imag
    logical                                       :: pbc_,iproject_
    character(len=256)                            :: xtics
    integer                                       :: Npts,Ndim
    integer                                       :: ipts,ik,ic,unit,iorb,ilat,io,nrot,u1,u2
    real(8)                                       :: coeff(Nsites*Nlso),klen,ktics(size(Kpath,1))
    type(rgb_color)                               :: corb(Nsites*Nlso),c(Nsites*Nlso)
    real(8),dimension(size(kpath,2))              :: kstart,kstop,kpoint,kdiff,bk_x,bk_y,bk_z
    complex(8),dimension(Nsites*Nlso,Nsites*Nlso) :: h,evec
    complex(8),dimension(Nsites*Nlso)             :: Eval
    integer,allocatable                           :: labels(:,:)
    real(8),allocatable                           :: kseg(:)
    complex(8),allocatable                        :: Ekval(:,:)
    integer,allocatable                           :: Ekcol(:,:)
    !
    !
    file_general    = "Eigenbands.tb";if(present(file))file_general=file
    file_real=reg(file_general)//".real"
    file_imag=reg(file_general)//".imag"
    iproject_= .false.        ;if(present(iproject))iproject_=iproject
    pbc_     = .true.         ;if(present(pbc))pbc_=pbc
    !
    Nilso  = Nsites*Nlso
    Npts  = size(kpath,1)
    Ndim  = size(kpath,2)
    Nktot = (Npts-1)*Nkpath
    !
    do ilat=1,Nsites
       do io=1,Nlso
          corb(io + (ilat-1)*Nlso) = colors_name(ilat,io)
       enddo
    enddo
    !
    !
    if(master)then
       write(*,*)"Solving model along the path:"
       do ipts=1,Npts
          write(*,"(A,10(A,1x),A1)")"Point"//str(ipts)//": [",(str(kpath(ipts,ic)),ic=1,size(kpath,2)),"]"
       enddo
    endif
    !
    ic=0
    allocate(kseg(Nktot))
    allocate(ekval(Nktot,Nilso))
    allocate(ekcol(Nktot,Nilso))
    klen = 0d0
    call start_timer()
    do ipts=1,Npts-1
       kstart = kpath(ipts,:)
       kstop  = kpath(ipts+1,:)
       kdiff  = (kstop-kstart)/Nkpath
       ktics(ipts)  = klen
       do ik=1,Nkpath
          ic=ic+1
          kpoint = kstart + (ik-1)*kdiff
          h = hkr_model(kpoint,Nsites,Nlso,pbc) - XMU*eye(Nilso)
          call eig(h,Eval,Evec)
          call sort_eigensystem(Eval,Evec)
          call eta(ic,Nktot)
          do io=1,Nilso
             coeff(:)=Evec(:,io)*conjg(Evec(:,io))
             c(io) = coeff.dot.corb
             Ekval(ic,io) = Eval(io)
             Ekcol(ic,io) = rgb(c(io))
          enddo
          kseg(ic) = klen
          klen = klen + sqrt(dot_product(kdiff,kdiff))
       enddo
    enddo
    ktics(Npts) = Kseg(ic-1)
    call stop_timer()
    !
     open(free_unit(unit),file=str(file_real))
     do io=1,Nilso
        do ic=1,Nktot
           write(unit,*)kseg(ic),real(Ekval(ic,io)),Ekcol(ic,io)
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
     write(unit,*)"plot '"//reg(file_real)//"' every :::0 u 1:2:3 w l lw 3 lc rgb variable"
     write(unit,*)"# to print from the i-th to the j-th block use every :::i::j"
     !
     close(unit)
     !
     call system("chmod +x "//reg(file_real)//".gp")
    
  end subroutine solve_nh_model
    
    
  subroutine sort_eigensystem(eigval,eigvec)
      integer                                :: ii
      integer,dimension(Nilso)               :: labels
      real(8),dimension(Nilso)               :: re_eigval
      complex(8),dimension(Nilso)            :: eigval,eigval_tmp
      complex(8),dimension(Nilso,Nilso)      :: eigvec,eigvec_tmp
      !
      labels=zero
      eigval_tmp=zero
      eigvec_tmp=zero
      !
      re_eigval=REAL(eigval)
      call sort_array(re_eigval,labels)
      !
      do ii=1,Nilso
        eigval_tmp(ii)=eigval(labels(ii))
        eigvec_tmp(:,ii)=eigvec(:,labels(ii))
      enddo
      !
      eigvec=eigvec_tmp
      eigval=eigval_tmp
      !
  end subroutine sort_eigensystem
  
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
  
    
end program ed_ti_slab
