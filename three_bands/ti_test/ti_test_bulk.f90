program ed_ti
  USE DMFT_ED
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  USE SF_MPI
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
  complex(8),allocatable                      :: Hk(:,:,:),tiHloc(:,:),sigmati(:,:),Zmats(:,:),impRho(:,:)
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
  complex(8),dimension(4,4)                   :: mat1,mat2,mat3
  complex(8),allocatable,dimension(:,:)       :: matMagX,matMagY,matMagZ,symBrk
  real(8),dimension(3)                        :: mag
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
  call parse_cmd_variable(finput,"FINPUT",default='inputED_ti.in')  
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
  call parse_input_variable(mag,"MAG",finput,default=[0d0,0d0,0d0])
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
  allocate(Weiss(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Weiss_(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Smats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Sreal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Greal(Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Gtest(Lmats))
  allocate(dens(Norb))
  allocate(Sigmati(Nso,Nso))
  allocate(Zmats(Nso,Nso))

  !SETUP THE GAMMA MATRICES:
  allocate(matMagX(Nso,Nso),matMagY(Nso,Nso),matMagZ(Nso,Nso),symBrk(Nso,Nso))
  
  
  !SET 4:
  
  mat1=kron_pauli( pauli_sigma_0, pauli_tau_x)
  mat2=kron_pauli( pauli_sigma_y, pauli_tau_0)
  mat3=kron_pauli( pauli_sigma_x, pauli_tau_0)
 
  matMagX=kron_pauli( pauli_sigma_x, pauli_tau_0 )
  matMagY=kron_pauli( pauli_sigma_y, pauli_tau_0 )
  matMagZ=kron_pauli( pauli_sigma_z, pauli_tau_0 )
  
  !SET 1: the original one (sigma/tau ordered conveniently for the code)
  
  !mat1=kron_pauli( pauli_sigma_0, pauli_tau_x)
  !mat2=kron_pauli( pauli_sigma_x, pauli_tau_y)
  !mat3=kron_pauli( pauli_sigma_y, pauli_tau_y)
  
  !matMagX=kron_pauli( pauli_sigma_x, pauli_tau_z )
  !matMagY=kron_pauli( pauli_sigma_y, pauli_tau_z )
  !matMagZ=kron_pauli( pauli_sigma_z, pauli_tau_z )

  
  !SET 2: problem: perturbation is proportional to tau_z in orbital space, this cannot be a self-energy imaginary part if Ua=Ub
  
  !mat1=kron_pauli( pauli_tau_x, pauli_sigma_0)
  !mat2=kron_pauli( pauli_tau_z, pauli_sigma_x)
  !mat3=kron_pauli( pauli_tau_z, pauli_sigma_y)
  
  !matMagX=kron_pauli( pauli_tau_z, pauli_sigma_x )
  !matMagY=kron_pauli( pauli_tau_z, pauli_sigma_y )
  !matMagZ=kron_pauli( pauli_tau_z, pauli_sigma_z )

  
  !SET 3: problem: in order to have the perturbation proportional to tau_0, we have to have xi*tau0 in the Hamiltonian (which otherwise is not TRS). But this influences the hermiticity of the noninteracting problem...
  
  !mat1=kron_pauli( pauli_sigma_0, pauli_tau_x)
  !mat2=kron_pauli( pauli_sigma_y, xi*pauli_tau_0)
  !mat3=kron_pauli( pauli_sigma_x, xi*pauli_tau_0)
 
  !matMagX=kron_pauli( pauli_sigma_x, pauli_tau_0 )
  !matMagY=kron_pauli( pauli_sigma_y, pauli_tau_0 )
  !matMagZ=kron_pauli( pauli_sigma_z, pauli_tau_0 )
  

  !Buil the Hamiltonian on a grid or on  path
  call set_sigmati()
  call build_hk(trim(hkfile))

  print_mode=1
  if(ed_mode=="nonsu2")print_mode=4

 
  !Setup solver
  Nb=ed_get_bath_dimension()

  allocate(Bath(Nb))
  allocate(Bath_(Nb))
  call ed_init_solver(comm,bath)

  call build_eigenbands()
  STOP

  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(comm,bath,Hloc=j2so(tiHloc))
     call ed_get_sigma_matsubara(Smats)
     call ed_get_sigma_realaxis(Sreal)
     call ed_get_dens(dens)


     !Get GLOC:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     call dmft_print_gf_matsubara(Gmats,"Gloc",iprint=print_mode)


     !Update WeissField:
     call dmft_self_consistency(Gmats,Smats,Weiss,j2so(tiHloc),cg_scheme)
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
        if(.not.spinsym)then
           call ed_chi2_fitgf(comm,Weiss,bath,ispin=2)
        else
           call ed_spin_symmetrize_bath(bath,save=.true.)
        endif
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


  !call dmft_gloc_realaxis(Hk,Greal,Sreal)
  !call dmft_print_gf_realaxis(Greal,"Gloc",iprint=print_mode)

  !call dmft_kinetic_energy(Hk,Smats)

  call build_eigenbands()

  !call save_array("Smats",Smats)
  !call save_array("Sreal",Sreal)

  call finalize_MPI()



contains


  !---------------------------------------------------------------------
  !PURPOSE: GET ti HAMILTONIAN (from the NonInteracting code)
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
    call TB_set_bk(bkx=[pi2,0d0],bky=[0d0,pi2])
    !
    if(master)write(LOGfile,*)"Build H(k) for ti:"
    Lk=Nk**2
    if(master)write(*,*)"# of k-points     :",Lk
    if(master)write(*,*)"# of SO-bands     :",Nso
    if(allocated(Hk))deallocate(Hk)
    if(allocated(wtk))deallocate(wtk)
    allocate(Hk(Nso,Nso,Lk)) ;Hk=zero
    allocate(wtk(Lk))
    !
    call TB_build_model(Hk,hk_ti,Nso,[Nk,Nk])
    wtk = 1d0/Lk
    !
    allocate(tiHloc(Nso,Nso))
    tiHloc = zero
    tiHloc = sum(Hk,dim=3)/Lk
    where(abs(dreal(tiHloc))<1d-6)tiHloc=zero
    if(master)  call TB_write_Hloc(tiHloc)
  end subroutine build_hk






  !--------------------------------------------------------------------!
  !PURPOSE: Set the Self-Energy
  !--------------------------------------------------------------------!
  subroutine set_Sigmati(sigma)
    complex(8),dimension(Nso,Nso),optional :: sigma(Nso,Nso)
    sigmati = zero;if(present(sigma))sigmati=sigma
  end subroutine set_Sigmati




  !--------------------------------------------------------------------!
  !ti HAMILTONIAN:
  !--------------------------------------------------------------------!
  function hk_ti(kvec,N) result(hk)
    integer                   :: N
    real(8),dimension(:)      :: kvec
    complex(8),dimension(N,N) :: hk
    real(8)                   :: ek,kx,ky
    integer                   :: ii
    if(N/=Nso)stop "hk_ti error: N != Nspin*Norb == 4"
    kx=kvec(1)
    ky=kvec(2)
    !
    ek = -1d0*(cos(kx)+cos(ky))
    Hk = (Mh+ek)*mat1 + lambda*sin(ky)*mat2 - lambda*sin(kx)*mat3 + xi*mag(1)*matMagX + xi*mag(2)*matMagY + xi*mag(3)*matMagZ
    !
    Hk = Hk + 0.0*SymBrk
    ! !add the Sigmati term to get Topologial Hamiltonian if required:
    Hk = Hk + dreal(Sigmati)
    if (usez) then
       Zmats=zero
       do ii=1,Nso
          Zmats(ii,ii)  = 1.d0/abs( 1.d0 +  abs(dimag(sigmati(ii,ii))/(pi/beta)) )
       end do
       Hk = matmul(Zmats,Hk)
    endif
  end function hk_ti




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
  
  
  
  subroutine build_eigenbands(kpath_)
    real(8),dimension(:,:),optional    :: kpath_
    real(8),dimension(:,:),allocatable :: kpath
    type(rgb_color),dimension(:),allocatable :: colors
    integer                            :: Npts
    character(len=64)                  :: file
    !
    !PRINT H(kx,Ry) ALONG A -pi:pi PATH
    if(master)write(LOGfile,*)"Solve H(kx,y,kz) along [-Z:Z]:"
    Npts=6
    allocate(Kpath(Npts,3))
    kpath(1,:)=[0,-1,0]*pi
    kpath(2,:)=[0,0,0]*pi
    kpath(3,:)=[0,1,0]*pi
    kpath(4,:)=[1,0,0]*pi
    kpath(5,:)=[0,0,0]*pi
    kpath(6,:)=[-1,0,0]*pi
    file="Eigenbands.nint"
    allocate(colors(Nso))
    colors = [red,blue,red,blue]
    call solve_nh_model(hk_ti,Nso,kpath,Nkpath,&
         colors_name=colors,&
         points_name=[character(len=10) :: "-Y","G","Y","X","G","-X"],&
         file="Eigenbands.nint")
         
    !call bands_3d(hk_ti,Nso,50)
  end subroutine build_eigenbands
  
  
  
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


end program ed_ti
