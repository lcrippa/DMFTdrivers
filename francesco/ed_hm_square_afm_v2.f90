program ed_hm_square
  USE EDIPACK2
  USE SCIFOR
  USE DMFT_TOOLS
  USE MPI
  implicit none
  integer                                       :: iloop,Nb,Lk,Nx,Nso,Nlso,ik,iorb,irepl,ispin
  logical                                       :: converged
  real(8)                                       :: wband,wmixing,onsite
  integer                                       :: Nlat,ilat,label,ilab,jlab,li,lj
  real(8)                                       :: ts
  real(8),dimension(:),allocatable              :: dens
  !Bath:
  real(8),allocatable                             :: Bath(:,:),Bath_(:,:)
  real(8),dimension(:,:,:),allocatable            :: lambdasym_vectors
  complex(8),dimension(:,:,:,:,:),allocatable     :: Hsym_basis
  !The local hybridization function:
  complex(8),allocatable,dimension(:,:,:,:,:)     :: Hloc
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Gmats
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Greal
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Smats
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Sreal
  complex(8),allocatable,dimension(:,:,:,:,:,:)   :: Weiss,Weiss_
  complex(8),allocatable,dimension(:,:,:,:,:,:,:) :: Gkmats
  complex(8),allocatable,dimension(:,:)           :: Gtest
  !Luttinger invariants:
  real(8),allocatable                           :: luttinger(:)
  !
  character(len=16)                             :: finput,foutput
  complex(8),allocatable                        :: Hk(:,:,:)
  !
  integer                                       :: comm,rank
  logical                                       :: master
  logical                                       :: mixG0,symOrbs
  !
  real(8)                                       :: Vel
  !
  real(8),dimension(:,:),allocatable            :: kpath
  integer                                       :: Npts,Nkpath
  type(rgb_color),dimension(:),allocatable      :: colors

  call init_MPI(comm,.true.)
  rank = get_Rank_MPI(comm)
  master = get_Master_MPI(comm)

  Nlat = 4                       !two sub-lattices

  call parse_cmd_variable(finput,"FINPUT",default='inputHM.in')
  call parse_input_variable(wmixing,"wmixing",finput,default=0.5d0,comment="Mixing bath parameter")
  call parse_input_variable(ts,"TS",finput,default=0.25d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=100,comment="Number of kx point for 2d BZ integration")
  call parse_input_variable(mixG0,"mixG0",finput,default=.false.)
  call parse_input_variable(symOrbs,"symOrbs",finput,default=.false.)
  call parse_input_variable(Vel,"Vel",finput,default=0d0,comment="hopping parameter")
  !
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(beta,"BETA")
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")
  call add_ctrl_var(nbath,"nbath")
  call add_ctrl_var(ed_hw_bath,"ed_hw_bath")

  if(Nspin/=2.OR.Norb>5)stop "Wrong setup from input file: Nspin/=2 OR Norb>5"
  Nso=Nspin*Norb
  Nlso=Nlat*Nspin*Norb

  !Allocate Fields:
  allocate(Weiss(Nlat,Nspin,Nspin,Norb,Norb,Lmats),Weiss_(Nlat,Nspin,Nspin,Norb,Norb,Lmats))
  allocate(Gmats(Nlat,Nspin,Nspin,Norb,Norb,Lmats),Greal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Smats(Nlat,Nspin,Nspin,Norb,Norb,Lmats),Sreal(Nlat,Nspin,Nspin,Norb,Norb,Lreal))
  allocate(Hloc(Nlat,Nspin,Nspin,Norb,Norb))
  allocate(dens(Norb))
  allocate(Gtest(Nlat,Lmats))

  !Build Hk
  call TB_set_bk(bkx=[pi2/2d0,0d0],bky=[0d0,pi2/2d0])
  Lk = Nx*Nx
  allocate(Hk(Nlso,Nlso,Lk))
  call TB_build_model(Hk(:,:,:),hk_model,Nlso,[Nx,Nx])
  Hloc   = zero
  Hloc = lso2nnn(sum(Hk,dim=3)/Lk)
  where(abs(dreal(Hloc))<1.d-6) Hloc=0d0


  !Printing the bands
  Npts=4
  Nkpath=500
  allocate(kpath(Npts,3))
  kpath(1,:)=kpoint_Gamma/2d0
  kpath(2,:)=kpoint_X1/2d0
  kpath(3,:)=kpoint_M1/2d0
  kpath(4,:)=kpoint_Gamma/2d0
  call TB_solve_model(Hk_model,Nlso,kpath,Nkpath,&
      colors_name=[red1,blue1,red1,blue1,red1,blue1,red1,blue1],&
      points_name=[character(len=20) :: 'G', 'X', 'M', 'G'],&
      file="bands_afm.nint")
  
  if(master)call TB_write_hk(Hk(:,:,:),"Hk2d.dat",Nlat=1,&
                             Nspin=Nspin,&
                             Norb=Norb,&
                             Nkvec=[Nx,Nx])

  !STOP

  !Setup Bath
  select case(bath_type)
      case default
         stop "Wrong setup from input file: bath_type has to be 'normal' or 'replica'"
      case("normal")
         !
         Nb=ed_get_bath_dimension()
         !
      case("replica")
         !
         allocate(lambdasym_vectors(Nlat,Nbath,1)) !Nsym=1
         allocate(Hsym_basis(Nspin,Nspin,Norb,Norb,1))
         !
         Hsym_basis(:,:,:,:,1) = so2nn(zeye(Nso)) !Replica onsite energy
         !
         write(*,*) "Replica initialization: ed_hw_bath="//str(ed_hw_bath)
         !
         do ilat=1,Nlat
           do irepl=1,Nbath
              onsite = irepl -1 - (Nbath-1)/2d0        ![-(Nbath-1)/2:(Nbath-1)/2]
              onsite = onsite * 2*ed_hw_bath/(Nbath-1) !P-H symmetric band, -ed_hw_bath:ed_hw_bath
              lambdasym_vectors(ilat,irepl,1) = onsite      !Multiplies the suitable identity 
           enddo
           !
           if(mod(Nbath,2)==0)then
              lambdasym_vectors(ilat,Nbath/2,1) = -1d-1    !Much needed small energies around
              lambdasym_vectors(ilat,Nbath/2+1,1) = 1d-1   !the fermi level. (for even Nbath)
           endif
         enddo
         !
         call ed_set_Hreplica(Hsym_basis,lambdasym_vectors)
         Nb=ed_get_bath_dimension()
         !
  end select
  !
  allocate(bath(Nlat,Nb))
  allocate(bath_(Nlat,Nb))
  bath_ = zero
  !
  !Setup solver
  call ed_set_hloc(Hloc,Nlat)
  call ed_init_solver(bath)
  
  !li=int(sqrt(real(Nlat)))
  !lj=ceiling(real(Nlat)/real(li))
  !label=Nlat
  !do ilab=1,li
  !   do jlab=lj,1,-1
  !      call ed_break_symmetry_bath(Bath(label,:),sb_field,(-1d0)**(ilab+jlab))
  !      label=label-1
  !   enddo
  !enddo

  !Break symmetry
  call ed_break_symmetry_bath(Bath(1,:),sb_field,-1d0)
  call ed_break_symmetry_bath(Bath(2,:),sb_field,1d0)
  call ed_break_symmetry_bath(Bath(3,:),sb_field,1d0)
  call ed_break_symmetry_bath(Bath(4,:),sb_field,-1d0)


  !DMFT loop
  iloop=0;converged=.false.
  do while(.not.converged.AND.iloop<nloop)
     iloop=iloop+1
     call start_loop(iloop,nloop,"DMFT-loop")

     !Solve the EFFECTIVE IMPURITY PROBLEM (first w/ a guess for the bath)
     call ed_solve(bath) 
     call ed_get_sigma(Smats,Nlat,axis='mats')
     call ed_get_sigma(Sreal,Nlat,axis='real')
     
     !Compute the local gfs on the imaginary axis:
     call dmft_gloc_matsubara(Hk,Gmats,Smats)
     do ilat=1,Nlat
      call dmft_write_gf(Gmats(ilat,:,:,:,:,:),"Gloc_"//trim(str(ilat))//"",axis='mats',iprint=3)
     enddo

     !Get the Weiss field/Delta function to be fitted
     if(cg_scheme=='delta')then        
        call dmft_self_consistency(Gmats,Smats,Weiss,Hloc)
     else
        call dmft_self_consistency(Gmats,Smats,Weiss)
     endif
     
     do ilat=1,Nlat
      call dmft_write_gf(Weiss(ilat,:,:,:,:,:),"Weiss_site"//trim(str(ilat))//"",axis='mats',iprint=3)
     enddo
     
     if(mixG0)then
        if(iloop>1)Weiss = wmixing*Weiss + (1.d0-wmixing)*Weiss_
        Weiss_=Weiss
     endif
     
     !Perform the SELF-CONSISTENCY by fitting the new bath
     call ed_chi2_fitgf(Weiss,bath,ispin=1)
     call ed_chi2_fitgf(Weiss,bath,ispin=Nspin)



     !MIXING:
     if(.not.mixG0)then
        if(iloop>1)Bath = wmixing*Bath + (1.d0-wmixing)*Bath_
        Bath_=Bath
     endif

     !Check convergence (if required change chemical potential)     
     Gtest=zero
     do ilat=1,Nlat
       do ispin=1,Nspin
         do iorb=1,Norb
            Gtest(ilat,:) = Gtest(ilat,:) + Weiss(ilat,ispin,ispin,iorb,iorb,:)/(Nspin*Norb)
         enddo
       enddo
     enddo
     converged = check_convergence(Gtest,dmft_error,nsuccess,nloop,reset=.false.)
     call ed_get_dens(dens)
     if(nread/=0d0)call ed_search_variable(xmu,sum(dens),converged)

     call end_loop
  enddo

  !Get kinetic energy:
  call dmft_kinetic_energy(Hk,Smats)

  !Compute the local gfs on the real axis:
  call dmft_gloc_realaxis(Hk,Greal,Sreal)
     do ilat=1,Nlat
      call dmft_write_gf(Greal(ilat,:,:,:,:,:),"Gloc_site"//trim(str(ilat))//"",axis='real',iprint=3)
     enddo


  call finalize_MPI()

contains

  !-------------------------------------------------------------------------------------------
  !PURPOSE:  Hk model for the 2d square lattice
  !-------------------------------------------------------------------------------------------
  function hk_model(kpoint,N) result(hk_out)
    real(8),dimension(:) :: kpoint
    integer              :: N,ih,ilat,jlat,ispin,iorb,distance
    real(8)              :: kx,ky
    complex(8)           :: hk_out(N,N)
    complex(8)           :: hk(Nlat,Nlat,Nspin,Nspin,Norb,Norb)
    kx=kpoint(1)
    ky=kpoint(2)
    Hk = zero
    do ispin=1,Nspin
      do iorb=1,Norb
        Hk(1,2,ispin,ispin,iorb,iorb)=ts*(1 + exp(xi*2*kx))
        Hk(2,1,ispin,ispin,iorb,iorb)=ts*(1 + exp(-xi*2*kx))
        Hk(2,4,ispin,ispin,iorb,iorb)=ts*(1 + exp(xi*2*ky))
        Hk(4,2,ispin,ispin,iorb,iorb)=ts*(1 + exp(-xi*2*ky))
        Hk(4,3,ispin,ispin,iorb,iorb)=ts*(1 + exp(-xi*2*kx))
        Hk(3,4,ispin,ispin,iorb,iorb)=ts*(1 + exp(xi*2*kx))
        Hk(3,1,ispin,ispin,iorb,iorb)=ts*(1 + exp(-xi*2*ky))
        Hk(1,3,ispin,ispin,iorb,iorb)=ts*(1 + exp(xi*2*ky))
      enddo
    enddo
    hk_out=lso2nnn_full(hk)
  end function hk_model


  !+---------------------------------------------------------------------------+
  !PURPOSE :
  ! > lso2nnn: [Nlso,Nlso] matrix -> [Nlat,Nspin,Nspin,Norb,Norb] array
  !+---------------------------------------------------------------------------+
   function lso2nnn(Fin) result(Fout)
      complex(8),dimension(Nlso,Nlso)                       :: Fin
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
   end function lso2nnn
   
   function lso2nnn_full(Fin) result(Fout)
      complex(8),dimension(Nlso,Nlso)                       :: Fout
      complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb) :: Fin
      integer                                               :: iorb,ispin,ilat,jlat,is
      integer                                               :: jorb,jspin,js
      Fout=zero
      do ilat=1,Nlat
        do jlat=1,Nlat
         do ispin=1,Nspin
            do jspin=1,Nspin
               do iorb=1,Norb
                  do jorb=1,Norb
                     is = iorb + (ispin-1)*Norb + (ilat-1)*Norb*Nspin !lattice-spin-orbit stride
                     js = jorb + (jspin-1)*Norb + (jlat-1)*Norb*Nspin !lattice-spin-orbit stride
                     Fout(is,js) = Fin(ilat,jlat,ispin,jspin,iorb,jorb)
                  enddo
               enddo
            enddo
         enddo
         enddo
      enddo
   end function lso2nnn_full

 !+---------------------------------------------------------------------------+
 !PURPOSE :
 ! > nnn2lso: [Nlat,Nspin,Nspin,Norb,Norb] array -> [Nlso,Nlso] matrix
 !+---------------------------------------------------------------------------+
   function nnn2lso(Fin) result(Fout)
      complex(8),dimension(Nlat,Nspin,Nspin,Norb,Norb)      :: Fin
      complex(8),dimension(Nlso,Nlso)                       :: Fout
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
                     Fout(is,js) = Fin(ilat,ispin,jspin,iorb,jorb)
                  enddo
               enddo
            enddo
         enddo
      enddo
   end function nnn2lso
   
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
    

end program ed_hm_square



