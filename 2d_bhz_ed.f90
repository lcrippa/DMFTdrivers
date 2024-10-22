program test
  USE EDLAT
  USE SCIFOR
  USE DMFT_TOOLS
  implicit none


  character(len=16)    :: finput
  integer              :: Nx,Ny,Nlat
  real(8)              :: ts,Mh,lambda
  call parse_cmd_variable(finput,"FINPUT",default='inputED.conf')
  call parse_input_variable(Mh,"Mh",finput,default=1.d0,comment="crystal field splitting")
  call parse_input_variable(lambda,"lambda",finput,default=0.3d0,comment="spin-orbit coupling")
  call parse_input_variable(ts,"TS",finput,default=-0.5d0,comment="hopping parameter")
  call parse_input_variable(Nx,"Nx",finput,default=2,comment="Number of cluster sites in x direction")
  call parse_input_variable(Ny,"Ny",finput,default=2,comment="Number of cluster sites in y direction")
  call ed_read_input(trim(finput))
  !
  call add_ctrl_var(Norb,"NORB")
  call add_ctrl_var(Nspin,"Nspin")
  call add_ctrl_var(xmu,"xmu")
  call add_ctrl_var(wini,"wini")
  call add_ctrl_var(wfin,"wfin")
  call add_ctrl_var(eps,"eps")


  Nlat=Nx*Ny
  Nsites(:)=Nlat

  call ed_Hij_init(Nsites)
  call set_hloc_model(Mh,ts,lambda)  

  call ed_Hij_info()
  call ed_Hij_write(file="hamiltonian.saved")


  call ed_init_solver()
  call ed_solve()

contains

subroutine set_hloc_model(Mh_,ts_,lambda_)
  integer                                                   :: ilat,jlat,ispin,iorb,jorb,ind1,ind2
  real(8)                                                   :: Mh_,ts_,lambda_
  complex(8),dimension(Nlat,Nlat,Nspin,Nspin,Norb,Norb)     :: hopping_matrix
  !
  hopping_matrix=zero
  !
  do ispin=1,Nspin
     do ilat=1,Nx
        do jlat=1,Ny
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
           if(jlat<Ny)then
              ind2=indices2N([ilat,jlat+1])
              hopping_matrix(ind2,ind1,ispin,ispin,:,:)= t_y(ts_,lambda_)
           endif
           if(jlat>1)then
              ind2=indices2N([ilat,jlat-1])
              hopping_matrix(ind2,ind1,ispin,ispin,:,:)= transpose(t_y(ts_,lambda_))
           endif
        enddo
     enddo
  enddo
  !
  do ispin=1,Nspin
    do ilat=1,Nlat
      do jlat=1,Nlat
        do iorb=1,Norb
          do jorb=1,Norb
            call ed_Hij_add_link(ilat,jlat,iorb,jorb,ispin,hopping_matrix(ilat,jlat,ispin,ispin,iorb,jorb))
          enddo
        enddo
      enddo
    enddo
  enddo   
  !
end subroutine set_hloc_model




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

subroutine naming_convention()
  integer                       :: i,j
  integer,dimension(Nx,Ny)      :: matrix
  !
  do j=1,Ny
     do i=1,Nx
        matrix(i,j)=indices2N([i,j])
     enddo
  enddo
  !
  write(LOGfile,"(A)")"The unique index of each site (on the cartesian plane) is as follows:"
  write(LOGfile,"(A)")" "
  do j=1,Ny
     write(LOGfile,"(20(I2,2x))")(matrix(i,Ny+1-j),i =1,Nx)
  enddo
  write(LOGfile,"(A)")" "
end subroutine naming_convention


end program test

