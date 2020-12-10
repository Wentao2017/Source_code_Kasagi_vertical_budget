!################################################################################
!This file is part of Incompact3d.
!
!Incompact3d
!Copyright (c) 2012 Eric Lamballais and Sylvain Laizet
!eric.lamballais@univ-poitiers.fr / sylvain.laizet@gmail.com
!
!    Incompact3d is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation.
!
!    Incompact3d is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with the code.  If not, see <http://www.gnu.org/licenses/>.
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!    We kindly request that you cite Incompact3d in your publications and 
!    presentations. The following citations are suggested:
!
!    1-Laizet S. & Lamballais E., 2009, High-order compact schemes for 
!    incompressible flows: a simple and efficient method with the quasi-spectral 
!    accuracy, J. Comp. Phys.,  vol 228 (15), pp 5989-6015
!
!    2-Laizet S. & Li N., 2011, Incompact3d: a powerful tool to tackle turbulence 
!    problems with up to 0(10^5) computational cores, Int. J. of Numerical 
!    Methods in Fluids, vol 67 (11), pp 1735-1757
!################################################################################

PROGRAM incompact3d

USE decomp_2d
USE decomp_2d_poisson
use decomp_2d_io
USE variables
USE param
USE var
USE MPI
USE IBM
USE derivX
USE derivZ

implicit none

integer :: code,nlock,i,j,k,ijk,ii,bcx,bcy,bcz,fh,ierror,nvect1
real(mytype) :: x,y,z,tmp1,phimax,phimin
double precision :: t1,t2
character(len=20) :: filename

TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

CALL MPI_INIT(code)
call decomp_2d_init(nx,ny,nz,p_row,p_col)
!start from 1 == true
call init_coarser_mesh_statS(nstat,nstat,nstat,.true.)
call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.)
call parameter()

call init_variables
!DEVELOPPEMENT A PLACER DANS incompact3d.prm     !Guo
iimplicit=1                                      !Guo 
if (nrank.eq.0) print *,'Parametre implicite : ',iimplicit  !Guo

call schemes()

if (nclx==0) then
   bcx=0
else
   bcx=1
endif
if (ncly==0) then
   bcy=0
else
   bcy=1
endif
if (nclz==0) then
   bcz=0
else
   bcz=1
endif

call decomp_2d_poisson_init(bcx,bcy,bcz)

call decomp_info_init(nxm,nym,nzm,phG)

!if you want to collect 100 snapshots randomly on XXXXX time steps
!call collect_data() !it will generate 100 random time steps

if (ilit==0) call init(ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)  
if (ilit==1) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,0)
call test_speed_min_max(ux1,uy1,uz1)
if (iscalar==1) call test_scalar_min_max(phi1)

!array for stat to zero
umean=0.;vmean=0.;wmean=0.
uumean=0.;vvmean=0.;wwmean=0.
uvmean=0.;uwmean=0.;vwmean=0.
phimean=0.;phiphimean=0.
utmean=0.;vtmean=0.;wtmean=0.
dudx=0.;dudy=0.;dudz=0.                                          !Budget
dvdx=0.;dvdy=0.;dvdz=0.                                          !Budget
dwdx=0.;dwdy=0.;dwdz=0.                                          !Budget
dudxdudx=0.;dudydudy=0.;dudzdudz=0.                              !Budget
dvdxdvdx=0.;dvdydvdy=0.;dvdzdvdz=0.                              !Budget
dwdxdwdx=0.;dwdydwdy=0.;dwdzdwdz=0.                              !Budget
uuvmean=0.;vvvmean=0.                                            !Budget
vwwmean=0.;pmean=0.;pvmean=0.                                            !Budget 
dphidx=0.;dphidxdphidx=0.;dphidy=0.;dphidydphidy=0.
dphidz=0.;dphidzdphidz=0.
dudxdphidx=0.;dudydphidy=0.;dudzdphidz=0.
uphiumean=0.;uphivmean=0.;uphiwmean=0.
dpdx=0.;phidpdx=0.
phidudx=0.;phidudy=0.;phidudz=0.
udphidx=0.;udphidy=0.;udphidz=0.
dudxdvdx=0.;dudydvdy=0.;dudzdvdz=0.
udpdy=0.;dpdy=0.;vdpdx=0.;uvvmean=0.
vttmean=0.;dvdxdphidx=0.
dvdydphidy=0.;dvdzdphidz=0.;phidpdy=0.
vphivmean=0.;vphiwmean=0.;phidvdx=0.
phidvdy=0.;phidvdz=0.;vdphidx=0.
vdphidy=0.;vdphidz=0.

t1 = MPI_WTIME()

!div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
call decomp_info_init(nxm, nym, nzm, ph1)
call decomp_info_init(nxm, ny, nz, ph4)

!gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
call decomp_info_init(nxm, ny, nz, ph2)  
call decomp_info_init(nxm, nym, nz, ph3) 


do itime=ifirst,ilast
   t=(itime-1)*dt
   if (nrank==0) then
      write(*,1001) itime,t
1001  format('Time step =',i7,', Time unit =',F9.3)
   endif
   
   do itr=1,iadvance_time

      if (nclx.eq.2) then
         call inflow (ux1,uy1,uz1,phi1) !X PENCILS
         call outflow(ux1,uy1,uz1,phi1) !X PENCILS 
      endif

!if (itr.eq.1)  then
!      print *, 'nrank, xstart(1:3): ', nrank, xstart(1), xstart(2), xstart(3)    
!      phimax=-1.0d6
!      phimin= 1.0d6
!      do k=1,xsize(3)
!      do j=1,xsize(2)
!      do i=1,xsize(1)
!        if (phi1(i,j,k).gt.phimax) phimax=phi1(i,j,k)
!        if (phi1(i,j,k).lt.phimin) phimin=phi1(i,j,k)
!      enddo
!      enddo
!      enddo
!      print *,'b-nrank:',nrank, phimax, phimin
!endif


     !X-->Y-->Z-->Y-->X
      call convdiff(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)


      do k=1,xsize(3)
      do j=1,xsize(2)
      do i=1,xsize(1)
         ta1(i,j,k)=ta1(i,j,k)+Betaexp*DeltaT*gravity*phi1(xstart(1)-1+i,xstart(2)-1+j,xstart(3)-1+k)
      enddo
      enddo
      enddo


!if (itr.eq.1)  then
!      phimax=-1.0d6
!      phimin= 1.0d6
!      do k=1,xsize(3)
!      do j=1,xsize(2)
!      do i=1,xsize(1)
!        if (phi1(i,j,k).gt.phimax) phimax=phi1(i,j,k)
!        if (phi1(i,j,k).lt.phimin) phimin=phi1(i,j,k)
!      enddo
!      enddo
!      enddo
!      print *,'a-nrank:',nrank, phimax, phimin
!endif

           
      if (iscalar==1) then
         if(iimplicit==0) then                         !Guo
         call scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
              uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,ep1) 
         else                                          !Guo
         call scalarimp(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&    !Guo
              uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3)           !Guo
         endif                                                                 !Guo 
      endif

      !X PENCILS
      if(iimplicit==0) then
          call intt (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1) 
      else ! d2/dy2 implicite                                                  !Guo
          call inttimp (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1,px1,py1,pz1,&  !Guo
               td1,te1,tf1,ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2)                        !Guo
      endif                                                                            !Guo


      call pre_correc(ux1,uy1,uz1)

      if (ivirt==1) then !solid body old school
         !we are in X-pencil
         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
         call body(ux1,uy1,uz1,ep1)
         call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
      endif

      !X-->Y-->Z
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)       

      !POISSON Z-->Z 
      call decomp_2d_poisson_stg(pp3,bcx,bcy,bcz)

      !Z-->Y-->X
      call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
           ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

      !X PENCILS
      call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1) 
      
     !does not matter -->output=DIV U=0 (in dv3)
      call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
           td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,dv3,&
           nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2)

!      if (nrank==0) ux1(12,64,:)=2.
!      if (nrank==1) ux1(12,64,:)=-1.
     !if (nrank==0) ux(12,64,42)=1.5
     !print *,nrank,xstart(2),xend(2),xstart(3),xend(3)

      call test_speed_min_max(ux1,uy1,uz1)
      if (iscalar==1) call test_scalar_min_max(phi1)

   enddo

   if(itime.gt.100000)then
    call STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
                   uvmean,uwmean,vwmean,phiphimean,tmean,utmean,vtmean,wtmean,dudx,dudy,&
                   dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dudxdudx,dudydudy,dudzdudz,&
                   dvdxdvdx,dvdydvdy,dvdzdvdz,dwdxdwdx,dwdydwdy,dwdzdwdz,uuvmean,&
                   vvvmean,vwwmean,pmean,pvmean,nxmsize,nymsize,nzmsize,&
                   phG,ph2,ph3,pp3,dphidx,dphidxdphidx,dphidy,dphidydphidy,dphidz,dphidzdphidz,&
                   dudxdphidx,dudydphidy,dudzdphidz,uphiumean,uphivmean,uphiwmean,phidpdx,dpdx,&
                   phidudx,phidudy,phidudz,udphidx,udphidy,udphidz,dudxdvdx,dudydvdy,dudzdvdz,&
                   udpdy,dpdy,vdpdx,uvvmean,vttmean,dvdxdphidx,dvdydphidy,dvdzdphidz,phidpdy,&
                   vphivmean,vphiwmean,phidvdx,phidvdy,phidvdz,vdphidx,vdphidy,vdphidz)                  !Budget

   endif

   if (mod(itime,isave)==0) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
        px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,1)
     
   if (mod(itime,imodulo)==0) then
      call VISU_INSTA(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
      call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
           ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
   endif
enddo

t2=MPI_WTIME()-t1
call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
                   MPI_COMM_WORLD,code)
if (nrank==0) print *,'time per time_step: ', &
     t1/float(nproc)/(ilast-ifirst+1),' seconds'
if (nrank==0) print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
if (nrank==0) print *,'Mapping p_row*p_col=',p_row,p_col


!call decomp_2d_poisson_finalize
call decomp_2d_finalize
CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
