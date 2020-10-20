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

!############################################################################
!
subroutine VISU_INSTA (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


!############################################################################
!VORTICITY
di1=0.
do ijk=1,nvect1
   di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
        (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
        (tb1(ijk,1,1)-td1(ijk,1,1))**2)
enddo
uvisu=0.
call fine_to_coarseV(1,di1,uvisu)
990 format('vort',I4.4)
write(filename, 990) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!     1,di1,filename)
!############################################################################

!############################################################################
!VELOCITY
uvisu=0.
call fine_to_coarseV(1,ux1,uvisu)
993 format('ux',I4.4)
      write(filename, 993) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ux1,filename)
uvisu=0.
call fine_to_coarseV(1,uy1,uvisu)
994 format('uy',I4.4)
      write(filename, 994) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,uy1,filename)
uvisu=0.
call fine_to_coarseV(1,uz1,uvisu)
995 format('uz',I4.4)
      write(filename, 995) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,uz1,filename)
!############################################################################

!############################################################################
!PASSIVE SCALAR
if (iscalar==1) then
uvisu=0.
call fine_to_coarseV(1,phi1,uvisu)
996 format('phi',I4.4)
   write(filename, 996) itime/imodulo
   call decomp_2d_write_one(1,uvisu,filename,2)
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!        1,phi1,filename)
endif
!############################################################################

!############################################################################
!PRESSURE
!IT IS IN A SEPARATE SUBROUTINE
!############################################################################
end subroutine VISU_INSTA

!############################################################################
!
subroutine STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
     uvmean,uwmean,vwmean,phiphimean,tmean,utmean,vtmean,dudy,uuvmean,vvvmean,vwwmean,&
     uiuiv,duiuivdy,k,d2kdy2,nxmsize,nymsize,nzmsize,phG,ph2,ph3,pp3)                  !Budget
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

integer :: nxmsize,nymsize,nzmsize                                                   !Budget
TYPE(DECOMP_INFO) :: phG,ph2,ph3                                                     !Budget  
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean,utmean,vtmean,&
                                                   dudy,uuvmean,vvvmean,vwwmean,uiuiv,duiuivdy,k,d2kdy2      !Budget
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: phimean, phiphimean               !Budget
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1, td1                       !Budget
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2, tb2, di2                  !Budget

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3         !Budget
!Z PENCILS NXM NYM NZM-->NXM NYM NZ                                                        !Budget 
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: tc3,dip3   !Budget
!Y PENCILS NXM NYM NZ -->NXM NY NZ                                                         !Budget 
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: td2                      !Budget
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: te2,dip2                !Budget
!X PENCILS NXM NY NZ  -->NX NY NZ                                                          !Budget  
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ti1                                   !Budget   
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tj1,tk1,dip1                         !Budget 


!umean=ux1
call fine_to_coarseS(1,ux1,tmean)
umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)

!dudy                                                                                !Budget
call transpose_x_to_y(ux1,ta2)                                                       !Budget
call dery (tb2,ta2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)           !Budget 
call transpose_y_to_x(tb2,td1)                                                       !Budget
call fine_to_coarseS(1,td1,tmean)                                                    !Budget
dudy(:,:,:)=dudy(:,:,:)+tmean(:,:,:)                                                 !Budget

!vmean=uy1
call fine_to_coarseS(1,uy1,tmean)
vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

!wmean=uz1
call fine_to_coarseS(1,uz1,tmean)
wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

!pmean=tj1                                                                           !Budget 
!WORK Z-PENCILS                                                                      !Budget
call interiz6(tc3,pp3,dip3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&            !Budget
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)         !Budget
!WORK Y-PENCILS                                                                      !Budget
call transpose_z_to_y(tc3,td2,ph3) !nxm nym nz                                       !Budget
call interiy6(te2,td2,dip2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&            !Budget
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)                          !Budget
!WORK X-PENCILS                                                                      !Budget
call transpose_y_to_x(te2,ti1,ph2) !nxm ny nz                                        !Budget
call interi6(tj1,ti1,dip1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&                !Budget     
     nxmsize,xsize(1),xsize(2),xsize(3),1)                                           !Budget
!The pressure field on the main mesh is in tj1                                       !Budget
!PRESSURE                                                                            !Budget
call fine_to_coarseS(1,tj1,tmean)                                                    !Budget
pmean(:,:,:)=pmean(:,:,:)+tmean(:,:,:)                                               !Budget

if (iscalar==1) then
   !phimean=phi1
   call fine_to_coarseS(1,phi1,tmean)
   phimean(:,:,:)=phimean(:,:,:)+tmean(:,:,:)
endif

!uumean=ux1*ux1
ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

!vvmean=uy1*uy1
ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

!wwmean=uz1*uz1
ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

!uvmean=ux1*uy1
ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

!uwmean=ux1*uz1
ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

!vwmean=uy1*uz1
ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

!uuvmean=ux1*ux1*uy1
ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)*uy1(:,:,:)   !Budget
call fine_to_coarseS(1,ta1,tmean)             !Budget
uuvmean(:,:,:)=uuvmean(:,:,:)+tmean(:,:,:)    !Budget  

!vvvmean=uy1*uy1*uy1
ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)*uy1(:,:,:)   !Budget
call fine_to_coarseS(1,ta1,tmean)             !Budget
vvvmean(:,:,:)=vvvmean(:,:,:)+tmean(:,:,:)    !Budget 

!vwwmean=uz1*uz1*uz1
ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)*uz1(:,:,:)   !Budget
call fine_to_coarseS(1,ta1,tmean)             !Budget
vwwmean(:,:,:)=vwwmean(:,:,:)+tmean(:,:,:)    !Budget

!utmean=ux1*phi1
ta1(:,:,:)=ux1(:,:,:)*phi1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
utmean(:,:,:)=utmean(:,:,:)+tmean(:,:,:)

!vtmean=uy1*phi1
ta1(:,:,:)=uy1(:,:,:)*phi1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vtmean(:,:,:)=vtmean(:,:,:)+tmean(:,:,:)

!uiuiv=uuv+vvv+vww=(uuvmean-uumean*vmean-2*uvmean*umean)+(vvvmean-3*vmean*vvmean)+(vwwmean-vmean*wwmean-2*wmean*vwmean)
ta1(:,:,:)=uuvmean(:,:,:)-uumean(:,:,:)*vmean(:,:,:)-2*uvmean(:,:,:)*umean(:,:,:)+ &            !Budget
           vvvmean(:,:,:)-3*vmean(:,:,:)*vvmean(:,:,:)+ &                                       !Budget
           vwwmean(:,:,:)-vmean(:,:,:)*wwmean(:,:,:)-2*wmean(:,:,:)*vwmean(:,:,:)               !Budget 
call fine_to_coarseS(1,ta1,tmean)                                                               !Budget 
uiuiv(:,:,:)=tmean(:,:,:)                                                                       !Budget

!duiuivdy
call transpose_x_to_y(uiuiv,ta2)                                                     !Budget
call dery (tb2,ta2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)           !Budget 
call transpose_y_to_x(tb2,td1)                                                       !Budget
call fine_to_coarseS(1,td1,tmean)                                                    !Budget
duiuivdy(:,:,:)=tmean(:,:,:)                                                         !Budget

!k                                                                                   !Budget
ta1(:,:,:)=0.5*(uumean(:,:,:)+vvmean(:,:,:)+wwmean(:,:,:))                           !Budget
call fine_to_coarseS(1,ta1,tmean)                                                    !Budget 
k(:,:,:)=tmean(:,:,:)                                                                !Budget

!d2kdy2                                                                              !Budget
call transpose_x_to_y(k,ta2)                                                         !Budget  
call deryy(tb2,ta2,di2,sy,sfyp,ssyp,swyp,ysize(1),ysize(2),ysize(3),1)               !Budget
call transpose_y_to_x(tb2,td1)                                                       !Budget
call fine_to_coarseS(1,td1,tmean)                                                    !Budget
d2kdy2(:,:,:)=tmean(:,:,:)                                                           !Budget 


if (iscalar==1) then
   !phiphimean=phi1*phi1
   ta1(:,:,:)=phi1(:,:,:)*phi1(:,:,:)
   call fine_to_coarseS(1,ta1,tmean)
   phiphimean(:,:,:)=phiphimean(:,:,:)+tmean(:,:,:)
endif

!for a verification
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ta1,'compa.dat')

if (mod(itime,isave)==0) then
   call decomp_2d_write_one(1,umean,'umean.dat',1)
   call decomp_2d_write_one(1,vmean,'vmean.dat',1)
   call decomp_2d_write_one(1,wmean,'wmean.dat',1)
   call decomp_2d_write_one(1,uumean,'uumean.dat',1)
   call decomp_2d_write_one(1,vvmean,'vvmean.dat',1)
   call decomp_2d_write_one(1,wwmean,'wwmean.dat',1)
   call decomp_2d_write_one(1,uvmean,'uvmean.dat',1)
   call decomp_2d_write_one(1,uwmean,'uwmean.dat',1)
   call decomp_2d_write_one(1,vwmean,'vwmean.dat',1)
   call decomp_2d_write_one(1,utmean,'utmean.dat',1)
   call decomp_2d_write_one(1,vtmean,'vtmean.dat',1)
   call decomp_2d_write_one(1,dudy,'dudy.dat',1)                   !Budget
   call decomp_2d_write_one(1,uuvmean,'uuvmean.dat',1)             !Budget
   call decomp_2d_write_one(1,vvvmean,'vvvmean.dat',1)             !Budget
   call decomp_2d_write_one(1,vwwmean,'vwwmean.dat',1)             !Budget
   call decomp_2d_write_one(1,duiuivdy,'duiuivdy.dat',1)           !Budget
   call decomp_2d_write_one(1,d2kdy2,'d2kdy2.dat',1)               !Budget

   if (nrank==0) print *,'write stat arrays velocity done!'
   if (iscalar==1) then
      call decomp_2d_write_one(1,phimean,'phimean.dat',1)
      call decomp_2d_write_one(1,phiphimean,'phiphimean.dat',1)
      if (nrank==0) print *,'write stat arrays scalar done!'
   endif
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,1,ux1,'compa.dat')
   
endif

end subroutine STATISTIC

!############################################################################
!
subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

integer :: nxmsize,nymsize,nzmsize
TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1 

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1
!PRESSURE
uvisu=0.
call fine_to_coarseV(1,tb1,uvisu)
990 format('pp',I3.3)
      write(filename, 990) itime/imodulo
call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,tb1,filename)

end subroutine VISU_PRE

!############################################################################
!
!subroutine BUDGET ()
!
!############################################################################

!USE param
!USE variables
!USE decomp_2d
!USE decomp_2d_io

!implicit none

!real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
!real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: uuv,uuvmean,vvv,vvvmean,vww,vwwmean 

!integer :: iperiod,ijk,nvect1,nvect2,nvect3,i,j,k

!iperiod=ilast-ifirst
!nvect1=xsize(1)*xsize(2)*xsize(3)

!Average in time
!uuv(:,:,:)=uuvmean(:,:,:)/iperiod
!vvv(:,:,:)=vvvmean(:,:,:)/iperiod
!vww(:,:,:)=vwwmean(:,:,:)/iperiod 







