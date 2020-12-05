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
     uvmean,uwmean,vwmean,phiphimean,tmean,utmean,vtmean,wtmean,dudx,dudy,dudz,dvdx,dvdy,dvdz,&
     dwdx,dwdy,dwdz,dudxdudx,dudydudy,dudzdudz,dvdxdvdx,dvdydvdy,dvdzdvdz,dwdxdwdx,&
     dwdydwdy,dwdzdwdz,uuvmean,vvvmean,vwwmean,pmean,pvmean,nxmsize,nymsize,nzmsize,&
     phG,ph2,ph3,pp3,dphidx,dphidxdphidx,dphidy,dphidydphidy,dphidz,dphidzdphidz,&
     dudxdphidx,dudydphidy,dudzdphidz,uphiumean,uphivmean,uphiwmean,phidpdx,dpdx,& 
     phidudx,phidudy,phidudz,udphidx,udphidy,udphidz,dudxdvdx,dudydvdy,dudzdvdz,&
     udpdy,dpdy,vdpdx,uvvmean)                  !Budget
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
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,&
                                                   uwmean,vwmean,tmean,utmean,vtmean,wtmean,dudx,dudy,&
                                                   dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dudxdudx,&
                                                   dudydudy,dudzdudz,dvdxdvdx,dvdydvdy,dvdzdvdz,&
                                                   dwdxdwdx,dwdydwdy,dwdzdwdz,uuvmean,vvvmean,vwwmean,&
                                                   pmean,pvmean,dphidx,dphidxdphidx,dphidy,dphidydphidy,&
                                                   dphidz,dphidzdphidz,dudxdphidx,dudydphidy,dudzdphidz,&
                                                   uphiumean,uphivmean,uphiwmean,phidpdx,dpdx,&
                                                   phidudx,phidudy,phidudz,udphidx,udphidy,udphidz,&
                                                   dudxdvdx,dudydvdy,dudzdvdz,udpdy,dpdy,vdpdx,uvvmean !Budget    
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: phimean, phiphimean                              !Budget
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,di1,td1,te1,tf1,tg1,th1               !Budget
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2, tb2,tc2,tf2,di2                  !Budget
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,td3,te3,di3                    !Budget

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

!dudx                                                                                !Budget
call derx (tb1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)           !Budget 
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dudx(:,:,:)=dudx(:,:,:)+tmean(:,:,:)                                                 !Budget

!phidudx=phi*dudx
tf1(:,:,:)=phi1(:,:,:)*tb1(:,:,:)
call fine_to_coarseS(1,tf1,tmean)
phidudx(:,:,:)=phidudx(:,:,:)+tmean(:,:,:)

!dudxdudx=dudx*dudx                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dudxdudx(:,:,:)=dudxdudx(:,:,:)+tmean(:,:,:)                                         !Budget


!dudy                                                                                !Budget
call transpose_x_to_y(ux1,ta2)                                                       !Budget
call dery (tb2,ta2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)           !Budget 
call transpose_y_to_x(tb2,tb1)                                                       !Budget
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dudy(:,:,:)=dudy(:,:,:)+tmean(:,:,:)                                                 !Budget

!phidudy=phi*dudy
tf1(:,:,:)=phi1(:,:,:)*tb1(:,:,:)
call fine_to_coarseS(1,tf1,tmean)
phidudy(:,:,:)=phidudy(:,:,:)+tmean(:,:,:)

!dudydudy=dudy*dudy                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dudydudy(:,:,:)=dudydudy(:,:,:)+tmean(:,:,:)                                         !Budget

!dudz                                                                                !Budget
call transpose_x_to_y(ux1,ta2)                                                       !Budget
call transpose_y_to_z(ta2,ta3)                                                       !Budget 
call derz (tb3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)               !Budget 
call transpose_z_to_y(tb3,tb2)                                                       !Budget
call transpose_y_to_x(tb2,tb1)                                                       !Budget
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dudz(:,:,:)=dudz(:,:,:)+tmean(:,:,:)                                                 !Budget

!phidudz=phi*dudz
tf1(:,:,:)=phi1(:,:,:)*tb1(:,:,:)
call fine_to_coarseS(1,tf1,tmean)
phidudz(:,:,:)=phidudz(:,:,:)+tmean(:,:,:)

!dudzdudz=dudz*dudz                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dudzdudz(:,:,:)=dudzdudz(:,:,:)+tmean(:,:,:)                                         !Budget

!dvdx                                                                                !Budget
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)           !Budget 
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dvdx(:,:,:)=dvdx(:,:,:)+tmean(:,:,:)                                                 !Budget

!dvdxdvdx=dvdx*dvdx                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dvdxdvdx(:,:,:)=dvdxdvdx(:,:,:)+tmean(:,:,:)                                         !Budget

!dvdy                                                                                !Budget
call transpose_x_to_y(uy1,ta2)                                                       !Budget
call dery (tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)           !Budget 
call transpose_y_to_x(tb2,tb1)                                                       !Budget
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dvdy(:,:,:)=dvdy(:,:,:)+tmean(:,:,:)                                                 !Budget

!dvdydvdy=dvdy*dvdy                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dvdydvdy(:,:,:)=dvdydvdy(:,:,:)+tmean(:,:,:)                                         !Budget

!dvdz                                                                                !Budget
call transpose_x_to_y(uy1,ta2)                                                       !Budget
call transpose_y_to_z(ta2,ta3)                                                       !Budget 
call derz (tb3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)               !Budget 
call transpose_z_to_y(tb3,tb2)                                                       !Budget
call transpose_y_to_x(tb2,tb1)                                                       !Budget
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dvdz(:,:,:)=dvdz(:,:,:)+tmean(:,:,:)                                                 !Budget

!dvdzdvdz=dvdz*dvdz                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dvdzdvdz(:,:,:)=dvdzdvdz(:,:,:)+tmean(:,:,:)                                         !Budget

!dwdx                                                                                !Budget
call derx (tb1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)           !Budget 
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dwdx(:,:,:)=dwdx(:,:,:)+tmean(:,:,:)                                                 !Budget

!dwdxdwdx=dwdx*dwdx                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dwdxdwdx(:,:,:)=dwdxdwdx(:,:,:)+tmean(:,:,:)                                         !Budget

!dwdy                                                                                !Budget
call transpose_x_to_y(uz1,ta2)                                                       !Budget
call dery (tb2,ta2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)           !Budget 
call transpose_y_to_x(tb2,tb1)                                                       !Budget
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dwdy(:,:,:)=dwdy(:,:,:)+tmean(:,:,:)                                                 !Budget

!dwdydwdy=dwdy*dwdy                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dwdydwdy(:,:,:)=dwdydwdy(:,:,:)+tmean(:,:,:)                                         !Budget

!dwdz                                                                                !Budget
call transpose_x_to_y(uz1,ta2)                                                       !Budget
call transpose_y_to_z(ta2,ta3)                                                       !Budget 
call derz (tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)               !Budget 
call transpose_z_to_y(tb3,tb2)                                                       !Budget
call transpose_y_to_x(tb2,tb1)                                                       !Budget
call fine_to_coarseS(1,tb1,tmean)                                                    !Budget
dwdz(:,:,:)=dwdz(:,:,:)+tmean(:,:,:)                                                 !Budget

!dwdzdwdz=dwdz*dwdz                                                                  !Budget  
tf1(:,:,:)=tb1(:,:,:)*tb1(:,:,:)                                                     !Budget
call fine_to_coarseS(1,tf1,tmean)                                                    !Budget 
dwdzdwdz(:,:,:)=dwdzdwdz(:,:,:)+tmean(:,:,:)                                         !Budget

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

!dpdx
call derx(ta1,tj1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call fine_to_coarseS(1,ta1,tmean)
dpdx(:,:,:)=dpdx(:,:,:)+tmean(:,:,:)

!phidpdx=phi*dpdx
tf1(:,:,:)=phi1(:,:,:)*ta1(:,:,:)
call fine_to_coarseS(1,tf1,tmean)
phidpdx(:,:,:)=phidpdx(:,:,:)+tmean(:,:,:)

!udpdy=u*dpdy
call transpose_x_to_y(tj1,ta2)
call dery(tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call transpose_y_to_x(tb2,td1)
te1(:,:,:)=ux1(:,:,:)*td1(:,:,:)
call fine_to_coarseS(1,te1,tmean)
udpdy(:,:,:)=udpdy(:,:,:)+tmean(:,:,:)

!dpdy
call fine_to_coarseS(1,td1,tmean)
dpdy(:,:,:)=dpdy(:,:,:)+tmean(:,:,:)

!vdpdx=v*dpdx
tf1(:,:,:)=uy1(:,:,:)*ta1(:,:,:)
call fine_to_coarseS(1,tf1,tmean)
vdpdx(:,:,:)=vdpdx(:,:,:)+tmean(:,:,:)

!pvmean=tj1*uy1                                                                      !Budget
ta1(:,:,:)=tj1(:,:,:)*uy1(:,:,:)                                                     !Budget                
call fine_to_coarseS(1,ta1,tmean)                                                    !Budget
pvmean(:,:,:)=pvmean(:,:,:)+tmean(:,:,:)                                             !Budget

if (iscalar==1) then
   !phimean=phi1
   call fine_to_coarseS(1,phi1,tmean)
   phimean(:,:,:)=phimean(:,:,:)+tmean(:,:,:)

   !dphidxdphidx=dphidx*dphidx
   call derx(ta1,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   tf1(:,:,:)=ta1(:,:,:)*ta1(:,:,:)
   call fine_to_coarseS(1,tf1,tmean)
   dphidxdphidx(:,:,:)=dphidxdphidx(:,:,:)+tmean(:,:,:)

   !dphidx
   call fine_to_coarseS(1,ta1,tmean)
   dphidx(:,:,:)=dphidx(:,:,:)+tmean(:,:,:)

   !udphidx=u*dphidx
   tf1(:,:,:)=ux1(:,:,:)*ta1(:,:,:)
   call fine_to_coarseS(1,tf1,tmean)
   udphidx(:,:,:)=udphidx(:,:,:)+tmean(:,:,:)

   !dphidydphidy=dphidy*dphidy
   call transpose_x_to_y(phi1,ta2)
   call dery(tb2,ta2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call transpose_y_to_x(tb2,td1)
   te1(:,:,:)=td1(:,:,:)*td1(:,:,:)
   call fine_to_coarseS(1,te1,tmean)
   dphidydphidy(:,:,:)=dphidydphidy(:,:,:)+tmean(:,:,:)

   !dphidy
   call fine_to_coarseS(1,td1,tmean)
   dphidy(:,:,:)=dphidy(:,:,:)+tmean(:,:,:)

   !udphidy=u*dphidy
   tf1(:,:,:)=ux1(:,:,:)*td1(:,:,:)
   call fine_to_coarseS(1,tf1,tmean)
   udphidy(:,:,:)=udphidy(:,:,:)+tmean(:,:,:)

   !dphidzdphidz=dphidz*dphidz
   call transpose_y_to_z(ta2,ta3)
   call derz(tb3,ta3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call transpose_z_to_y(tb3,tc2)
   call transpose_y_to_x(tc2,tg1)
   th1(:,:,:)=tg1(:,:,:)*tg1(:,:,:)
   call fine_to_coarseS(1,th1,tmean)
   dphidzdphidz(:,:,:)=dphidzdphidz(:,:,:)+tmean(:,:,:)

   !dphidz
   call fine_to_coarseS(1,tg1,tmean)
   dphidz(:,:,:)=dphidz(:,:,:)+tmean(:,:,:)

   !udphidz=u*dphidz
   tf1(:,:,:)=ux1(:,:,:)*tg1(:,:,:)
   call fine_to_coarseS(1,tf1,tmean)
   udphidz(:,:,:)=udphidz(:,:,:)+tmean(:,:,:)

   !dudxdphidx=dudx*dphidx
   call derx(tb1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx(ta1,phi1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0) 
   td1(:,:,:)=tb1(:,:,:)*ta1(:,:,:) 
   call fine_to_coarseS(1,td1,tmean)
   dudxdphidx(:,:,:)=dudxdphidx(:,:,:)+tmean(:,:,:) 

   !dudydphidy=dudy*dphidy
   call transpose_x_to_y(ux1,ta2)
   call transpose_x_to_y(phi1,tb2)
   call dery (tc2,ta2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1) 
   call dery (tf2,tb2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call transpose_y_to_x(tc2,ta1) 
   call transpose_y_to_x(tf2,tb1) 
   td1(:,:,:)=ta1(:,:,:)*tb1(:,:,:)
   call fine_to_coarseS(1,td1,tmean)
   dudydphidy(:,:,:)=dudydphidy(:,:,:)+tmean(:,:,:)

   !dudzdphidz=dudz*dphidz
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call derz (td3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) 
   call derz (te3,tb3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call transpose_z_to_y(td3,ta2)
   call transpose_z_to_y(te3,tb2)
   call transpose_y_to_x(ta2,ta1)
   call transpose_y_to_x(tb2,tb1)
   td1(:,:,:)=ta1(:,:,:)*tb1(:,:,:)
   call fine_to_coarseS(1,td1,tmean)
   dudzdphidz(:,:,:)=dudzdphidz(:,:,:)+tmean(:,:,:)

   !dudxdvdx=dudx*dvdx
   call derx(tb1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   call derx(ta1,uy1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
   td1(:,:,:)=tb1(:,:,:)*ta1(:,:,:)
   call fine_to_coarseS(1,td1,tmean)
   dudxdvdx(:,:,:)=dudxdvdx(:,:,:)+tmean(:,:,:)

   !dudydvdy=dudy*dvdy
   call transpose_x_to_y(ux1,ta2)
   call transpose_x_to_y(uy1,tb2)
   call dery (tc2,ta2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
   call dery (tf2,tb2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
   call transpose_y_to_x(tc2,ta1)
   call transpose_y_to_x(tf2,tb1)
   td1(:,:,:)=ta1(:,:,:)*tb1(:,:,:)
   call fine_to_coarseS(1,td1,tmean)
   dudydvdy(:,:,:)=dudydvdy(:,:,:)+tmean(:,:,:)

   !dudzdvdz=dudz*dvdz
   call transpose_y_to_z(ta2,ta3)
   call transpose_y_to_z(tb2,tb3)
   call derz (td3,ta3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1) 
   call derz (te3,tb3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
   call transpose_z_to_y(td3,ta2)
   call transpose_z_to_y(te3,tb2)
   call transpose_y_to_x(ta2,ta1)
   call transpose_y_to_x(tb2,tb1)
   td1(:,:,:)=ta1(:,:,:)*tb1(:,:,:)
   call fine_to_coarseS(1,td1,tmean)
   dudzdvdz(:,:,:)=dudzdvdz(:,:,:)+tmean(:,:,:)
 
   !uphiumean=ux1*phi1*ux1
   ta1(:,:,:)=ux1(:,:,:)*phi1(:,:,:)*ux1(:,:,:)      !Budget
   call fine_to_coarseS(1,ta1,tmean)                 !Budget
   uphiumean(:,:,:)=uphiumean(:,:,:)+tmean(:,:,:)    !Budget  

   !uphivmean=ux1*phi1*uy1
   ta1(:,:,:)=ux1(:,:,:)*phi1(:,:,:)*uy1(:,:,:)      !Budget
   call fine_to_coarseS(1,ta1,tmean)                 !Budget
   uphivmean(:,:,:)=uphivmean(:,:,:)+tmean(:,:,:)    !Budget

   !uphiwmean=ux1*phi1*uz1
   ta1(:,:,:)=ux1(:,:,:)*phi1(:,:,:)*uz1(:,:,:)      !Budget
   call fine_to_coarseS(1,ta1,tmean)                 !Budget
   uphiwmean(:,:,:)=uphiwmean(:,:,:)+tmean(:,:,:)    !Budget
     
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

!uvvmean=ux1*uy1*uy1
ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)*uy1(:,:,:)   !Budget
call fine_to_coarseS(1,ta1,tmean)             !Budget
uvvmean(:,:,:)=uvvmean(:,:,:)+tmean(:,:,:)    !Budget

!utmean=ux1*phi1
ta1(:,:,:)=ux1(:,:,:)*phi1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
utmean(:,:,:)=utmean(:,:,:)+tmean(:,:,:)

!vtmean=uy1*phi1
ta1(:,:,:)=uy1(:,:,:)*phi1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vtmean(:,:,:)=vtmean(:,:,:)+tmean(:,:,:)

!wtmean=uz1*phi1
ta1(:,:,:)=uz1(:,:,:)*phi1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
wtmean(:,:,:)=wtmean(:,:,:)+tmean(:,:,:)

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
   call decomp_2d_write_one(1,wtmean,'wtmean.dat',1)
   call decomp_2d_write_one(1,dudx,'dudx.dat',1)                   !Budget
   call decomp_2d_write_one(1,dudy,'dudy.dat',1)                   !Budget
   call decomp_2d_write_one(1,dudz,'dudz.dat',1)                   !Budget
   call decomp_2d_write_one(1,dvdx,'dvdx.dat',1)                   !Budget
   call decomp_2d_write_one(1,dvdy,'dvdy.dat',1)                   !Budget
   call decomp_2d_write_one(1,dvdz,'dvdz.dat',1)                   !Budget
   call decomp_2d_write_one(1,dwdx,'dwdx.dat',1)                   !Budget
   call decomp_2d_write_one(1,dwdy,'dwdy.dat',1)                   !Budget
   call decomp_2d_write_one(1,dwdz,'dwdz.dat',1)                   !Budget
   call decomp_2d_write_one(1,dudxdudx,'dudxdudx.dat',1)                   !Budget
   call decomp_2d_write_one(1,dudydudy,'dudydudy.dat',1)                   !Budget
   call decomp_2d_write_one(1,dudzdudz,'dudzdudz.dat',1)                   !Budget
   call decomp_2d_write_one(1,dvdxdvdx,'dvdxdvdx.dat',1)                   !Budget
   call decomp_2d_write_one(1,dvdydvdy,'dvdydvdy.dat',1)                   !Budget
   call decomp_2d_write_one(1,dvdzdvdz,'dvdzdvdz.dat',1)                   !Budget
   call decomp_2d_write_one(1,dwdxdwdx,'dwdxdwdx.dat',1)                   !Budget
   call decomp_2d_write_one(1,dwdydwdy,'dwdydwdy.dat',1)                   !Budget
   call decomp_2d_write_one(1,dwdzdwdz,'dwdzdwdz.dat',1)                   !Budget
   call decomp_2d_write_one(1,uuvmean,'uuvmean.dat',1)             !Budget
   call decomp_2d_write_one(1,vvvmean,'vvvmean.dat',1)             !Budget
   call decomp_2d_write_one(1,vwwmean,'vwwmean.dat',1)             !Budget
   call decomp_2d_write_one(1,pmean,'pmean.dat',1)               !Budget
   call decomp_2d_write_one(1,pvmean,'pvmean.dat',1)               !Budget
   call decomp_2d_write_one(1,dpdx,'dpdx.dat',1)               !Budget
   call decomp_2d_write_one(1,dudxdvdx,'dudxdvdx.dat',1)               !Budget
   call decomp_2d_write_one(1,dudydvdy,'dudydvdy.dat',1)               !Budget
   call decomp_2d_write_one(1,dudzdvdz,'dudzdvdz.dat',1)               !Budget
   call decomp_2d_write_one(1,udpdy,'udpdy.dat',1)               !Budget
   call decomp_2d_write_one(1,dpdy,'dpdy.dat',1)               !Budget
   call decomp_2d_write_one(1,vdpdx,'vdpdx.dat',1)               !Budget
   call decomp_2d_write_one(1,uvvmean,'uvvmean.dat',1)               !Budget

   if (nrank==0) print *,'write stat arrays velocity done!'
   if (iscalar==1) then
      call decomp_2d_write_one(1,phimean,'phimean.dat',1)
      call decomp_2d_write_one(1,phiphimean,'phiphimean.dat',1)
      call decomp_2d_write_one(1,dphidx,'dphidx.dat',1)
      call decomp_2d_write_one(1,dphidxdphidx,'dphidxdphidx.dat',1)
      call decomp_2d_write_one(1,dphidy,'dphidy.dat',1)
      call decomp_2d_write_one(1,dphidydphidy,'dphidydphidy.dat',1)
      call decomp_2d_write_one(1,dphidz,'dphidz.dat',1)
      call decomp_2d_write_one(1,dphidzdphidz,'dphidzdphidz.dat',1)
      call decomp_2d_write_one(1,dudxdphidx,'dudxdphidx.dat',1)
      call decomp_2d_write_one(1,dudydphidy,'dudydphidy.dat',1)
      call decomp_2d_write_one(1,dudzdphidz,'dudzdphidz.dat',1)
      call decomp_2d_write_one(1,uphiumean,'uphiumean.dat',1)
      call decomp_2d_write_one(1,uphivmean,'uphivmean.dat',1)
      call decomp_2d_write_one(1,uphiwmean,'uphiwmean.dat',1)
      call decomp_2d_write_one(1,phidpdx,'phidpdx.dat',1)
      call decomp_2d_write_one(1,phidudx,'phidudx.dat',1)
      call decomp_2d_write_one(1,phidudy,'phidudy.dat',1)
      call decomp_2d_write_one(1,phidudz,'phidudz.dat',1)
      call decomp_2d_write_one(1,udphidx,'udphidx.dat',1)
      call decomp_2d_write_one(1,udphidy,'udphidy.dat',1)
      call decomp_2d_write_one(1,udphidz,'udphidz.dat',1)
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







