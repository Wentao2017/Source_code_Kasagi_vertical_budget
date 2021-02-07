!***********************************************************
!
program post
!
!***********************************************************

implicit none

integer,parameter :: nzb=1
integer,parameter :: nx1=512,ny1=257,nz1=512
character(len=3) suffix,suffix1
character(len=20) nfichier,nfichier1,nchamp
integer :: num,i,j,k,nb_fich,nb_proc_total,nb_proc,longueur,long,itime
integer :: nxyz,count,itime1,itime2,itime3


real(8),dimension(nx1,ny1,nz1) :: umean1,vmean1,wmean1,phimean1,uumean1,vvmean1,wwmean1,phiphimean1
real(8),dimension(nx1,ny1,nz1) :: uvmean1,uwmean1,vwmean1,utmean1,vtmean1,wtmean1,dudx1,dudy1,dudz1,dvdx1, &
                                  dvdy1,dvdz1,dwdx1,dwdy1,dwdz1,dudxdudx1,dudydudy1,dudzdudz1, &
                                  dvdxdvdx1,dvdydvdy1,dvdzdvdz1,dwdxdwdx1,dwdydwdy1,dwdzdwdz1,uuvmean1, &
                                  vvvmean1,vwwmean1,pmean1,pvmean1,dphidx1,dphidxdphidx1,dphidy1, &
                                  dphidydphidy1,dphidz1,dphidzdphidz1,dudxdphidx1,dudydphidy1,dudzdphidz1, &
                                  uphiumean1,uphivmean1,uphiwmean1,phidpdx1,dpdx1,phidudx1,phidudy1,phidudz1,&
                                  udphidx1,udphidy1,udphidz1,uvvmean1,dudxdvdx1,dudydvdy1,dudzdvdz1,&
                                  udpdy1,dpdy1,vdpdx1,vttmean1,dvdxdphidx1,dvdydphidy1,dvdzdphidz1,phidpdy1,vphivmean1,&
                                  vphiwmean1,phidvdx1,phidvdy1,phidvdz1,vdphidx1,vdphidy1,vdphidz1
real(8),dimension(nx1,ny1,nz1) :: dudx2,dvdx2,dwdx2,dudy2,dvdy2,dwdy2,dudz2,dvdz2,dwdz2 
real(4),dimension(ny1) :: yp,ypi,ypii
real(8),dimension(nx1,ny1,nz1) :: umean2,vmean2,wmean2,uumean2,vvmean2,wwmean2
!real(8),dimension(nx1,ny1,nz1) :: uvmean2,uwmean2,vwmean2
real(4) :: u_to,u_to1,u_to2,re,xl2,xl3,xnu,x14,x15,pr,alpha,dxp,dzp,Re_tau_A,Re_tau_O
real(4) :: xlx=31.41593,zlz=12.56637,Gr=21344400
real(8),dimension(ny1,89) :: q_stat


open (15,file='yp.dat',form='formatted',status='unknown')
do j=1,ny1
   read(15,*) yp(j)
print *,yp(j)
enddo
close(15)
print *,'read yp'

itime1=100000


OPEN(11,FILE='umean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) umean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UX',umean1(nx1/2,ny1/2,k)/itime
    ! print *,k,'UX',umean1(nx1/2,ny1/2,k)/itime1
  ENDDO
  CLOSE(11)
OPEN(11,FILE='vmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UY',vmean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)
OPEN(11,FILE='wmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) wmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
    ! print *,k,'UZ',wmean1(nx1/2,ny1/2,k)/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uumean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uumean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UXUX',uumean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UYUY',vvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='wwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) wwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wwmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uuvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uuvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wwmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vvvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vvvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wwmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vwwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vwwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wwmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wwmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phimean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phimean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phiphimean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phiphimean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
  ENDDO
  CLOSE(11)

OPEN(11,FILE='utmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) utmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',utmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vtmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vtmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',vtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='wtmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) wtmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',vtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)


OPEN(11,FILE='dvdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dwdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dwdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dwdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dwdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dwdz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dwdz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudxdudx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudxdudx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudydudy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudydudy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudzdudz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudzdudz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdxdvdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdxdvdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdydvdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdydvdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdzdvdz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdzdvdz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dwdxdwdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dwdxdwdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dwdydwdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dwdydwdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dwdzdwdz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dwdzdwdz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='pmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) pmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='pvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) pvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',uvmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dphidx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dphidx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dphidxdphidx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dphidxdphidx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dphidy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dphidy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dphidydphidy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dphidydphidy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dphidz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dphidz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dphidzdphidz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dphidzdphidz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudxdphidx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudxdphidx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudydphidy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudydphidy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudzdphidz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudzdphidz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uphiumean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uphiumean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uphivmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uphivmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uphiwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uphiwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phidpdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidpdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dpdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dpdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11) 

OPEN(11,FILE='phidudx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidudx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phidudy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidudy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)
 
OPEN(11,FILE='phidudz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidudz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='udphidx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) udphidx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='udphidy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) udphidy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='udphidz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) udphidz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='uvvmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) uvvmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudxdvdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudxdvdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudydvdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudydvdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dudzdvdz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dudzdvdz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='udpdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) udpdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dpdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dpdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vdpdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vdpdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vttmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vttmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdxdphidx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdxdphidx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdydphidy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdydphidy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='dvdzdphidz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) dvdzdphidz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phidpdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidpdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vphivmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vphivmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vphiwmean.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vphiwmean1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phidvdx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidvdx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phidvdy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidvdy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='phidvdz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) phidvdz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vdphidx.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vdphidx1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vdphidy.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vdphidy1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

OPEN(11,FILE='vdphidz.dat',FORM='UNFORMATTED',&
       ACCESS='DIRECT', RECL=8, STATUS='OLD')
  COUNT = 1
  DO K=1,nz1
     DO J=1,ny1
        DO I=1,nx1
           READ(11,REC=COUNT) vdphidz1(I,J,K)
           COUNT = COUNT + 1
        ENDDO
     ENDDO
     !print *,k,'UZUZ',wtmean1(nx1/2,ny1/2,k)!/itime
  ENDDO
  CLOSE(11)

  print *,'READ DATA DONE 1'

do j=1,ny1-1
   ypi(j)=(yp(j)+yp(j+1))/2.
enddo


do j=1,ny1-2
   ypii(j)=(ypi(j)+ypi(j+1))/2.
enddo

umean1=umean1/itime1
vmean1=vmean1/itime1
wmean1=wmean1/itime1
uumean1=uumean1/itime1
vvmean1=vvmean1/itime1
wwmean1=wwmean1/itime1
uuvmean1=uuvmean1/itime1
vvvmean1=vvvmean1/itime1
vwwmean1=vwwmean1/itime1
vwmean1=vwmean1/itime1
phimean1=phimean1/itime1
phiphimean1=phiphimean1/itime1
utmean1=utmean1/itime1
vtmean1=vtmean1/itime1
wtmean1=wtmean1/itime1
uvmean1=uvmean1/itime1
uwmean1=uwmean1/itime1
dudx1=dudx1/itime1
dudy1=dudy1/itime1
dudz1=dudz1/itime1
dvdx1=dvdx1/itime1
dvdy1=dvdy1/itime1
dvdz1=dvdz1/itime1
dwdx1=dwdx1/itime1
dwdy1=dwdy1/itime1
dwdz1=dwdz1/itime1
dudxdudx1=dudxdudx1/itime1
dudydudy1=dudydudy1/itime1
dudzdudz1=dudzdudz1/itime1
dvdxdvdx1=dvdxdvdx1/itime1
dvdydvdy1=dvdydvdy1/itime1
dvdzdvdz1=dvdzdvdz1/itime1
dwdxdwdx1=dwdxdwdx1/itime1
dwdydwdy1=dwdydwdy1/itime1
dwdzdwdz1=dwdzdwdz1/itime1
pmean1=pmean1/itime1
pvmean1=pvmean1/itime1
dphidx1=dphidx1/itime1
dphidxdphidx1=dphidxdphidx1/itime1
dphidy1=dphidy1/itime1
dphidydphidy1=dphidydphidy1/itime1
dphidz1=dphidz1/itime1
dphidzdphidz1=dphidzdphidz1/itime1
dudxdphidx1=dudxdphidx1/itime1
dudydphidy1=dudydphidy1/itime1
dudzdphidz1=dudzdphidz1/itime1
uphiumean1=uphiumean1/itime1
uphivmean1=uphivmean1/itime1
uphiwmean1=uphiwmean1/itime1
phidpdx1=phidpdx1/itime1
dpdx1=dpdx1/itime1
phidudx1=phidudx1/itime1
phidudy1=phidudy1/itime1
phidudz1=phidudz1/itime1
udphidx1=udphidx1/itime1
udphidy1=udphidy1/itime1
udphidz1=udphidz1/itime1
uvvmean1=uvvmean1/itime1
dudxdvdx1=dudxdvdx1/itime1
dudydvdy1=dudydvdy1/itime1
dudzdvdz1=dudzdvdz1/itime1
udpdy1=udpdy1/itime1
dpdy1=dpdy1/itime1
vdpdx1=vdpdx1/itime1
vttmean1=vttmean1/itime1
dvdxdphidx1=dvdxdphidx1/itime1
dvdydphidy1=dvdydphidy1/itime1
dvdzdphidz1=dvdzdphidz1/itime1
phidpdy1=phidpdy1/itime1
vphivmean1=vphivmean1/itime1
vphiwmean1=vphiwmean1/itime1
phidvdx1=phidvdx1/itime1
phidvdy1=phidvdy1/itime1
phidvdz1=phidvdz1/itime1
vdphidx1=vdphidx1/itime1
vdphidy1=vdphidy1/itime1
vdphidz1=vdphidz1/itime1 

!DO AN AVERAGE IN X AND Z

q_stat=0.
do j=1,ny1
   do k=1,nz1
   do i=1,nx1
      q_stat(j,1)=q_stat(j,1)+umean1(i,j,k)
      q_stat(j,2)=q_stat(j,2)+vmean1(i,j,k)
      q_stat(j,3)=q_stat(j,3)+wmean1(i,j,k)
      q_stat(j,7)=q_stat(j,7)+phimean1(i,j,k)
!      q_stat(j,7)=q_stat(j,7)+phimean1(i,j,k)
      ! Direct Reynolds stresses
      q_stat(j,4)=q_stat(j,4)+uumean1(i,j,k)
      q_stat(j,5)=q_stat(j,5)+vvmean1(i,j,k)
      q_stat(j,6)=q_stat(j,6)+wwmean1(i,j,k)
      q_stat(j,8)=q_stat(j,8)+phiphimean1(i,j,k)
      q_stat(j,9)=q_stat(j,9)+utmean1(i,j,k)
      q_stat(j,10)=q_stat(j,10)+vtmean1(i,j,k)
      q_stat(j,11)=q_stat(j,11)+uvmean1(i,j,k)
      q_stat(j,12)=q_stat(j,12)+dudy1(i,j,k)
      q_stat(j,13)=q_stat(j,13)+uuvmean1(i,j,k)
      q_stat(j,14)=q_stat(j,14)+vvvmean1(i,j,k)
      q_stat(j,15)=q_stat(j,15)+vwwmean1(i,j,k)
      q_stat(j,16)=q_stat(j,16)+vwmean1(i,j,k)
      q_stat(j,18)=q_stat(j,18)+pmean1(i,j,k)
      q_stat(j,19)=q_stat(j,19)+pvmean1(i,j,k)
      q_stat(j,21)=q_stat(j,21)+dudx1(i,j,k)
      q_stat(j,22)=q_stat(j,22)+dudz1(i,j,k)
      q_stat(j,23)=q_stat(j,23)+dvdx1(i,j,k)
      q_stat(j,24)=q_stat(j,24)+dvdy1(i,j,k)
      q_stat(j,25)=q_stat(j,25)+dvdz1(i,j,k)
      q_stat(j,26)=q_stat(j,26)+dwdx1(i,j,k)
      q_stat(j,27)=q_stat(j,27)+dwdy1(i,j,k)
      q_stat(j,28)=q_stat(j,28)+dwdz1(i,j,k)
      q_stat(j,29)=q_stat(j,29)+dudxdudx1(i,j,k)
      q_stat(j,30)=q_stat(j,30)+dudydudy1(i,j,k)
      q_stat(j,31)=q_stat(j,31)+dudzdudz1(i,j,k)
      q_stat(j,32)=q_stat(j,32)+dvdxdvdx1(i,j,k)
      q_stat(j,33)=q_stat(j,33)+dvdydvdy1(i,j,k)
      q_stat(j,34)=q_stat(j,34)+dvdzdvdz1(i,j,k)
      q_stat(j,35)=q_stat(j,35)+dwdxdwdx1(i,j,k)
      q_stat(j,36)=q_stat(j,36)+dwdydwdy1(i,j,k)
      q_stat(j,37)=q_stat(j,37)+dwdzdwdz1(i,j,k)
      q_stat(j,39)=q_stat(j,39)+dphidx1(i,j,k)
      q_stat(j,40)=q_stat(j,40)+dphidxdphidx1(i,j,k)
      q_stat(j,41)=q_stat(j,41)+dphidy1(i,j,k)
      q_stat(j,42)=q_stat(j,42)+dphidydphidy1(i,j,k)
      q_stat(j,43)=q_stat(j,43)+dphidz1(i,j,k)
      q_stat(j,44)=q_stat(j,44)+dphidzdphidz1(i,j,k)
      q_stat(j,47)=q_stat(j,47)+dudxdphidx1(i,j,k)
      q_stat(j,48)=q_stat(j,48)+dudydphidy1(i,j,k)
      q_stat(j,49)=q_stat(j,49)+dudzdphidz1(i,j,k)
      q_stat(j,51)=q_stat(j,51)+uphiumean1(i,j,k)
      q_stat(j,52)=q_stat(j,52)+uphivmean1(i,j,k)
      q_stat(j,53)=q_stat(j,53)+uphiwmean1(i,j,k)
      q_stat(j,54)=q_stat(j,54)+wtmean1(i,j,k)
      q_stat(j,55)=q_stat(j,55)+uwmean1(i,j,k)
      q_stat(j,57)=q_stat(j,57)+phidpdx1(i,j,k)
      q_stat(j,58)=q_stat(j,58)+dpdx1(i,j,k)
      q_stat(j,60)=q_stat(j,60)+phidudx1(i,j,k)
      q_stat(j,61)=q_stat(j,61)+phidudy1(i,j,k)
      q_stat(j,62)=q_stat(j,62)+phidudz1(i,j,k)
      q_stat(j,63)=q_stat(j,63)+udphidx1(i,j,k)
      q_stat(j,64)=q_stat(j,64)+udphidy1(i,j,k)
      q_stat(j,65)=q_stat(j,65)+udphidz1(i,j,k)
      q_stat(j,68)=q_stat(j,68)+uvvmean1(i,j,k)
      q_stat(j,70)=q_stat(j,70)+dudxdvdx1(i,j,k)
      q_stat(j,71)=q_stat(j,71)+dudydvdy1(i,j,k)
      q_stat(j,72)=q_stat(j,72)+dudzdvdz1(i,j,k)
      q_stat(j,73)=q_stat(j,73)+udpdy1(i,j,k)
      q_stat(j,74)=q_stat(j,74)+dpdy1(i,j,k)
      q_stat(j,75)=q_stat(j,75)+vdpdx1(i,j,k)
      q_stat(j,77)=q_stat(j,77)+vttmean1(i,j,k)
      q_stat(j,78)=q_stat(j,78)+dvdxdphidx1(i,j,k)
      q_stat(j,79)=q_stat(j,79)+dvdydphidy1(i,j,k)
      q_stat(j,80)=q_stat(j,80)+dvdzdphidz1(i,j,k)
      q_stat(j,81)=q_stat(j,81)+phidpdy1(i,j,k)
      q_stat(j,82)=q_stat(j,82)+vphivmean1(i,j,k)
      q_stat(j,83)=q_stat(j,83)+vphiwmean1(i,j,k)
      q_stat(j,84)=q_stat(j,84)+phidvdx1(i,j,k)
      q_stat(j,85)=q_stat(j,85)+phidvdy1(i,j,k)
      q_stat(j,86)=q_stat(j,86)+phidvdz1(i,j,k)
      q_stat(j,87)=q_stat(j,87)+vdphidx1(i,j,k)
      q_stat(j,88)=q_stat(j,88)+vdphidy1(i,j,k)
      q_stat(j,89)=q_stat(j,89)+vdphidz1(i,j,k)
    

   enddo
   enddo
!   print *,j,umean1(nx1/2,j,nz1/2)
!   print 100,q_stat(j,1)/nx1/nz1,q_stat(j,4)/nx1/nz1,q_stat(j,4)/nx1/nz1-q_stat(j,1)/nx1/nz1*q_stat(j,1)/nx1/nz1
enddo
q_stat(:,:)=q_stat(:,:)/nx1/nz1
!   print *,q_stat(:,1)

xnu=1./3500!10322.629 !**********no s'e qu'e son estos num creo que la visco
re=3500!10322.629 !Creo qeu es el reynold usadopost_log_ut_budget.txt
pr=0.025
alpha=xnu/pr

!!!!!!!!!!!!
!Aiding flow
!!!!!!!!!!!!

   xl3=xnu*(q_stat(2,1)-q_stat(1,1))/(yp(2)-yp(1))
   xl2=sqrt(xl3)/xnu
   print *,'Re_tau_Aiding = ',xl2
   Re_tau_A=xl2

   if ((q_stat(ny1-1,1)-q_stat(ny1,1))>0) then
       Re_tau_O=sqrt(xnu*(q_stat(ny1-1,1)-q_stat(ny1,1))/(yp(ny1)-yp(ny1-1)))/xnu
   else
       Re_tau_O=-sqrt(-xnu*(q_stat(ny1-1,1)-q_stat(ny1,1))/(yp(ny1)-yp(ny1-1)))/xnu
   endif

   Re_tau_A=(Re_tau_A+Re_tau_O)/2.0  
   u_to1=xl2
   u_to=u_to1/re
   dxp=xlx/(nx1-1)*xl2
   dzp=zlz/(nz1-1)*xl2

!friction temperature
   x14=-(q_stat(2,7)-q_stat(1,7))/(yp(2)-yp(1))
   x15=x14/(pr*re*u_to) !friction temperature  
   print *,'Theta_tau_Aiding = ',x15

!do j=1,ny1
!   q_stat(j,56)=-(q_stat(j,51)-2*q_stat(j,1)*q_stat(j,9)-q_stat(j,7)*q_stat(j,4)+ &
!                 q_stat(j,52)-q_stat(j,1)*q_stat(j,10)-q_stat(j,7)*q_stat(j,11)-q_stat(j,2)*q_stat(j,9)+ &
!                 q_stat(j,53)-q_stat(j,1)*q_stat(j,54)-q_stat(j,7)*q_stat(j,55)-q_stat(j,3)*q_stat(j,9))/(u_to*x15*u_to)  !-utui 
!enddo

do j=1,ny1
    q_stat(j,56)=(q_stat(j,82)-2*q_stat(j,2)*q_stat(j,10)-q_stat(j,7)*q_stat(j,5))/(u_to*x15*u_to)    !vtv
enddo

do j=1,ny1-1
   q_stat(j,56)=-(q_stat(j+1,56)-q_stat(j,56))/((yp(j+1)-yp(j))*xl2)             !Turbulent diffusion  
enddo

do j=1,ny1
   q_stat(j,66)=-((q_stat(j,84)-q_stat(j,7)*q_stat(j,23)+q_stat(j,85)-q_stat(j,7)*q_stat(j,24)+q_stat(j,86)-q_stat(j,7)*q_stat(j,25))/(x15*u_to*xl2)+ &
                1/pr*(q_stat(j,87)-q_stat(j,2)*q_stat(j,39)+q_stat(j,88)-q_stat(j,2)*q_stat(j,41)+q_stat(j,89)-q_stat(j,2)*q_stat(j,43))/(u_to*x15*xl2))     !-phidudxi+(1/pr)*udphidxi
enddo

do j=1,ny1-1
   q_stat(j,66)=-(q_stat(j+1,66)-q_stat(j,66))/((yp(j+1)-yp(j))*xl2)             !Molecular diffusion  
enddo


do j=1,ny1
   q_stat(j,59)=-1000*(q_stat(j,81)-q_stat(j,7)*q_stat(j,74))/(u_to*u_to*xl2*x15)   !Temperature pressure-gradient correlation
enddo

do j=1,ny1
   print 100,q_stat(j,1),q_stat(j,4),q_stat(j,4)-q_stat(j,1)*q_stat(j,1),q_stat(j,8)-q_stat(j,7)*q_stat(j,7)
   q_stat(j,4)=sqrt((q_stat(j,4)-q_stat(j,1)*q_stat(j,1)))/u_to
   q_stat(j,5)=sqrt((q_stat(j,5)-q_stat(j,2)*q_stat(j,2)))/u_to
   q_stat(j,6)=sqrt((q_stat(j,6)-q_stat(j,3)*q_stat(j,3)))/u_to
   q_stat(j,8)=sqrt((q_stat(j,8)-q_stat(j,7)*q_stat(j,7)))/x15
   q_stat(j,9)=(q_stat(j,9)-q_stat(j,1)*q_stat(j,7))/(-u_to*x15)
   q_stat(j,10)=(q_stat(j,10)-q_stat(j,2)*q_stat(j,7))/(-u_to*x15)
   q_stat(j,11)=(q_stat(j,11)-q_stat(j,1)*q_stat(j,2))/(u_to*u_to)  
   q_stat(j,38)=-q_stat(j,5)**2*q_stat(j,41)*xnu/(u_to*x15)   !production
   q_stat(j,20)=(Gr/(2*xl2)**3)*q_stat(j,8)**2*x15        !Production by buoyancy
   q_stat(j,50)=-(1+1/Pr)*(q_stat(j,78)+q_stat(j,79)+q_stat(j,80)- &
                           q_stat(j,23)*q_stat(j,39)-q_stat(j,24)*q_stat(j,41)- &
                           q_stat(j,25)*q_stat(j,43))/(u_to*x15*xl2**2)           !Dissipation
   q_stat(j,45)=alpha*(q_stat(j,40)+q_stat(j,42)+q_stat(j,44)-q_stat(j,39)*q_stat(j,39)-q_stat(j,41)*&
                 q_stat(j,41)-q_stat(j,43)*q_stat(j,43))/(x15*x15)
   q_stat(j,46)=q_stat(j,11)*q_stat(j,41)/(q_stat(j,10)*q_stat(j,12))*u_to/x15

enddo

100 format(2F12.6,2F12.6,2F12.6,2F12.6,2F12.6,2F12.6,3x,2F12.6,3x,2F12.6)
!100 format(13F12.6)

  open (144,file='Aiding_flow_cf_Nu.dat',form='formatted',status='unknown')
  do j=1,ny1
      write(144,100) yp(j),yp(j)*xl2,((q_stat(j,1)/u_to)+(q_stat(j,1)/u_to))/2.,&
           (q_stat(j,4)+q_stat(j,4))/2.,(q_stat(j,5)+q_stat(j,5))/2.,&
           (q_stat(j,6)+q_stat(j,6))/2.,((q_stat(j,7)/x15)+(q_stat(j,7)/x15))/2.,&
           (q_stat(j,8)+q_stat(j,8))/2.,(q_stat(j,9)+q_stat(j,9))/2.,(q_stat(j,10)+q_stat(j,10))/2.,&
           (q_stat(j,11)+q_stat(j,11))/2.,(q_stat(j,45)+q_stat(j,45))/2.,(q_stat(j,46)+q_stat(j,46))/2.
   enddo
   close(144)

  open (144,file='Aiding_flow_budget_vt.dat',form='formatted',status='unknown')
  do j=1,ny1/2
      write(144,100) yp(j),yp(j)*xl2,q_stat(j,38),q_stat(j,20),q_stat(j,50),q_stat(j,59)
   enddo
   close(144)

  open (144,file='Aiding_flow_budget2_vt.dat',form='formatted',status='unknown')
  do j=1,ny1/2
      write(144,100) ypi(j),ypi(j)*xl2,q_stat(j,56),q_stat(j,66)
   enddo
   close(144)

  open (144,file='Aiding_flow_budget3_vt.dat',form='formatted',status='unknown')
  do j=1,ny1/2
      write(144,100) ypii(j),ypii(j)*xl2
   enddo
   close(144)

!!!!!!!!!!!!!!
!Opposing flow
!!!!!!!!!!!!!!

   xl3=xnu*(q_stat(ny1-1,1)-q_stat(ny1,1))/(yp(ny1)-yp(ny1-1))
   if (xl3>0) then
           xl2=sqrt(xl3)/xnu
           print *,'Re_tau_Opposing = ',xl2
   else
           xl2=-sqrt(-xl3)/xnu
           print *,'Re_tau_Opposing = ',xl2
   endif
   u_to1=xl2
   u_to=u_to1/re

!friction temperature
   x14=(q_stat(ny1-1,7)-q_stat(ny1,7))/(yp(ny1)-yp(ny1-1))
   x15=x14/(pr*re*u_to) !friction temperature  
   print *,'Theta_tau_Opposing = ',x15

q_stat=0.
do j=1,ny1
   do k=1,nz1
   do i=1,nx1
      q_stat(j,1)=q_stat(j,1)+umean1(i,j,k)
      q_stat(j,2)=q_stat(j,2)+vmean1(i,j,k)
      q_stat(j,3)=q_stat(j,3)+wmean1(i,j,k)
      q_stat(j,7)=q_stat(j,7)+phimean1(i,j,k)
!      q_stat(j,7)=q_stat(j,7)+phimean1(i,j,k)
      ! Direct Reynolds stresses
      q_stat(j,4)=q_stat(j,4)+uumean1(i,j,k)
      q_stat(j,5)=q_stat(j,5)+vvmean1(i,j,k)
      q_stat(j,6)=q_stat(j,6)+wwmean1(i,j,k)
      q_stat(j,8)=q_stat(j,8)+phiphimean1(i,j,k)
      q_stat(j,9)=q_stat(j,9)+utmean1(i,j,k)
      q_stat(j,10)=q_stat(j,10)+vtmean1(i,j,k)
      q_stat(j,11)=q_stat(j,11)+uvmean1(i,j,k)
      q_stat(j,12)=q_stat(j,12)+dudy1(i,j,k)
      q_stat(j,13)=q_stat(j,13)+uuvmean1(i,j,k)
      q_stat(j,14)=q_stat(j,14)+vvvmean1(i,j,k)
      q_stat(j,15)=q_stat(j,15)+vwwmean1(i,j,k)
      q_stat(j,16)=q_stat(j,16)+vwmean1(i,j,k)
      q_stat(j,18)=q_stat(j,18)+pmean1(i,j,k)
      q_stat(j,19)=q_stat(j,19)+pvmean1(i,j,k)
      q_stat(j,21)=q_stat(j,21)+dudx1(i,j,k)
      q_stat(j,22)=q_stat(j,22)+dudz1(i,j,k)
      q_stat(j,23)=q_stat(j,23)+dvdx1(i,j,k)
      q_stat(j,24)=q_stat(j,24)+dvdy1(i,j,k)
      q_stat(j,25)=q_stat(j,25)+dvdz1(i,j,k)
      q_stat(j,26)=q_stat(j,26)+dwdx1(i,j,k)
      q_stat(j,27)=q_stat(j,27)+dwdy1(i,j,k)
      q_stat(j,28)=q_stat(j,28)+dwdz1(i,j,k)
      q_stat(j,29)=q_stat(j,29)+dudxdudx1(i,j,k)
      q_stat(j,30)=q_stat(j,30)+dudydudy1(i,j,k)
      q_stat(j,31)=q_stat(j,31)+dudzdudz1(i,j,k)
      q_stat(j,32)=q_stat(j,32)+dvdxdvdx1(i,j,k)
      q_stat(j,33)=q_stat(j,33)+dvdydvdy1(i,j,k)
      q_stat(j,34)=q_stat(j,34)+dvdzdvdz1(i,j,k)
      q_stat(j,35)=q_stat(j,35)+dwdxdwdx1(i,j,k)
      q_stat(j,36)=q_stat(j,36)+dwdydwdy1(i,j,k)
      q_stat(j,37)=q_stat(j,37)+dwdzdwdz1(i,j,k)
      q_stat(j,39)=q_stat(j,39)+dphidx1(i,j,k)
      q_stat(j,40)=q_stat(j,40)+dphidxdphidx1(i,j,k)
      q_stat(j,41)=q_stat(j,41)+dphidy1(i,j,k)
      q_stat(j,42)=q_stat(j,42)+dphidydphidy1(i,j,k)
      q_stat(j,43)=q_stat(j,43)+dphidz1(i,j,k)
      q_stat(j,44)=q_stat(j,44)+dphidzdphidz1(i,j,k)
      q_stat(j,47)=q_stat(j,47)+dudxdphidx1(i,j,k)
      q_stat(j,48)=q_stat(j,48)+dudydphidy1(i,j,k)
      q_stat(j,49)=q_stat(j,49)+dudzdphidz1(i,j,k)
      q_stat(j,51)=q_stat(j,51)+uphiumean1(i,j,k)
      q_stat(j,52)=q_stat(j,52)+uphivmean1(i,j,k)
      q_stat(j,53)=q_stat(j,53)+uphiwmean1(i,j,k)
      q_stat(j,54)=q_stat(j,54)+wtmean1(i,j,k)
      q_stat(j,55)=q_stat(j,55)+uwmean1(i,j,k)
      q_stat(j,57)=q_stat(j,57)+phidpdx1(i,j,k)
      q_stat(j,58)=q_stat(j,58)+dpdx1(i,j,k)
      q_stat(j,60)=q_stat(j,60)+phidudx1(i,j,k)
      q_stat(j,61)=q_stat(j,61)+phidudy1(i,j,k)
      q_stat(j,62)=q_stat(j,62)+phidudz1(i,j,k)
      q_stat(j,63)=q_stat(j,63)+udphidx1(i,j,k)
      q_stat(j,64)=q_stat(j,64)+udphidy1(i,j,k)
      q_stat(j,65)=q_stat(j,65)+udphidz1(i,j,k)
      q_stat(j,68)=q_stat(j,68)+uvvmean1(i,j,k)
      q_stat(j,70)=q_stat(j,70)+dudxdvdx1(i,j,k)
      q_stat(j,71)=q_stat(j,71)+dudydvdy1(i,j,k)
      q_stat(j,72)=q_stat(j,72)+dudzdvdz1(i,j,k)
      q_stat(j,73)=q_stat(j,73)+udpdy1(i,j,k)
      q_stat(j,74)=q_stat(j,74)+dpdy1(i,j,k)
      q_stat(j,75)=q_stat(j,75)+vdpdx1(i,j,k)
      q_stat(j,77)=q_stat(j,77)+vttmean1(i,j,k)
      q_stat(j,78)=q_stat(j,78)+dvdxdphidx1(i,j,k)
      q_stat(j,79)=q_stat(j,79)+dvdydphidy1(i,j,k)
      q_stat(j,80)=q_stat(j,80)+dvdzdphidz1(i,j,k)
      q_stat(j,81)=q_stat(j,81)+phidpdy1(i,j,k)
      q_stat(j,82)=q_stat(j,82)+vphivmean1(i,j,k)
      q_stat(j,83)=q_stat(j,83)+vphiwmean1(i,j,k)
      q_stat(j,84)=q_stat(j,84)+phidvdx1(i,j,k)
      q_stat(j,85)=q_stat(j,85)+phidvdy1(i,j,k)
      q_stat(j,86)=q_stat(j,86)+phidvdz1(i,j,k)
      q_stat(j,87)=q_stat(j,87)+vdphidx1(i,j,k)
      q_stat(j,88)=q_stat(j,88)+vdphidy1(i,j,k)
      q_stat(j,89)=q_stat(j,89)+vdphidz1(i,j,k)

   enddo
   enddo
!   print *,j,umean1(nx1/2,j,nz1/2)
!   print 100,q_stat(j,1)/nx1/nz1,q_stat(j,4)/nx1/nz1,q_stat(j,4)/nx1/nz1-q_stat(j,1)/nx1/nz1*q_stat(j,1)/nx1/nz1
enddo
q_stat(:,:)=q_stat(:,:)/nx1/nz1

!do j=1,ny1
!   q_stat(j,56)=-(q_stat(j,51)-2*q_stat(j,1)*q_stat(j,9)-q_stat(j,7)*q_stat(j,4)+ &
!                 q_stat(j,52)-q_stat(j,1)*q_stat(j,10)-q_stat(j,7)*q_stat(j,11)-q_stat(j,2)*q_stat(j,9)+ &
!                 q_stat(j,53)-q_stat(j,1)*q_stat(j,54)-q_stat(j,7)*q_stat(j,55)-q_stat(j,3)*q_stat(j,9))/(u_to*x15*u_to)  !-utui 
!enddo

do j=1,ny1
    q_stat(j,56)=(q_stat(j,82)-2*q_stat(j,2)*q_stat(j,10)-q_stat(j,7)*q_stat(j,5))/(u_to*x15*u_to)    !vtv
enddo

do j=1,ny1-1
   q_stat(j,56)=-(q_stat(j+1,56)-q_stat(j,56))/((yp(j+1)-yp(j))*xl2)             !Turbulent diffusion  
enddo

do j=1,ny1
   q_stat(j,66)=-((q_stat(j,84)-q_stat(j,7)*q_stat(j,23)+q_stat(j,85)-q_stat(j,7)*q_stat(j,24)+q_stat(j,86)-q_stat(j,7)*q_stat(j,25))/(x15*u_to*xl2)+ &
                1/pr*(q_stat(j,87)-q_stat(j,2)*q_stat(j,39)+q_stat(j,88)-q_stat(j,2)*q_stat(j,41)+q_stat(j,89)-q_stat(j,2)*q_stat(j,43))/(u_to*x15*xl2))     !-phidudxi+(1/pr)*udphidxi
enddo

do j=1,ny1-1
   q_stat(j,66)=-(q_stat(j+1,66)-q_stat(j,66))/((yp(j+1)-yp(j))*xl2)             !Molecular diffusion  
enddo

do j=1,ny1
   q_stat(j,59)=-1000*(q_stat(j,81)-q_stat(j,7)*q_stat(j,74))/(u_to*u_to*xl2*x15)   !Temperature pressure-gradient correlation
enddo

do j=1,ny1
   print 100,q_stat(j,1),q_stat(j,4),q_stat(j,4)-q_stat(j,1)*q_stat(j,1),q_stat(j,8)-q_stat(j,7)*q_stat(j,7)
   q_stat(j,4)=sqrt((q_stat(j,4)-q_stat(j,1)*q_stat(j,1)))/u_to
   q_stat(j,5)=sqrt((q_stat(j,5)-q_stat(j,2)*q_stat(j,2)))/u_to
   q_stat(j,6)=sqrt((q_stat(j,6)-q_stat(j,3)*q_stat(j,3)))/u_to
   q_stat(j,8)=sqrt((q_stat(j,8)-q_stat(j,7)*q_stat(j,7)))/x15
   q_stat(j,9)=(q_stat(j,9)-q_stat(j,1)*q_stat(j,7))/(-u_to*x15)
   q_stat(j,10)=(q_stat(j,10)-q_stat(j,2)*q_stat(j,7))/(-u_to*x15)
   q_stat(j,11)=(q_stat(j,11)-q_stat(j,1)*q_stat(j,2))/(u_to*u_to)  
   q_stat(j,38)=-q_stat(j,5)**2*q_stat(j,41)*xnu/(u_to*x15)   !production
   q_stat(j,20)=(Gr/(2*xl2)**3)*q_stat(j,8)**2*x15       !Production by buoyancy
   q_stat(j,50)=-(1+1/Pr)*(q_stat(j,78)+q_stat(j,79)+q_stat(j,80)- &
                           q_stat(j,23)*q_stat(j,39)-q_stat(j,24)*q_stat(j,41)- &
                           q_stat(j,25)*q_stat(j,43))/(u_to*x15*xl2**2)           !Dissipation
   q_stat(j,45)=alpha*(q_stat(j,40)+q_stat(j,42)+q_stat(j,44)-q_stat(j,39)*q_stat(j,39)-q_stat(j,41)*&
                 q_stat(j,41)-q_stat(j,43)*q_stat(j,43))/(x15*x15)
   q_stat(j,46)=q_stat(j,11)*q_stat(j,41)/(q_stat(j,10)*q_stat(j,12))*u_to/x15
enddo

  open (144,file='Opposing_flow_cf_Nu.dat',form='formatted',status='unknown')
  do j=ny1,1,-1
      write(144,100) (2.0-yp(j)),(2.0-yp(j))*ABS(xl2),((q_stat(j,1)/u_to)+(q_stat(j,1)/u_to))/2.,&
           (q_stat(j,4)+q_stat(j,4))/2.,(q_stat(j,5)+q_stat(j,5))/2.,&
           (q_stat(j,6)+q_stat(j,6))/2.,((q_stat(j,7)/x15)+(q_stat(j,7)/x15))/2.,&
           (q_stat(j,8)+q_stat(j,8))/2.,(q_stat(j,9)+q_stat(j,9))/2.,(q_stat(j,10)+q_stat(j,10))/2.,&
           (q_stat(j,11)+q_stat(j,11))/2.,(q_stat(j,45)+q_stat(j,45))/2.,(q_stat(j,46)+q_stat(j,46))/2.
   enddo
   close(144)


  open (144,file='Opposing_flow_budget_vt.dat',form='formatted',status='unknown')
  do j=ny1,ny1/2+1,-1
      write(144,100) (2.0-yp(j)),(2.0-yp(j))*ABS(xl2),q_stat(j,38),q_stat(j,20),q_stat(j,50),q_stat(j,59)
   enddo
   close(144)
   
  open (144,file='Opposing_flow_budget2_vt.dat',form='formatted',status='unknown')
  do j=ny1-1,ny1/2+1,-1
      write(144,100) (2.0-ypi(j)),(2.0-ypi(j))*ABS(xl2),q_stat(j,56),q_stat(j,66)
   enddo
   close(144)

  open (144,file='Opposing_flow_budget3_vt.dat',form='formatted',status='unknown')
  do j=ny1-2,ny1/2+1,-1
      write(144,100) (2.0-ypii(j)),(2.0-ypii(j))*ABS(xl2)
   enddo
   close(144)

end program post
!******************************************************************
!
subroutine numcar (num,car)
!
!******************************************************************!

character(len=3) car

if (num.ge.100) then
   write (car,1) num
1  format (i3)
else
   if (num.ge.10) then
      write (car,2) num
2     format ('0',i2)
   else
      write (car,3) num
3     format ('00',i1)
   endif
endif

return
end subroutine numcar

!1. y
!2. y+
!3. u+
!4. u_rms
!5. v_rms
!6. w_rms
!7. t+
!8. t_rms
!9. ut
!10. vt
!11. uv
!12. wt
!13. Epsilon_theta_wrong
!14. Epsilon_theta_correct
!15. Prt
