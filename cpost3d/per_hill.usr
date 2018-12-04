C-----------------------------------------------------------------------
C  nek5000 user-file template
C
C  user specified routines:
C     - userbc : boundary conditions
C     - useric : initial conditions
C     - uservp : variable properties
C     - userf  : local acceleration term for fluid
C     - userq  : local source term for scalars
C     - userchk: general purpose routine for checking errors etc. 
C
C-----------------------------------------------------------------------

      subroutine uservp(ix,iy,iz,eg) ! set variable properties
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      udiff  = 0.0
      utrans = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userf(ix,iy,iz,ieg) 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      ffx = 0.0 
      ffy = 0.0
      ffz = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userq(ix,iy,iz,eg) ! set source term
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      integer e,f,eg
c     e = gllel(eg)

      qvol   = 0.0
      source = 0.0

      return
      end

c-----------------------------------------------------------------------

      subroutine userbc(ix,iy,iz,iside,ieg) 
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      ux = 0.0
      uy = 0.0
      uz = 0.0 
      
      return
      end

c-----------------------------------------------------------------------

      subroutine useric(ix,iy,iz,ieg) ! set up initial conditions
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      ux = 0.0
      uy = 0.0
      uz = 0.0
      
      return
      end

c-----------------------------------------------------------------------

      subroutine userchk()

      implicit none

! !      include 'SIZE_DEF'      
      include 'SIZE'
! !      include 'GEOM_DEF'
      include 'GEOM'                    ! xm1, ym1, zm1
! !      include 'SOLN_DEF'
      include 'SOLN'                    ! T
! !      include 'MASS_DEF'
      include 'MASS'                    !BM1 for lambda2
! !      include 'TSTEP_DEF'
      include 'TSTEP'                   ! ISTEP
! !      include 'INPUT_DEF'
      include 'INPUT'                   ! PARAM(12) (DT)
      include 'USERPAR'                 ! l2freq, FIXGEOM, NEW_DT, STATS3D
      include 'STATS'
      
      character*132 inputname1,hdr,pippo,pippa
      integer nfiles,ssf,sfn,ntot,ifld,j
      real nu,rho,mtimes,ftime,ltime,fdt,fttime,alpha
      real dum1(lx1*ly1*lz1*lelt),dum2(lx1*ly1*lz1*lelt)
      real U(lx1*ly1*lz1*lelt),V(lx1*ly1*lz1*lelt),W(lx1*ly1*lz1*lelt)
      real P(lx1*ly1*lz1*lelt),uu(lx1*ly1*lz1*lelt),vv(lx1*ly1*lz1*lelt)
      real ww(lx1*ly1*lz1*lelt),uv(lx1*ly1*lz1*lelt)
      real uw(lx1*ly1*lz1*lelt),vw(lx1*ly1*lz1*lelt)
      real pp(lx1*ly1*lz1*lelt),ppp(lx1*ly1*lz1*lelt)
      real pppp(lx1*ly1*lz1*lelt),uuu(lx1*ly1*lz1*lelt)
      real vvv(lx1*ly1*lz1*lelt),www(lx1*ly1*lz1*lelt)
      real uuv(lx1*ly1*lz1*lelt),uuw(lx1*ly1*lz1*lelt)
      real uvv(lx1*ly1*lz1*lelt),vvw(lx1*ly1*lz1*lelt)
      real uww(lx1*ly1*lz1*lelt),vww(lx1*ly1*lz1*lelt)
      real uvw(lx1*ly1*lz1*lelt),dUdx(lx1*ly1*lz1*lelt)
      real dUdy(lx1*ly1*lz1*lelt),dUdz(lx1*ly1*lz1*lelt)
      real dVdx(lx1*ly1*lz1*lelt),dVdy(lx1*ly1*lz1*lelt)
      real dVdz(lx1*ly1*lz1*lelt),dWdx(lx1*ly1*lz1*lelt)
      real dWdy(lx1*ly1*lz1*lelt),dWdz(lx1*ly1*lz1*lelt)
      real Pxx(lx1*ly1*lz1*lelt),Pyy(lx1*ly1*lz1*lelt)
      real Pzz(lx1*ly1*lz1*lelt),Pxy(lx1*ly1*lz1*lelt)
      real Pyz(lx1*ly1*lz1*lelt),Pxz(lx1*ly1*lz1*lelt)
      real Dxx(lx1*ly1*lz1*lelt),Dyy(lx1*ly1*lz1*lelt)
      real Dzz(lx1*ly1*lz1*lelt),Dxy(lx1*ly1*lz1*lelt)
      real Dyz(lx1*ly1*lz1*lelt),Dxz(lx1*ly1*lz1*lelt)
      real duudx(lx1*ly1*lz1*lelt),duudy(lx1*ly1*lz1*lelt)
      real duudz(lx1*ly1*lz1*lelt),dvvdx(lx1*ly1*lz1*lelt)
      real dvvdy(lx1*ly1*lz1*lelt),dvvdz(lx1*ly1*lz1*lelt)
      real dwwdx(lx1*ly1*lz1*lelt),dwwdy(lx1*ly1*lz1*lelt)
      real dwwdz(lx1*ly1*lz1*lelt),duvdx(lx1*ly1*lz1*lelt)
      real duvdy(lx1*ly1*lz1*lelt),duvdz(lx1*ly1*lz1*lelt)
      real duwdx(lx1*ly1*lz1*lelt),duwdy(lx1*ly1*lz1*lelt)
      real duwdz(lx1*ly1*lz1*lelt),dvwdx(lx1*ly1*lz1*lelt)
      real dvwdy(lx1*ly1*lz1*lelt),dvwdz(lx1*ly1*lz1*lelt)
      real Cxx(lx1*ly1*lz1*lelt),Cyy(lx1*ly1*lz1*lelt)
      real Czz(lx1*ly1*lz1*lelt),Cxy(lx1*ly1*lz1*lelt)
      real Cyz(lx1*ly1*lz1*lelt),Cxz(lx1*ly1*lz1*lelt)
      real d2uudx2(lx1*ly1*lz1*lelt),d2uudy2(lx1*ly1*lz1*lelt)
      real d2uudz2(lx1*ly1*lz1*lelt),d2vvdx2(lx1*ly1*lz1*lelt)
      real d2vvdy2(lx1*ly1*lz1*lelt),d2vvdz2(lx1*ly1*lz1*lelt)
      real d2wwdx2(lx1*ly1*lz1*lelt),d2wwdy2(lx1*ly1*lz1*lelt)
      real d2wwdz2(lx1*ly1*lz1*lelt),d2uvdx2(lx1*ly1*lz1*lelt)
      real d2uvdy2(lx1*ly1*lz1*lelt),d2uvdz2(lx1*ly1*lz1*lelt)
      real d2uwdx2(lx1*ly1*lz1*lelt),d2uwdy2(lx1*ly1*lz1*lelt)
      real d2uwdz2(lx1*ly1*lz1*lelt),d2vwdx2(lx1*ly1*lz1*lelt)
      real d2vwdy2(lx1*ly1*lz1*lelt),d2vwdz2(lx1*ly1*lz1*lelt)
      real VDxx(lx1*ly1*lz1*lelt),VDyy(lx1*ly1*lz1*lelt)
      real VDzz(lx1*ly1*lz1*lelt),VDxy(lx1*ly1*lz1*lelt)
      real VDyz(lx1*ly1*lz1*lelt),VDxz(lx1*ly1*lz1*lelt)
      real duuudx(lx1*ly1*lz1*lelt),duvvdx(lx1*ly1*lz1*lelt)
      real duvvdy(lx1*ly1*lz1*lelt),duwwdx(lx1*ly1*lz1*lelt)
      real duwwdz(lx1*ly1*lz1*lelt),duuvdx(lx1*ly1*lz1*lelt)
      real duuvdy(lx1*ly1*lz1*lelt),duuwdx(lx1*ly1*lz1*lelt)
      real duuwdz(lx1*ly1*lz1*lelt),duvwdx(lx1*ly1*lz1*lelt)
      real duvwdy(lx1*ly1*lz1*lelt),duvwdz(lx1*ly1*lz1*lelt)
      real dvvvdy(lx1*ly1*lz1*lelt),dvwwdy(lx1*ly1*lz1*lelt)
      real dvwwdz(lx1*ly1*lz1*lelt),dvvwdy(lx1*ly1*lz1*lelt)
      real dvvwdz(lx1*ly1*lz1*lelt),dwwwdz(lx1*ly1*lz1*lelt)
      real Txx(lx1*ly1*lz1*lelt),Tyy(lx1*ly1*lz1*lelt)
      real Tzz(lx1*ly1*lz1*lelt),Txy(lx1*ly1*lz1*lelt)
      real Tyz(lx1*ly1*lz1*lelt),Txz(lx1*ly1*lz1*lelt)
      real PTxx(lx1*ly1*lz1*lelt),PTyy(lx1*ly1*lz1*lelt)
      real PTzz(lx1*ly1*lz1*lelt),PTxy(lx1*ly1*lz1*lelt)
      real PTyz(lx1*ly1*lz1*lelt),PTxz(lx1*ly1*lz1*lelt)
      real dpudx(lx1*ly1*lz1*lelt),dpudy(lx1*ly1*lz1*lelt)
      real dpudz(lx1*ly1*lz1*lelt),dpvdx(lx1*ly1*lz1*lelt)
      real dpvdy(lx1*ly1*lz1*lelt),dpvdz(lx1*ly1*lz1*lelt)
      real dpwdx(lx1*ly1*lz1*lelt),dpwdy(lx1*ly1*lz1*lelt)
      real dpwdz(lx1*ly1*lz1*lelt),dPdx(lx1*ly1*lz1*lelt)
      real dPdy(lx1*ly1*lz1*lelt),dPdz(lx1*ly1*lz1*lelt)
      real PSxx(lx1*ly1*lz1*lelt),PSyy(lx1*ly1*lz1*lelt)
      real PSzz(lx1*ly1*lz1*lelt),PSxy(lx1*ly1*lz1*lelt)
      real PSyz(lx1*ly1*lz1*lelt),PSxz(lx1*ly1*lz1*lelt)
      real pdudx(lx1*ly1*lz1*lelt),pdudy(lx1*ly1*lz1*lelt)
      real pdudz(lx1*ly1*lz1*lelt),pdvdx(lx1*ly1*lz1*lelt)
      real pdvdy(lx1*ly1*lz1*lelt),pdvdz(lx1*ly1*lz1*lelt)
      real pdwdx(lx1*ly1*lz1*lelt),pdwdy(lx1*ly1*lz1*lelt)
      real pdwdz(lx1*ly1*lz1*lelt)
      real Pixx(lx1*ly1*lz1*lelt),Piyy(lx1*ly1*lz1*lelt)
      real Pizz(lx1*ly1*lz1*lelt),Pixy(lx1*ly1*lz1*lelt)
      real Piyz(lx1*ly1*lz1*lelt),Pixz(lx1*ly1*lz1*lelt)
      real Pk(lx1*ly1*lz1*lelt),Dk(lx1*ly1*lz1*lelt)
      real Tk(lx1*ly1*lz1*lelt),VDk(lx1*ly1*lz1*lelt)
      real Pik(lx1*ly1*lz1*lelt),Ck(lx1*ly1*lz1*lelt)
      real Resk(lx1*ly1*lz1*lelt)
      real STATS_TEMP(lx1*ly1*lz1*lelt,STAT_LVAR)

      nfiles=param(68)          ! Number of stat tiles
      nu=param(2)               ! Kinematic viscosity
      rho=1.0                   ! Fluid density
      mtimes=0.0                ! Starting time CHANGE
      
      ntot=lx1*ly1*lz1*lelt
      alpha=1.0
      
      call rzero(STAT,lx1*ly1*lz1*lelt*STAT_LVAR)
      call rzero(STATS_TEMP,lx1*ly1*lz1*lelt*STAT_LVAR)

      do sfn=1,11
         write(pippa,'(i2.2)') sfn
         ltime=mtimes           ! Time of last field
         ftime=0.0              ! Time of current field
         fttime=0.0             ! Total accumulated time
         do ssf = 1,nfiles
         
            write(pippo,'(i5.5)') ssf
            inputname1 = 'STAT3D/s'//trim(pippa)//'per_hill0.f'//
     &           trim(pippo)
            call read_hdr(inputname1,ftime) ! We read header to get the times

            fdt=ftime-ltime     ! Time span of this field
            fttime=fttime+fdt   ! Total averaging time
         
            if (nid.eq.0) write(*,*) '**FIELD,Ts,Tf,Ta',ssf,ltime,ftime,
     &           fdt
         
            ltime=ftime         ! Update last field time
         
            call load_field(inputname1)

            call add2sxy(STATS_TEMP(1,4*(sfn-1)+1),alpha,vx,fdt,ntot)
            call add2sxy(STATS_TEMP(1,4*(sfn-1)+2),alpha,vy,fdt,ntot)
            call add2sxy(STATS_TEMP(1,4*(sfn-1)+3),alpha,vz,fdt,ntot)
            call add2sxy(STATS_TEMP(1,4*(sfn-1)+4),alpha,t,fdt,ntot)

         enddo
      enddo

!     Reorder fields. No change for the first 26 variables
      do ifld=1,26
         do j=1,ntot
            STAT(j,ifld) = STATS_TEMP(j,ifld)
         enddo
      enddo

!     Various shifts in position
      do ifld=27,32
         do j=1,ntot
            STAT(j,ifld+1) = STATS_TEMP(j,ifld)
         enddo
      enddo
      
      do j=1,ntot
         STAT(j,27) = STATS_TEMP(j,33)
      enddo

      do j=1,ntot
         STAT(j,38) = STATS_TEMP(j,34)
      enddo
      
      do ifld=35,38
         do j=1,ntot
            STAT(j,ifld-1) = STATS_TEMP(j,ifld)
         enddo
      enddo
      
!     No change for the last variables
      do ifld=39,44
         do j=1,ntot
            STAT(j,ifld) = STATS_TEMP(j,ifld)
         enddo
      enddo
      
!     Divide by the total averaging time
      do ifld=1,STAT_LVAR
         do j=1,ntot
            STAT(j,ifld) = STAT(j,ifld)/fttime
         enddo
      enddo

!     Mean velocities. Tensors of Rank 1.
      do j=1,ntot
         U(j)=STAT(j,1)
         V(j)=STAT(j,2)
         W(j)=STAT(j,3)
      enddo

!     Reynolds-stress tensor. Tensor of Rank 2.
      do j=1,ntot
         uu(j)=STAT(j,5)-U(j)*U(j)
         vv(j)=STAT(j,6)-V(j)*V(j)
         ww(j)=STAT(j,7)-W(j)*W(j)
         uv(j)=STAT(j,9)-U(j)*V(j)
         uw(j)=STAT(j,11)-U(j)*W(j)
         vw(j)=STAT(j,10)-V(j)*W(j)
      enddo

!     Mean, RMS, skewness and flatness of pressure
      do j=1,ntot
         P(j)=STAT(j,4)
         pp(j)=STAT(j,8)-P(j)*P(j)
         ppp(j)=STAT(j,27)-3*P(j)*pp(j)-P(j)*P(j)*P(j)
         pppp(j)=STAT(j,38)-4*P(j)*ppp(j)-6*P(j)*P(j)*pp(j)-
     &        P(j)*P(j)*P(j)*P(j)
      enddo

!     Skewness tensor. Tensor of Rank 3.
!     Form of the tensor.
!     [ uuu uuv uuw ] [ uuv uvv uvw ] [ uuw uvw uww ]
!     [ uuv uvv uvw ] [ uvv vvv vvw ] [ uvw vvw vww ]
!     [ uuw uvw uww ] [ uvw vvw vww ] [ uww vww www ]
      do j=1,ntot
         uuu(j)=STAT(j,24)-3*U(j)*uu(j)-U(j)*U(j)*U(j)
         vvv(j)=STAT(j,25)-3*V(j)*vv(j)-V(j)*V(j)*V(j)
         www(j)=STAT(j,26)-3*W(j)*ww(j)-W(j)*W(j)*W(j)
         uuv(j)=STAT(j,28)-2*U(j)*uv(j)-V(j)*uu(j)-U(j)*U(j)*V(j)
         uuw(j)=STAT(j,29)-2*U(j)*uw(j)-W(j)*uu(j)-U(j)*U(j)*W(j)
         uvv(j)=STAT(j,30)-2*V(j)*uv(j)-U(j)*vv(j)-V(j)*V(j)*U(j)
         vvw(j)=STAT(j,31)-2*V(j)*vw(j)-W(j)*vv(j)-V(j)*V(j)*W(j)
         uww(j)=STAT(j,32)-2*W(j)*uw(j)-U(j)*ww(j)-W(j)*W(j)*U(j)
         vww(j)=STAT(j,33)-2*W(j)*vw(j)-V(j)*ww(j)-W(j)*W(j)*V(j)
         uvw(j)=STAT(j,34)-U(j)*vw(j)-V(j)*uw(j)-W(j)*uv(j)-
     &        U(j)*V(j)*W(j)
      enddo

!     Velocity gradient tensor. Tensor of Rank 2.
      call gradm1(dUdx,dUdy,dUdz,U)
      call gradm1(dVdx,dVdy,dVdz,V)
      call gradm1(dWdx,dWdy,dWdz,W)

!     Production tensor. Tensor of Rank 2.
      do j=1,ntot
         Pxx(j)=-2*(uu(j)*dUdx(j)+uv(j)*dUdy(j)+uw(j)*dUdz(j))
         Pyy(j)=-2*(uv(j)*dVdx(j)+vv(j)*dVdy(j)+vw(j)*dVdz(j))
         Pzz(j)=-2*(uw(j)*dWdx(j)+vw(j)*dWdy(j)+ww(j)*dWdz(j))
         Pxy(j)=-(uu(j)*dVdx(j)+uv(j)*dVdy(j)+uv(j)*dUdx(j)+
     &        vv(j)*dUdy(j)+uw(j)*dVdz(j)+vw(j)*dUdz(j))
         Pxz(j)=-(uu(j)*dWdx(j)+uv(j)*dWdy(j)+uw(j)*dUdx(j)+
     &        vw(j)*dUdy(j)+uw(j)*dWdz(j)+ww(j)*dUdz(j))
         Pyz(j)=-(uv(j)*dWdx(j)+vv(j)*dWdy(j)+uw(j)*dVdx(j)+
     &        vw(j)*dVdy(j)+vw(j)*dWdz(j)+ww(j)*dVdz(j))
      enddo

!     Dissipation tensor. Tensor of Rank 2.
      do j=1,ntot
         Dxx(j)=-2*nu*(STAT(j,39)-dUdx(j)*dUdx(j)-dUdy(j)*dUdy(j)-
     &        dUdz(j)*dUdz(j))
         Dyy(j)=-2*nu*(STAT(j,40)-dVdx(j)*dVdx(j)-dVdy(j)*dVdy(j)-
     &        dVdz(j)*dVdz(j))
         Dzz(j)=-2*nu*(STAT(j,41)-dWdx(j)*dWdx(j)-dWdy(j)*dWdy(j)-
     &        dWdz(j)*dWdz(j))
         Dxy(j)=-2*nu*(STAT(j,42)-dUdx(j)*dVdx(j)-dUdy(j)*dVdy(j)-
     &        dUdz(j)*dVdz(j))
         Dxz(j)=-2*nu*(STAT(j,43)-dUdx(j)*dWdx(j)-dUdy(j)*dWdy(j)-
     &        dUdz(j)*dWdz(j))
         Dyz(j)=-2*nu*(STAT(j,44)-dVdx(j)*dWdx(j)-dVdy(j)*dWdy(j)-
     &        dVdz(j)*dWdz(j))
      enddo

!     Derivatives of the Reynolds-stress tensor components
      call gradm1(duudx,duudy,duudz,uu)
      call gradm1(dvvdx,dvvdy,dvvdz,vv)
      call gradm1(dwwdx,dwwdy,dwwdz,ww)
      call gradm1(duvdx,duvdy,duvdz,uv)
      call gradm1(duwdx,duwdy,duwdz,uw)
      call gradm1(dvwdx,dvwdy,dvwdz,vw)

!     Mean convection tensor. Tensor of Rank 2.
      do j=1,ntot
         Cxx(j)=U(j)*duudx(j)+V(j)*duudy(j)+W(j)*duudz(j)
         Cyy(j)=U(j)*dvvdx(j)+V(j)*dvvdy(j)+W(j)*dvvdz(j)
         Czz(j)=U(j)*dwwdx(j)+V(j)*dwwdy(j)+W(j)*dwwdz(j)
         Cxy(j)=U(j)*duvdx(j)+V(j)*duvdy(j)+W(j)*duvdz(j)
         Cxz(j)=U(j)*duwdx(j)+V(j)*duwdy(j)+W(j)*duwdz(j)
         Cyz(j)=U(j)*dvwdx(j)+V(j)*dvwdy(j)+W(j)*dvwdz(j)
      enddo
      
!     Second derivatives of the Reynolds-stress tensor components
      call gradm1(d2uudx2,dum1,dum2,duudx)
      call gradm1(dum1,d2uudy2,dum2,duudy)
      call gradm1(dum1,dum2,d2uudz2,duudz)
      call gradm1(d2vvdx2,dum1,dum2,dvvdx)
      call gradm1(dum1,d2vvdy2,dum2,dvvdy)
      call gradm1(dum1,dum2,d2vvdz2,dvvdz)
      call gradm1(d2wwdx2,dum1,dum2,dwwdx)
      call gradm1(dum1,d2wwdy2,dum2,dwwdy)
      call gradm1(dum1,dum2,d2wwdz2,dwwdz)
      call gradm1(d2uvdx2,dum1,dum2,duvdx)
      call gradm1(dum1,d2uvdy2,dum2,duvdy)
      call gradm1(dum1,dum2,d2uvdz2,duvdz)
      call gradm1(d2uwdx2,dum1,dum2,duwdx)
      call gradm1(dum1,d2uwdy2,dum2,duwdy)
      call gradm1(dum1,dum2,d2uwdz2,duwdz)
      call gradm1(d2vwdx2,dum1,dum2,dvwdx)
      call gradm1(dum1,d2vwdy2,dum2,dvwdy)
      call gradm1(dum1,dum2,d2vwdz2,dvwdz)
      
!     Viscous diffusion tensor. Tensor of Rank 2.
      do j=1,ntot
         VDxx(j)=nu*(d2uudx2(j)+d2uudy2(j)+d2uudz2(j))
         VDyy(j)=nu*(d2vvdx2(j)+d2vvdy2(j)+d2vvdz2(j))
         VDzz(j)=nu*(d2wwdx2(j)+d2wwdy2(j)+d2wwdz2(j))
         VDxy(j)=nu*(d2uvdx2(j)+d2uvdy2(j)+d2uvdz2(j))
         VDxz(j)=nu*(d2uwdx2(j)+d2uwdy2(j)+d2uwdz2(j))
         VDyz(j)=nu*(d2vwdx2(j)+d2vwdy2(j)+d2vwdz2(j))
      enddo
      
!     Derivatives of the triple-product terms
      call gradm1(duuudx,dum1,dum2,uuu)
      call gradm1(duvvdx,duvvdy,dum2,uvv)
      call gradm1(duwwdx,dum1,duwwdz,uww)
      call gradm1(duuvdx,duuvdy,dum2,uuv)
      call gradm1(duuwdx,dum1,duuwdz,uuw)
      call gradm1(duvwdx,duvwdy,duvwdz,uvw)
      call gradm1(dum1,dvvvdy,dum2,vvv)
      call gradm1(dum1,dvwwdy,dvwwdz,vww)
      call gradm1(dum1,dvvwdy,dvvwdz,vvw)
      call gradm1(dum1,dum2,dwwwdz,www)

!     Turbulent transport tensor. Tensor of Rank 2.
      do j=1,ntot
         Txx(j)=-(duuudx(j)+duuvdy(j)+duuwdz(j))
         Tyy(j)=-(duvvdx(j)+dvvvdy(j)+dvvwdz(j))
         Tzz(j)=-(duwwdx(j)+dvwwdy(j)+dwwwdz(j))
         Txy(j)=-(duuvdx(j)+duvvdy(j)+duvwdz(j))
         Txz(j)=-(duuwdx(j)+duvwdy(j)+duwwdz(j))
         Tyz(j)=-(duvwdx(j)+dvvwdy(j)+dvwwdz(j))
      enddo

!     Derivatives of the pressure-velocity products
      call gradm1(dpudx,dpudy,dpudz,STAT(1,12))
      call gradm1(dpvdx,dpvdy,dpvdz,STAT(1,13))
      call gradm1(dpwdx,dpwdy,dpwdz,STAT(1,14))

!     Derivatives of the mean pressure field
      call gradm1(dPdx,dPdy,dPdz,P)
      
!     Pressure transport tensor. Tensor of Rank 2.
      do j=1,ntot
         dpudx(j)=dpudx(j)-P(j)*dUdx(j)-U(j)*dPdx(j)
         dpvdx(j)=dpvdx(j)-P(j)*dVdx(j)-V(j)*dPdx(j)
         dpwdx(j)=dpwdx(j)-P(j)*dWdx(j)-W(j)*dPdx(j)
         dpudy(j)=dpudy(j)-P(j)*dUdy(j)-U(j)*dPdy(j)
         dpvdy(j)=dpvdy(j)-P(j)*dVdy(j)-V(j)*dPdy(j)
         dpwdy(j)=dpwdy(j)-P(j)*dWdy(j)-W(j)*dPdy(j)
         dpudz(j)=dpudz(j)-P(j)*dUdz(j)-U(j)*dPdz(j)
         dpvdz(j)=dpvdz(j)-P(j)*dVdz(j)-V(j)*dPdz(j)
         dpwdz(j)=dpwdz(j)-P(j)*dWdz(j)-W(j)*dPdz(j)
      enddo

      do j=1,ntot
         PTxx(j)=-2./rho*dpudx(j)
         PTyy(j)=-2./rho*dpvdy(j)
         PTzz(j)=-2./rho*dpwdz(j)
         PTxy(j)=-1./rho*(dpudy(j)+dpvdx(j))
         PTxz(j)=-1./rho*(dpudz(j)+dpwdx(j))
         PTyz(j)=-1./rho*(dpvdz(j)+dpwdy(j))
      enddo

!     Pressure strain tensor. Tensor of Rank 2.
      do j=1,ntot
         pdudx(j)=STAT(j,15)-P(j)*dUdx(j)
         pdudy(j)=STAT(j,16)-P(j)*dUdy(j)
         pdudz(j)=STAT(j,17)-P(j)*dUdz(j)
         pdvdx(j)=STAT(j,18)-P(j)*dVdx(j)
         pdvdy(j)=STAT(j,19)-P(j)*dVdy(j)
         pdvdz(j)=STAT(j,20)-P(j)*dVdz(j)
         pdwdx(j)=STAT(j,21)-P(j)*dWdx(j)
         pdwdy(j)=STAT(j,22)-P(j)*dWdy(j)
         pdwdz(j)=STAT(j,23)-P(j)*dWdz(j)
      enddo

      do j=1,ntot
         PSxx(j)=-2./rho*pdudx(j)
         PSyy(j)=-2./rho*pdvdy(j)
         PSzz(j)=-2./rho*pdwdz(j)
         PSxy(j)=-1./rho*(pdudy(j)+pdvdx(j))
         PSxz(j)=-1./rho*(pdudz(j)+pdwdx(j))
         PSyz(j)=-1./rho*(pdvdz(j)+pdwdy(j))
      enddo

!     Velocity-pressure-gradient tensor. Tensor of Rank 2.
      do j=1,ntot
         Pixx(j)=PTxx(j)-PSxx(j)
         Piyy(j)=PTyy(j)-PSyy(j)
         Pizz(j)=PTzz(j)-PSzz(j)
         Pixy(j)=PTxy(j)-PSxy(j)
         Pixz(j)=PTxz(j)-PSxz(j)
         Piyz(j)=PTyz(j)-PSyz(j)
      enddo

!     Calculation of TKE budget
      do j=1,ntot
         Pk(j)=0.5*(Pxx(j)+Pyy(j)+Pzz(j))
         Dk(j)=0.5*(Dxx(j)+Dyy(j)+Dzz(j))
         Tk(j)=0.5*(Txx(j)+Tyy(j)+Tzz(j))
         VDk(j)=0.5*(VDxx(j)+VDyy(j)+VDzz(j))
         Pik(j)=0.5*(Pixx(j)+Piyy(j)+Pizz(j))
         Ck(j)=0.5*(Cxx(j)+Cyy(j)+Czz(j))
         Resk(j)=Pk(j)+Dk(j)+Tk(j)+VDk(j)+Pik(j)-Ck(j)
      enddo
      
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      ifpo=.FALSE.
      
      call outpost(U,V,W,pr,uu,'a01')
      call outpost(vv,ww,uv,pr,uw,'a02')
      call outpost(vw,P,pp,pr,ppp,'a03')
      call outpost(pppp,uuu,vvv,pr,www,'a04')
      call outpost(uuv,uuw,uvv,pr,vvw,'a05')
      call outpost(uww,vww,uvw,pr,Pxx,'a06')
      call outpost(Pyy,Pzz,Pxy,pr,Pxz,'a07')
      call outpost(Pyz,Dxx,Dyy,pr,Dzz,'a08')
      call outpost(Dxy,Dxz,Dyz,pr,Txx,'a09')
      call outpost(Tyy,Tzz,Txy,pr,Txz,'a10')
      call outpost(Tyz,VDxx,VDyy,pr,VDzz,'a11')
      call outpost(VDxy,VDxz,VDyz,pr,Pixx,'a12')
      call outpost(Piyy,Pizz,Pixy,pr,Pixz,'a13')
      call outpost(Piyz,Cxx,Cyy,pr,Czz,'a14')
      call outpost(Cxy,Cxz,Cyz,pr,Pk,'a15')
      call outpost(Dk,Tk,VDk,pr,Pik,'a16')
      call outpost(Ck,Resk,PTxx,pr,PTyy,'a17')
      call outpost(PTzz,PTxy,PTxz,pr,PTyz,'a18')
      call outpost(PSxx,PSyy,PSzz,pr,PSxy,'a19')
      call outpost(PSxz,PTyz,dUdx,pr,dUdy,'a20')
      call outpost(dUdz,dVdx,dVdy,pr,dVdz,'a21')
      call outpost(dWdx,dWdy,dWdz,pr,Tk,'a22')
      
      ifpo=.TRUE.
      
      return
      end


c-----------------------------------------------------------------------
      
      subroutine load_field(field)

      implicit none

! !      include 'SIZE_DEF'
      include 'SIZE'
      include 'USERPAR'                 ! l2freq, FIXGEOM, NEW_DT, STATS3D
      include 'STATS'           ! Statistics specific variables
! !      include 'INPUT_DEF'
      include 'INPUT'           ! if3d
! !      include 'SOLN_DEF'
      include 'SOLN'
! !      include 'TSTEP_DEF'
      include 'TSTEP'

      character*132 field
      
      call load_fld(field)
      
c      call outpost(vx,vy,vz,pr,t,'new')
      
      return
      end

c-----------------------------------------------------------------------
      subroutine read_hdr(field,mtimee)

      implicit none

! !      include 'SIZE_DEF'
      include 'SIZE'
      include 'USERPAR'                 ! l2freq, FIXGEOM, NEW_DT, STATS3D
      include 'STATS'           ! Statistics specific variables
! !      include 'INPUT_DEF'
      include 'INPUT'           ! if3d
! !      include 'SOLN_DEF'
      include 'SOLN'
! !      include 'TSTEP_DEF'
      include 'TSTEP'

      character*132 field,hdr,fmt1
      character*10 tmpchar
      integer twdsize,mnelx,mnely,mnelz,nelo,isteps,fid0,nfileoo
!      integer stat_gnum
      real mtimee

      open(unit=33,file=field,form='unformatted')
      read(33) hdr
      close(33)

      fmt1 = '(1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,
     &1x,i9,1x,i6,1x,i6,1x,10a)'

      read(hdr,fmt1) twdsize,mnelx,mnely,mnelz,nelo,
     &     stat_gnum,mtimee,isteps,fid0,nfileoo,tmpchar

      return
      end
c-----------------------------------------------------------------------
      
      subroutine usrdat()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      include 'ZPER'            ! For nelx,nely,nelz - needed for z_average

      return
      end
      
c-----------------------------------------------------------------------

      subroutine usrdat2()  
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
     
      return
      end

c-----------------------------------------------------------------------
      
      subroutine usrdat3()
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'

      return
      end
c---------------------------------------------------------------------- 

c automatically added by makenek
      subroutine usrsetvert(glo_num,nel,nx,ny,nz) ! to modify glo_num
      integer*8 glo_num(1)

      return
      end

c automatically added by makenek
      subroutine userqtl

      call userqtl_scig

      return
      end
