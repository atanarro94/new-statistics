c========================================================================
      subroutine uservp (ix,iy,iz,IEG)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      UDIFF = 0.0
      UTRANS = 0.0
      return
      end
c========================================================================
      subroutine userf  (ix,iy,iz,IEG)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      
      ffx = 0.0  ! Force through my_vol_flow
      ffy = 0.0
      ffz = 0.0

      return
      end
c========================================================================

      subroutine userq  (ix,iy,iz,IEG)
      include 'SIZE'
      include 'TOTAL'
      include 'NEKUSE'
      QVOL = 0.0
      SOURCE = 0.0
      return
      end
c========================================================================

      subroutine userchk
      include 'SIZE'
      include 'TOTAL'
      include 'RESTART'

!     Define  variables
      integer ss, ssf, ssfold
      integer nsmax
      parameter (nsmax = max(nstats,nderiv))
      real stat_2D1(lx1*lelx*ly1*lely,nsmax)
      real stat_2D (lx1*lelx*ly1*lely,nsmax)
      real stat_R(lx1*lelx*ly1*lely)

      real x_pts(lhis),y_pts(lhis)
      real xx(lhis), yy(lhis)      

      integer nelx, nely, ng 

      real mRey
      real mdomain_x, mdomain_y, mdomain_z
      real mtimes   , mtimee, mDT, matime_tot,matime
      real Ttimes   , Tatime, Tint      

      integer mnelx, mnely, mnelz
      integer polyx, polyy, polyz
      integer nstat, nrec , Tnrec 

      character*80 pippo
      character*132 inputname1, inputname2 
      character*132 hdr, inputname3
      character*80 val1, val2, val3, val4, val5, val6
      character*80 val7, val8, val9, val10, val11, val12, val13

!     Necessary to write averaged 2D fields in Nek
      real stat_temp1(lx1*lelx*ly1*lely)
      real stat_temp2(lx1*lelx*ly1*lely)

!     Variables related to 2D derivatives
      real duidxj(lx1*ly1*lz1,lelt,2*ldim)
      real ur(lx1*ly1*lz1),us(lx1*ly1*lz1)
      real vr(lx1*ly1*lz1),vs(lx1*ly1*lz1)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      parameter(nfldm=ldim+ldimt+1)

      real    pts(ldim,lhis)
      real    fieldout(lhis)

      real    dist(lhis)

      real    rst(lhis*ldim)
      integer rcode(lhis),elid(lhis),proc(lhis)
      common /hpts_r/ rst
      common /hpts_i/ rcode,elid,proc

      common /scrcg/  pm1 (lx1,ly1,lz1,lelv) ! mapped pressure
      common /outtmp/ wrk(lx1*ly1*lz1*lelt,nfldm)

      logical iffind

      integer icalld,npoints
      save    icalld,npoints
      data    icalld  /0/
      data    npoints /0/

      integer icall2
      save    icall2
      data    icall2 /-9/

      save    inth_hpts

!----------------------------------------------------------------------
!     For new version of statistics 
      integer nelo,stat_gnum,isteps,nfileoo,twdsize
      integer ierr,mfid

      integer fld_new(nstats)

      character*10 tmpchar
      character*100 fmt1

      logical IFINTP 

!----------------------------------------------------------------------  

      ifintp = .TRUE.

!      write(*,*) ifgetu,ifgetx,ifgetp,ifgtps
      ifgetx=.TRUE.
      do j=1,ldimt
          ifgtps(j)=.TRUE.
      enddo

!     Total number of grid points per element
      nxyz  = nx1*ny1*nz1

!     Total number of grid points (we run in 1 core)
      ntot  = nxyz*nelt

!     Number of elements in x and y directions
      nelx = lelx      
      nely = lely      

!     Number of stat files to process
      nfiles=param(68)

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
!     Read stat files and average them. Avearaging is based on the 
!     number of statistically independent samples per stat file
!     or number of records nrec
      
!     Number of grid-points in x, y and in the 2D plane
      mx=nx1*nelx
      my=ny1*nely
      m=mx*my

!     Initialize effective averaging time nrec*dt, time interval
!     and number of records
      Tatime = 0.
      Tint   = 0.
      Tnrec  = 0

!     Initialize accumulated (_2D1) and averaged 2D statistics
 
      call rzero(stat_2D1,m*nsmax)

!      do ifld=1,nstats
!         do j=1,m
!            stat_2D1(j,ifld) = 0.
!         enddo
!      enddo

      call rzero(stat_2D,m*nsmax)

!      do ifld=1,nftot
!         do j=1,m
!            stat_2D(j,ifld) = 0.
!         enddo
!      enddo

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!           Go through all stat files and read them
      mtimes = 0.0;   ! #change. Right now this is the start time

      do ssf = 1,nfiles
         ! Define name of corresponding stat file
         write(pippo,'(i5.5)') ssf
         inputname1 = 'ZSTAT/stsper_hill0.f'//trim(pippo)
         call load_fld(inputname1)
          open(unit=33,file=inputname1,form='unformatted')
          read(33) hdr
          close(33)
!         just get the header
        fmt1 = '(1x,i1,1x,i2,1x,i2,1x,i2,1x,i10,1x,i10,1x,e20.13,
     &1x,i9,1x,i6,1x,i6,1x,10a)'
!           write(*,*) fmt1
          read(hdr,fmt1) twdsize,mnelx,mnely,mnelz,nelo,
     $          stat_gnum,mtimee,isteps,fid0,nfileoo,tmpchar 
          read(tmpchar,'(2x,i2)') nstat
!                  
!          write(*,*) mnelx, mnely, mnelz, nelo, stat_gnum, mtimee,
!     $         isteps,fid0,nfileoo,tmpchar,nstat

          call my_stat_file(inputname3,ssf)
          call mfi(inputname3,33)
!          call full_restart(inputname3,1)        ! this interpolates data to the old mesh

          call geom_reset(1)
!          call gengeom(2)
!      Find inverse (jacmi) of the Jacobian array (jacm1)
          if (istep.ne.icall2) then
             call invers2(jacmi,jacm1,ntot)
             icall2=istep
          endif
  
          polyx = mnelx-1
          polyy = mnely-1
          polyz = mnelz-1

          matime_tot = (mtimee - mtimes)          
          matime = matime_tot/10
          mDT = 0.001
          nrec = isteps/10

          mdomain_x = 9.0 
          mdomain_y = 3.035
          mdomain_z = 4.5
          mnelx = nelo 
          mnely = 1
          mnelz = 1
          mRey  = 700. 

         ! Read stat file parameters
!        read(33) mRey,                        ! Reynolds number
!    &        mdomain_x, mdomain_y, mdomain_z, ! domain size
!    &        mnelx    , mnely    , mnelz,     ! number of elements 
!    &        polyx    , polyy    , polyz,     ! polynomial order
!    &        nstat,                           ! number of saved statistics 
!    &        mtimes,                          ! start time 
!    &        mtimee,                          ! end time
!    &        matime_tot,                      ! time interval
!    &        matime,                          ! effective averge time
!    &        mDT,                             ! time step
!    &        nrec                             ! number of time records

       
         ! Print starting, ending, total time of the stat file,
         ! as well as the number of records
         write(*,*) 'Ts,Te,Tint,Ta,nv,dt,nrec',
     &        mtimes,mtimee,matime_tot,matime,matime_tot/matime,mDT,nrec

         ! Update total effective time, time interval and
         ! total number of records
         Tatime = Tatime + matime
         Tint   = Tint + matime_tot
         Tnrec  = Tnrec  + nrec

         ! If it is the first stat file, its start time will be the
         ! start time of the averaged field
         if(ssf.eq.1) then
            Ttimes = mtimes
         endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

         ! If the number of fields in the stat files does not match
         ! the one allocated for postprocessing, exit
      if (nstat.ne.nstats) then
            write(6,*) 'Change nstats in SIZE to be equal to nstat'
            call exitt 
      endif 

         ! If the x-number of elements in the stat files does not
         ! match the one allocated for postprocessing, exit
      if (nelx.ne.mnelx) then
         write(*,*) nelx, mnelx
            write(6,*) 'change lelx to be equal to that of the run'
            call exitt
      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C
         ! Go through all the fields from the current stat file
!
!********** Shuffle arrangement of variables *************************  
!********** to match old post processing implementation ************** 

      do ifld=1,nstats
           fld_new(ifld)=ifld
      enddo

!         No change for first 26 variables
      do ifld=1,26
           fld_new(ifld)=ifld
      enddo
!         Changes here          
      do ifld=27,32
           fld_new(ifld) = ifld+1
      enddo
      fld_new(33) = 27
      fld_new(34) = 38
      do ifld=35,38
           fld_new(ifld) = ifld-1
      enddo
!         No changes from 39-44
      do ifld=39,44
           fld_new(ifld) = ifld
      enddo

      write(6,*) 'Stats order re-shuffled'

!********** End of Shuffling ******************************************   
!********************************************************************** 
               
         do ifld=1,nstat
            ! Read field and store it in temporary variable
               call copy(stat_R,t(1,1,1,1,ifld+1),m) 
            ! Update accumulated sum of that particular field.
            ! Fields are weighted with their time interval
            do j=1,m
               stat_2D1(j,fld_new(ifld)) = stat_2D1(j,fld_new(ifld)) +
     $              matime_tot*stat_R(j)
            enddo
         enddo
            mtimes = mtimee 
      ! finish loop over all stat files
      enddo  

      write(6,*) 'Accumulation of stats done'
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!     Divide accumulated fields by the total time interval
      do ifld=1,nstat
         do j=1,m
            stat_2D(j,ifld) = stat_2D1(j,ifld)/Tint
         enddo
      enddo

!     Print total number of records and time interval considered
!     for the statistics
      write(*,*) 'Total number of records:', Tnrec
      write(*,*) 'Total time interval:', Tint

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!     Write averaged 2D fields in Nek
!     Store fields 1 and 2 in stat_temp1 and stat_temp2
      do j=1,m
         stat_temp1(j)=stat_2D(j,1)
         stat_temp2(j)=stat_2D(j,2)
      enddo

!     Copy stat_temp1 and stat_temp2 to the vx and vy array
!      call opcopy(vx,vy,pr,stat_temp1,stat_temp2,stat_temp1)
      call copy(vx,stat_temp1,m)
      call copy(vy,stat_temp2,m)

!     Outpost velocity field, to visualize stat_temp1 and stat_temp2
      call outpost(vx,vy,vz,pr,t,'UV_')

!     Store fields 3 and 4 in stat_temp1 and stat_temp2
      do j=1,m
         stat_temp1(j)=stat_2D(j,3)
         stat_temp2(j)=stat_2D(j,4)
      enddo

!     Copy stat_temp1 and stat_temp2 to the vx and vy array
!      call opcopy(vx,vy,pr,stat_temp1,stat_temp2,stat_temp1)
      call copy(vx,stat_temp1,m)
      call copy(vy,stat_temp2,m)

!     Outpost velocity field, to visualize stat_temp1 and stat_temp2
      call outpost(vx,vy,vz,pr,t,'WP_')

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

      call rzero(stat_2D1,m*nsmax) 
      call copy(stat_2D1,stat_2D,m*nstat)
      call rzero(stat_2D,m*nsmax)


!     Compute 2D derivatives

c-----------------------------------------------------------------------

!     Tensor 1. Compute 2D derivative tensor with U and V
      call comp_derivat(duidxj,stat_2D1(1,1),stat_2D1(1,2),ur,us,vr,vs)
      call deriv1(stat_2D(1, 1),duidxj,ntot,'e11')                      ! dU/dx
      call deriv1(stat_2D(1,2),duidxj,ntot,'e12')                       ! dU/dy
      call deriv1(stat_2D(1,3),duidxj,ntot,'e21')                       ! dV/dx
      call deriv1(stat_2D(1,4),duidxj,ntot,'e22')                       ! dV/dy

c-----------------------------------------------------------------------

!     Tensor 2. Compute 2D derivative tensor with W and P
      call comp_derivat(duidxj,stat_2D1(1,3),stat_2D1(1,4),ur,us,vr,vs)
      call deriv1(stat_2D(1,5),duidxj,ntot,'e11')                       ! dW/dx
      call deriv1(stat_2D(1,6),duidxj,ntot,'e12')                       ! dW/dy
      call deriv1(stat_2D(1,7),duidxj,ntot,'e21')                       ! dP/dx
      call deriv1(stat_2D(1,8),duidxj,ntot,'e22')                       ! dP/dy

c-----------------------------------------------------------------------

!     Tensor 3. Compute 2D derivative tensor with <uu> and <vv>
      call comp_derivat(duidxj,stat_2D1(1,5),stat_2D1(1,6),ur,us,vr,vs)
      call deriv1(stat_2D(1,9),duidxj,ntot,'e11')                       ! d<uu>/dx
      call deriv1(stat_2D(1,10),duidxj,ntot,'e12')                      ! d<uu>/dy
      call deriv1(stat_2D(1,11),duidxj,ntot,'e21')                      ! d<vv>/dx
      call deriv1(stat_2D(1,12),duidxj,ntot,'e22')                      ! d<vv>/dy

c-----------------------------------------------------------------------

!     Tensor 4. Compute 2D derivative tensor with <ww> and <pp>
      call comp_derivat(duidxj,stat_2D1(1,7),stat_2D1(1,8),ur,us,vr,vs)
      call deriv1(stat_2D(1,13),duidxj,ntot,'e11')                      ! d<ww>/dx
      call deriv1(stat_2D(1,14),duidxj,ntot,'e12')                      ! d<ww>/dy
      call deriv1(stat_2D(1,15),duidxj,ntot,'e21')                      ! d<pp>/dx
      call deriv1(stat_2D(1,16),duidxj,ntot,'e22')                      ! d<pp>/dy

c-----------------------------------------------------------------------

!     Tensor 5. Compute 2D derivative tensor with <uv> and <vw>
      call comp_derivat(duidxj,stat_2D1(1,9),stat_2D1(1,10),ur,us,vr,vs)
      call deriv1(stat_2D(1,17),duidxj,ntot,'e11')                      ! d<uv>/dx
      call deriv1(stat_2D(1,18),duidxj,ntot,'e12')                      ! d<uv>/dy
      call deriv1(stat_2D(1,19),duidxj,ntot,'e21')                      ! d<vw>/dx
      call deriv1(stat_2D(1,20),duidxj,ntot,'e22')                      ! d<vw>/dy

c-----------------------------------------------------------------------

!     Tensor 6. Compute 2D derivative tensor with <uw> and <uuu>
      call comp_derivat(duidxj,stat_2D1(1,11),stat_2D1(1,24),
     &                   ur,us,vr,vs)
      call deriv1(stat_2D(1,21),duidxj,ntot,'e11')                      !  d<uw>/dx
      call deriv1(stat_2D(1,22),duidxj,ntot,'e12')                      !  d<uw>/dy
      call deriv1(stat_2D(1,23),duidxj,ntot,'e21')                      ! d<uuu>/dx
      call deriv1(stat_2D(1,24),duidxj,ntot,'e22')                      ! d<uuu>/dy

c-----------------------------------------------------------------------

!     Tensor 7. Compute 2D derivative tensor with <vvv> and <www>
      call comp_derivat(duidxj,stat_2D1(1,25),stat_2D1(1,26),
     &                   ur,us,vr,vs)
      call deriv1(stat_2D(1,25),duidxj,ntot,'e11')                      ! d<vvv>/dx
      call deriv1(stat_2D(1,26),duidxj,ntot,'e12')                      ! d<vvv>/dy
      call deriv1(stat_2D(1,27),duidxj,ntot,'e21')                      ! d<www>/dx
      call deriv1(stat_2D(1,28),duidxj,ntot,'e22')                      ! d<www>/dy

c-----------------------------------------------------------------------
!!
!     Tensor 8. Compute 2D derivative tensor with <ppp> and <uuv>
      call comp_derivat(duidxj,stat_2D1(1,27),stat_2D1(1,28),
     &                   ur,us,vr,vs)
      call deriv1(stat_2D(1,29),duidxj,ntot,'e11')                      ! d<ppp>/dx
      call deriv1(stat_2D(1,30),duidxj,ntot,'e12')                      ! d<ppp>/dy
      call deriv1(stat_2D(1,31),duidxj,ntot,'e21')                      ! d<uuv>/dx
      call deriv1(stat_2D(1,32),duidxj,ntot,'e22')                      ! d<uuv>/dy

c-----------------------------------------------------------------------
!!
!     Tensor 9. Compute 2D derivative tensor with <uuw> and <vvu>
      call comp_derivat(duidxj,stat_2D1(1,29),stat_2D1(1,30),
     &                   ur,us,vr,vs)
      call deriv1(stat_2D(1,33),duidxj,ntot,'e11')                      ! d<uuw>/dx
      call deriv1(stat_2D(1,34),duidxj,ntot,'e12')                      ! d<uuw>/dy
      call deriv1(stat_2D(1,35),duidxj,ntot,'e21')                      ! d<vvu>/dx
      call deriv1(stat_2D(1,36),duidxj,ntot,'e22')                      ! d<vvu>/dy

c-----------------------------------------------------------------------
!!
!     Tensor 10. Compute 2D derivative tensor with <vvw> and <wwu>
      call comp_derivat(duidxj,stat_2D1(1,31),stat_2D1(1,32),
     &                   ur,us,vr,vs)
      call deriv1(stat_2D(1,37),duidxj,ntot,'e11')                      ! d<vvw>/dx
      call deriv1(stat_2D(1,38),duidxj,ntot,'e12')                      ! d<vvw>/dy
      call deriv1(stat_2D(1,39),duidxj,ntot,'e21')                      ! d<wwu>/dx
      call deriv1(stat_2D(1,40),duidxj,ntot,'e22')                      ! d<wwu>/dy

c-----------------------------------------------------------------------
!!
!     Tensor 11. Compute 2D derivative tensor with <wwv> and <uvw>
      call comp_derivat(duidxj,stat_2D1(1,33),stat_2D1(1,34),
     &                   ur,us,vr,vs)
      call deriv1(stat_2D(1,41),duidxj,ntot,'e11')                      ! d<wwv>/dx
      call deriv1(stat_2D(1,42),duidxj,ntot,'e12')                      ! d<wwv>/dy
      call deriv1(stat_2D(1,43),duidxj,ntot,'e21')                      ! d<uvw>/dx
      call deriv1(stat_2D(1,44),duidxj,ntot,'e22')                      ! d<uvw>/dy

c-----------------------------------------------------------------------

!     Tensor 12. Compute 2D derivative tensor with dU/dx and dU/dy
      call comp_derivat(duidxj,stat_2D(1,1),stat_2D(1,2),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,45),duidxj,ntot,'e11')                      ! d2U/dx2
      call deriv1(stat_2D(1,46),duidxj,ntot,'e22')                      ! d2U/dy2

c-----------------------------------------------------------------------

!     Tensor 13. Compute 2D derivative tensor with dV/dx and dV/dy
      call comp_derivat(duidxj,stat_2D(1,3),stat_2D(1,4),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,47),duidxj,ntot,'e11')                      ! d2V/dx2
      call deriv1(stat_2D(1,48),duidxj,ntot,'e22')                      ! d2V/dy2

c-----------------------------------------------------------------------

!     Tensor 14. Compute 2D derivative tensor with dW/dx and dW/dy
      call comp_derivat(duidxj,stat_2D(1,5),stat_2D(1,6),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,49),duidxj,ntot,'e11')                      ! d2W/dx2
      call deriv1(stat_2D(1,50),duidxj,ntot,'e22')                      ! d2W/dy2

c-----------------------------------------------------------------------

!     Tensor 15. Compute 2D derivative tensor with d<uu>/dx and d<uu>/dy
      call comp_derivat(duidxj,stat_2D(1,9),stat_2D(1,10),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,51),duidxj,ntot,'e11')                      ! d2<uu>/dx2
      call deriv1(stat_2D(1,52),duidxj,ntot,'e22')                      ! d2<uu>/dy2

c-----------------------------------------------------------------------

!     Tensor 16. Compute 2D derivative tensor with d<vv>/dx and d<vv>/dy
      call comp_derivat(duidxj,stat_2D(1,11),stat_2D(1,12),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,53),duidxj,ntot,'e11')                      ! d2<vv>/dx2
      call deriv1(stat_2D(1,54),duidxj,ntot,'e22')                      ! d2<vv>/dy2

c-----------------------------------------------------------------------

!     Tensor 17. Compute 2D derivative tensor with d<ww>/dx and d<ww>/dy
      call comp_derivat(duidxj,stat_2D(1,13),stat_2D(1,14),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,55),duidxj,ntot,'e11')                      ! d2<ww>/dx2
      call deriv1(stat_2D(1,56),duidxj,ntot,'e22')                      ! d2<ww>/dy2

c-----------------------------------------------------------------------

!     Tensor 18. Compute 2D derivative tensor with d<uv>/dx and d<uv>/dy
      call comp_derivat(duidxj,stat_2D(1,17),stat_2D(1,18),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,57),duidxj,ntot,'e11')                      ! d2<uv>/dx2
      call deriv1(stat_2D(1,58),duidxj,ntot,'e22')                      ! d2<uv>/dy2

c-----------------------------------------------------------------------

!     Tensor 19. Compute 2D derivative tensor with d<uw>/dx and d<uw>/dy
      call comp_derivat(duidxj,stat_2D(1,21),stat_2D(1,22),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,59),duidxj,ntot,'e11')                      ! d2<uw>/dx2
      call deriv1(stat_2D(1,60),duidxj,ntot,'e22')                      ! d2<uw>/dy2

c-----------------------------------------------------------------------

!     Tensor 20. Compute 2D derivative tensor with d<vw>/dx and d<vw>/dy
      call comp_derivat(duidxj,stat_2D(1,19),stat_2D(1,20),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,61),duidxj,ntot,'e11')                      ! d2<vw>/dx2
      call deriv1(stat_2D(1,62),duidxj,ntot,'e22')                      ! d2<vw>/dy2

c-----------------------------------------------------------------------

!     Tensor 21. Compute 2D derivative tensor with <pu> and <pv>
      call comp_derivat(duidxj,stat_2D1(1,12),stat_2D1(1,13),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,63),duidxj,ntot,'e11')                      ! d<pu>/dx
      call deriv1(stat_2D(1,64),duidxj,ntot,'e12')                      ! d<pu>/dy
      call deriv1(stat_2D(1,65),duidxj,ntot,'e21')                      ! d<pv>/dx
      call deriv1(stat_2D(1,66),duidxj,ntot,'e22')                      ! d<pv>/dy

c-----------------------------------------------------------------------

!     Tensor 22. Compute 2D derivative tensor with <pw> and 
      call comp_derivat(duidxj,stat_2D1(1,14),stat_2D1(1,1),
     &     ur,us,vr,vs)
      call deriv1(stat_2D(1,67),duidxj,ntot,'e11')                      ! d<pw>/dx
      call deriv1(stat_2D(1,68),duidxj,ntot,'e12')                      ! d<pw>/dy

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!    RMS terms
!!    u_rms
!      do j=1,m
!          stat_R(j)=stat_2D(j,1)**2
!      enddo
!     
!      call sub3(stat_2D(1,nstat+69),stat_2D(1,5),stat_R,m)
!!      call copy(stat_2D(1,nstat+69),stat_R,m)
!
!      do j=1,m
!          stat_2D(j,nstat+69)=sqrt(abs(stat_2D(j,nstat+69)))
!      enddo
!
!!      call vsqrt(stat_2D(1,nstat+69),m)
!
!!---------------------------------------------------------------------- 
!
!!    v_rms
!      do j=1,m
!          stat_R(j)=stat_2D(j,2)**2
!      enddo
!     
!      call sub3(stat_2D(1,nstat+70),stat_2D(1,6),stat_R,m)
!
!      do j=1,m
!          stat_2D(j,nstat+70)=sqrt(abs(stat_2D(j,nstat+70)))
!      enddo
!
!
!!---------------------------------------------------------------------- 
!
!!    w_rms
!      do j=1,m
!          stat_R(j)=stat_2D(j,3)**2
!      enddo
!     
!      call sub3(stat_2D(1,nstat+71),stat_2D(1,7),stat_R,m)
!      do j=1,m
!          stat_2D(j,nstat+71)=sqrt(abs(stat_2D(j,nstat+71)))
!      enddo
!
!
C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C


      ifto = .FALSE.

      call copy(T(1,1,1,1,2),stat_2D1,m*nstats)
      call outpost2(vx,vy,vz,pr,T,nstats+1,'st1')

      call copy(T(1,1,1,1,2),stat_2D,m*nderiv)
      call outpost2(vx,vy,vz,pr,T,nderiv+1,'st2')

      if (ifintp) then

!     Read the grid in the x direction
      open(unit=20,form='unformatted',file='ZSTAT/x.fort')
      read(20) ng
      read(20) (x_pts(i),i=1,ng)
      close(20)

      write(*,*) 'ngx',ng

!     Read the grid in the y direction
      open(unit=212,form='unformatted',file='ZSTAT/y.fort')
      read(212) ng
      read(212) (y_pts(i),i=1,ng)
      close(212)

      write(*,*) 'ngy',ng

!     Generate mesh for interpolation
      do i=1,ng
         npoints=npoints+1
         pts(1,npoints)=x_pts(i)
         pts(2,npoints)=y_pts(i)
      enddo

      write(*,*) 'Check npoints',ng,npoints

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!     Generate file to store the interpolated and averaged field
!     Name of the file is int_fld
      inputname2 = 'ZSTAT/int_fld'

!     Open / Create int_fld file
      open(unit=37,form='unformatted',file=inputname2) 

!     Parameters to write on the header of the int_fld file
      write(val1,'(1p15e17.9)') mRey                          ! Reynolds number	  
      write(val2,'(1p15e17.9)') mdomain_x,mdomain_y,mdomain_z ! domain size
      write(val3,'(9i9)') mnelx,mnely,mnelz                   ! number of elements 
      write(val4,'(9i9)') polyx,polyy,polyz                   ! polynomial order
      write(val5,'(9i9)')       nstat                         ! number of fields on stat files
      write(val6,'(9i9)')       nderiv                        ! number of derivative fields
      write(val7,'(1p15e17.9)') Ttimes                        ! start time of first stat file
      write(val8,'(1p15e17.9)') mtimee                        ! end time of last stat file
      write(val9,'(1p15e17.9)') Tatime                        ! effective average time
      write(val10,'(1p15e17.9)') mDT                          ! time step
      write(val11,'(9i9)')      Tnrec                         ! total number of records
      write(val12,'(1p15e17.9)') Tint                         ! time interval
      write(val13,'(9i9)')      npoints                       ! number of points in int. mesh 


!     Write header
      write(37) '(Re ='//trim(val1)
     &     //') (Lx, Ly, Lz ='//trim(val2)
     &     //') (nelx, nely, nelz ='//trim(val3)
     &     //') (Polynomial order ='//trim(val4)
     &     //') (Nstat ='//trim(val5)
     &     //') (Nderiv ='//trim(val6)
     &     //') (start time ='//trim(val7)
     &     //') (end time ='//trim(val8)
     &     //') (effective average time ='//trim(val9)
     &     //') (time step ='//trim(val10)
     &     //') (nrec ='//trim(val11)
     &     //') (time interval ='//trim(val12)
     &     //') (npoints ='//trim(val13)
     &     //')'

!     Write values corresponding to the header
      write(37) mRey,
     &     mdomain_x, mdomain_y, mdomain_z,
     &     mnelx    , mnely    , mnelz,
     &     polyx    , polyy    , polyz,
     &     nstat,
     &     nderiv,
     &     Ttimes,
     &     mtimee,
     &     Tatime,
     &     mDT,
     &     Tnrec,
     &     Tint,
     &     npoints

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!     Do necessary checks before interpolation
      if(icalld.eq.0) then
         if(nid.eq.0) then
            ! Compare number of points in interpolating mesh npoints
            ! with maximum defined in SIZE file lhis
            if(npoints.gt.lhis) then
               write(6,*) 'Increase lhis to npoints',lhis,npoints
               call exitt
            endif
            write(6,*) 'found ', npoints, ' points in interp. mesh'
         endif 
         ! Setup for interpolation tool. Use default tolerance of -1
         call intpts_setup(-1.0,inth_hpts) 
      endif

!     Compare number of points in interpolating mesh npoints
!     with maximum defined in SIZE file lhis
      if(npoints.gt.lhis) then
         if(nid.eq.0) write(6,*) 
     &        'ABORT: lhis too low, increase in SIZE', npoints, lhis
         call exitt
      endif

!     Start interpolation
      if(icalld.eq.0) then

!     Conceptually, locate npt points. The data corresponding to each point
!     is whether it is inside an element, closesst to a border, not found.
!     Also identify the processor where the point was found, the element
!     where the point was found, parametric coordinates of the point and
!     the distance squared from found to sought point
         call fgslib_findpts(inth_hpts,rcode,1,
     &        proc,1,
     &        elid,1,
     &        rst,ndim,
     &        dist,1,
     &        pts(1,1),ndim,
     &        pts(2,1),ndim,
     &        pts(3,1),ndim,npoints)

         
!     Check the return code
         do i=1,npoints
            ! Interpolating point is on boundary or outside the SEM mesh
            if(rcode(i).eq.1) then
               if(dist(i).gt.1e-12) then
                  write(6,'(A,4E15.7)') 
     &                 ' WARNING: point on boundary or outside
     &the mesh xy[z]d^2:',
     &                 (pts(k,i),k=1,ndim),dist(i)
               endif   
            ! Interpolating ponit is not in the SEM mesh   
            elseif(rcode(i).eq.2) then
               nfail = nfail + 1
c               write(6,'(A,3E15.7)') 
c     &              ' WARNING: point not within mesh xy[z]: !',
c     &              (pts(k,i),k=1,ndim)
            endif
         enddo
         icalld = 1
      endif

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!     Once the points from the interpolating mesh are located in the
!     SEM mesh, we proceed to interpolate
      do ifld=1,nstat
         ! Evaluate the input field at the given points
         fieldout = 0.

         call fgslib_findpts_eval(inth_hpts,fieldout,1,
     &        rcode,1,
     &        proc,1,
     &        elid,1,
     &        rst,ndim,npoints,
     &        stat_2D1(1,ifld))         ! write NEK calculated fields
         
         ! Write out the interpolated field
         if(nid.eq.0) then
            write(37) (fieldout(i),i=1,npoints)
         endif

      enddo

      do ifld=1,nderiv
         ! Evaluate the input field at the given points
         fieldout = 0.

         call fgslib_findpts_eval(inth_hpts,fieldout,1,
     &        rcode,1,
     &        proc,1,
     &        elid,1,
     &        rst,ndim,npoints,
     &        stat_2D(1,ifld))          ! write post processed fields
         
         ! Write out the interpolated field
         if(nid.eq.0) then
            write(37) (fieldout(i),i=1,npoints)
         endif

      ! Conclude interpolating and writing fields         
      enddo

!     Close the int_fld file after finish writing
      close(37)
      write(6,*) 'Output interporlated field in int_fld' 

      endif    ! ifintp

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

!     Finish postprocessing code
      call exitt

      return
      end subroutine userchk

C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%C

c========================================================================

      subroutine userbc (ix,iy,iz,iside,IEG)
      include 'SIZE'
      include 'TSTEP'
      include 'NEKUSE'
      ux=1.
      uy=1.
      return
      end

c========================================================================

      subroutine useric (ix,iy,iz,IEG)
      include 'SIZE'
      include 'TSTEP'
      include 'NEKUSE'
      ux=1.
      uy=1.
      return
      end

c========================================================================

      subroutine usrdat

      include 'SIZE'
      include 'TOTAL'
      include 'ZPER'    

      return
      end

c========================================================================

      subroutine usrdat2
      include 'SIZE'
      include 'TOTAL'

c      call fix_mygll
      return
      end

c========================================================================

      subroutine usrdat3
      include 'SIZE'
      include 'TOTAL'
      return
      end

c-----------------------------------------------------------------------

      subroutine comp_derivat(duidxj,u,v,ur,us,vr,vs)
      include 'SIZE'
      include 'TOTAL'

      integer e

      real duidxj(lx1*ly1*lz1,lelt,2*ldim)    ! 4 terms
      real u  (lx1*ly1*lz1,lelt)
      real v  (lx1*ly1*lz1,lelt)
      real ur (1) , us (1) 
      real vr (1) , vs (1) 
c
c      common /dudxyj/ jacmi(lx1*ly1*lz1,lelt)
c      real jacmi
c
      n    = nx1-1                          ! Polynomial degree
      nxyz = nx1*ny1*nz1

      do e=1,nelv
         call local_grad2(ur,us,u,N,e,dxm1,dxtm1)
         call local_grad2(vr,vs,v,N,e,dxm1,dxtm1)

!     Derivative tensor computed by using the inverse of
!     the Jacobian array jacmi
      do k=1,nxyz
         duidxj(k,e,1) = jacmi(k,e)*(ur(k)*rxm1(k,1,1,e)+
     $        us(k)*sxm1(k,1,1,e))
         duidxj(k,e,2) = jacmi(k,e)*(ur(k)*rym1(k,1,1,e)+
     $        us(k)*sym1(k,1,1,e))
         duidxj(k,e,3) = jacmi(k,e)*(vr(k)*rxm1(k,1,1,e)+
     $        vs(k)*sxm1(k,1,1,e))
         duidxj(k,e,4) = jacmi(k,e)*(vr(k)*rym1(k,1,1,e)+
     $        vs(k)*sym1(k,1,1,e))

      enddo
      enddo

      return
      end

c-----------------------------------------------------------------------

      subroutine deriv1(avg,duidxj,n,name)
      include 'SIZE'

      integer n,k
      real duidxj(lx1*ly1*lz1,lelt,2*ldim)
      real avg(n)
      character*3 name

      if (name .eq. 'e11') then
         do k=1,n
            avg(k) = duidxj(k,1,1)
         enddo
      elseif (name .eq. 'e12') then
         do k=1,n
            avg(k) = duidxj(k,1,2)
         enddo
      elseif (name .eq. 'e21') then
         do k=1,n
            avg(k) = duidxj(k,1,3)
         enddo
      elseif (name .eq. 'e22') then
         do k=1,n
            avg(k) = duidxj(k,1,4)
         enddo
      endif

      return
      end

c-----------------------------------------------------------------------

      subroutine stat_mfi_prepare(hname,hdr)

      character*132 hname

      include 'SIZE'
      include 'PARALLEL'
      include 'RESTART'

      integer stride
      character*132 hdr
      logical if_byte_swap_test
      real*4 bytetest

      integer*8 offs0,offs

      integer sum

      ierr = 0
#ifndef MPIIO
      ! rank0 (i/o master) will do a pre-read to get some infos 
      ! we need to have in advance
      if (nid.eq.0) then
!          call mbyte_open(hname,0,ierr) ! open  blah000.fldnn
         call byte_open(hname,ierr)
         if(ierr.ne.0) goto 101
         call blank     (hdr,iHeaderSize)
         call byte_read (hdr, iHeaderSize/4,ierr)
          
         if(ierr.ne.0) goto 101
         call byte_read (bytetest,1,ierr)
         if(ierr.ne.0) goto 101
         if_byte_sw = if_byte_swap_test(bytetest,ierr) ! determine endianess
        
         if(ierr.ne.0) goto 101
         call mfi_parse_hdr(hdr,ierr)
      endif

 101  continue
      call err_chk(ierr,'Error reading restart header in mfi_prepare$')

      call bcast(if_byte_sw,lsize) 
      call bcast(hdr,iHeaderSize)  
      if(nid.ne.0) call mfi_parse_hdr(hdr,ierr)

      stride = np / nfiler
      if (stride.lt.1) then
         write(6,*) nfiler,np,'  TOO MANY FILES, mfi_prepare'
         call exitt
      endif

      if (mod(nid,stride).eq.0) then ! i/o clients
         pid0r = nid
         pid1r = nid + stride
         fid0r = nid / stride
         if (nid.ne.0) then ! don't do it again for rank0
            call blank     (hdr,iHeaderSize)
!             call mbyte_open(hname,fid0r,ierr) ! open  blah000.fldnn
            call byte_open(hname,ierr)
            if(ierr.ne.0) goto 102
            call byte_read (hdr, iHeaderSize/4,ierr)  
            if(ierr.ne.0) goto 102
            call byte_read (bytetest,1,ierr) 
            if(ierr.ne.0) goto 102
            call mfi_parse_hdr (hdr,ierr)  ! replace hdr with correct one 
         endif
         call byte_read (er,nelr,ierr)     ! get element mapping
         if (if_byte_sw) call byte_reverse(er,nelr,ierr)
      endif
#else
      pid0r = nid
      pid1r = nid
      offs0 = iHeaderSize + 4
!      call mbyte_open(hname,0,ierr)
      call byte_open(hname,ierr)
      ierr=iglmax(ierr,1)
      if(ierr.ne.0) goto 103

      call byte_read_mpi(hdr,iHeaderSize/4,pid00,ifh_mbyte,ierr)
      ierr=iglmax(ierr,1)
      if(ierr.ne.0) goto 103

      call byte_read_mpi(bytetest,1,pid00,ifh_mbyte,ierr)

 103  continue 
      call err_chk(ierr,'Error reading header/element map.$')
      
      call bcast(hdr,iHeaderSize) 
      call bcast(bytetest,4) 

      if_byte_sw = if_byte_swap_test(bytetest,ierr) ! determine endianess
      call mfi_parse_hdr(hdr,ierr)
      if(nfiler.ne.1) then
        if(nid.eq.0) write(6,*) 'ABORT: too many restart files!'
        call exitt
      endif
      nfiler = np

      ! number of elements to read 
      nelr = nelgr/np
      do i = 0,mod(nelgr,np)-1
         if(i.eq.nid) nelr = nelr + 1
      enddo
      nelBr = igl_running_sum(nelr) - nelr 
      offs = offs0 + nelBr*isize

      call byte_set_view(offs,ifh_mbyte)
      call byte_read_mpi(er,nelr,-1,ifh_mbyte,ierr)
      if (if_byte_sw) call byte_reverse(er,nelr,ierr)
#endif
 102  continue
      call err_chk(ierr,'Error reading header/element map.$')


      return
      end subroutine stat_mfi_prepare
!---------------------------------------------------------------------- 

      subroutine my_stat_file(inputfile1,ssf)

      include 'SIZE'            ! NID
      include 'TSTEP'           ! IOSTEP, ISTEP
      include 'INPUT'           ! SESSION, IFREGUO
      include 'RESTART'         ! IFRIRO, NFILEO
      include 'PARALLEL'        ! ISIZE

      integer nfil, len, k, ndigit, i, ierr, ssf
      real time_l, rfileo

      character*132 inputfile1
      character*132 fname
      character*1 fnam1(132)
      equivalence (fnam1,fname)

      character*6  six,str
      save         six
      data         six / "??????" /

      character*1 slash,dot
      save        slash,dot
      data        slash,dot  / '/' , '.' /

      character*17 kst
      save         kst
      data         kst / '0123456789abcdefx' /
      character*1  ks1(0:16)
      equivalence (ks1,kst)

      integer i_set, nrsf
      save    i_set
      parameter  (nrsf=1)

!      logical ifrestart
!      save    ifrestart
      character*80 s80, rstfile
      save         s80,rstfile

c     create stsfile name (SESSION.restart)
         call blank(fname,132)
         call blank(rstfile,80)
         
c     create names acording to mfo_open_files
          len=0
          call bcast(len ,ISIZE)
          i_set=len
          
          call blank(fname,132)

#ifndef NOMPIIO
          rfileo = 1
#else
          rfileo = NFILEO
#endif
          ndigit = log10(rfileo) + 1
  
          IFDIRO = .TRUE.
          k = 1
          if (IFDIRO) then !  Add directoy
             call chcopy(fnam1(1),'ZSTAT',5)
             k=k+5
!             call chcopy(fnam1(6),six,ndigit) ! put ???? in string
!             k = k + ndigit
             call chcopy(fnam1(k),slash,1)
             k = k+1
          endif
  
          call chcopy(fnam1(k),'sts',3) !  Add prefix
          k = k+3

          len = ltrunc(SESSION,132) !  Add SESSION
!          len = len-2               !  remove "2D"
          call chcopy(fnam1(k),SESSION,len)
          k = k+len

!          if (IFREGUO) then
!             len=4
!             call chcopy(fnam1(k),'_reg',len)
!             k = k+len
!          endif
          write(*,*) 'k',k, ndigit
          call chcopy(fnam1(k),six,ndigit) !  Add file-id holder
          k = k + ndigit
          write(*,*) 'TEST0', k, ndigit
          call chcopy(fnam1(k),dot,1) !  Add .f appendix
          write(*,*) 'TEST1'
          call chcopy(fnam1(k+1),'f',1)
          k = k + 2
c     is fname too long?
         if ((k+5).gt.80) then
            if(NID.eq.0) write(6,*) 'my_stat_file: k too big'
            call exitt
         endif
c     fill array with file names
         call blank(s80,80)
         write(str,54) nrsf*i_set+ssf
54       format(i5.5)
         call chcopy(fnam1(k),str,5)
         call chcopy(s80,fname,k+5)

!         write(*,*) 'filename:', fname
!         write(*,*) 'filename:', s80

         inputfile1 = s80
      return
      end subroutine my_stat_file

c-----------------------------------------------------------------------
cc copied from nek1093/trunk/nek/postpro.f
c-----------------------------------------------------------------------
      subroutine intpts_setup(tolin,ih)
c
c setup routine for interpolation tool
c tolin ... stop point seach interation if 1-norm of the step in (r,s,t) 
c           is smaller than tolin 
c
      include 'SIZE'
      include 'GEOM'

      common /nekmpi/ nidd,npp,nekcomm,nekgroup,nekreal

      tol = tolin
      if (tolin.lt.0) tol = 1e-13 ! default tolerance 

      n       = lx1*ly1*lz1*lelt 
      npt_max = 256
      nxf     = 2*nx1 ! fine mesh for bb-test
      nyf     = 2*ny1
      nzf     = 2*nz1
      bb_t    = 0.1 ! relative size to expand bounding boxes by
c
      if(nidd.eq.0) write(6,*) 'initializing intpts(), tol=', tol
      call fgslib_findpts_setup(ih,nekcomm,npp,ndim,
     &                     xm1,ym1,zm1,nx1,ny1,nz1,
     &                     nelt,nxf,nyf,nzf,bb_t,n,n,
     &                     npt_max,tol)
c       
      return
      end
c-----------------------------------------------------------------------
      
!====================================================================== 

c
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
