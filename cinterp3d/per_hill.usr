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
      include 'SIZE'
      include 'GEOM'                    ! xm1, ym1, zm1
      include 'SOLN'                    ! T
      include 'MASS'                    !BM1 for lambda2
      include 'TSTEP'                   ! ISTEP
      include 'INPUT'                   ! PARAM(12) (DT)
      include 'USERPAR'                 ! l2freq, FIXGEOM, NEW_DT
      include 'STATS'
      include 'mpif.h'


      character*132 inputname1,hdr,pippo,pippa
      character*132 inputnameu,inputnamew,inputnamel,inputnamev
      character*80 val1,val2
      integer ntot,nfiles,sfn,nface,e,f,i,j,count,countb
      
      real ftime
      integer n_g,ngx,ngy,ngz
      real buff_pts(IntSize),pts(ldim,lhis)
      integer nnhis,npr,i_point,i_buff,bnhis
      
      real    dist(lhis)
      real    rst(lhis*ldim)
      integer rcode(lhis),elid(lhis),proc(lhis)
      common /hpts_r/ rst
      common /hpts_i/ rcode,elid,proc
      integer icalld
      save    icalld
      data    icalld  /0/
      integer inth_hpts
      save    inth_hpts
      integer nfail
      integer i_field
      real    buff_fieldout(lhis,n_out_fields)
      real struct(lx1*ly1*lz1*lelt,n_out_fields)

      real    b_u(lhis),b_w(lhis),b_l(lhis),b_v(lhis)
      real    bb_u(lhis),bb_w(lhis),bb_l(lhis),bb_v(lhis)
            
      common /nekmpi/ nid_,np_,nekcomm,nekgroup,nekreal
      integer dest,tag0,rnk,ierr
      integer tag1,tag2,tag3,tag4
      integer status(mpi_status_size)

      nfiles=param(68)          ! Number of stat tiles
      ntot=lx1*ly1*lz1*lelt     ! Number of grid points per element
      ftime=0.0

!     Check distribution of points
      if (nid.eq.0) write(*,*) 'IntSize,lhis*lp',IntSize,lp*lhis

!     Read interpolating mesh in x
      open(unit=20,form='unformatted',file='ZSTAT/x.fort')
      read(20) ngx
      read(20) (buff_pts(i),i=1,ngx)
      close(20)

      n_g=ngx
      if (mod(n_g,lp).eq.0) then
         nnhis=n_g/lp
         npr=0
      else
         npr=lp-mod(n_g,lp)
         if (nid.lt.lp-npr) then
            nnhis = n_g/lp+1
         else
            nnhis = n_g/lp
         endif
      endif

      do i_point=1,nnhis
         if (mod(n_g,lp).eq.0) then
            i_buff = nnhis*nid+i_point
         else
            if (nid.lt.lp-npr) then
               i_buff = nnhis*nid+i_point
            else
               i_buff = nnhis*nid+lp-npr+i_point
            endif
         endif
         pts(1,i_point) = buff_pts(i_buff)
      enddo

!     Read interpolating mesh in y
      open(unit=21,form='unformatted',file='ZSTAT/y.fort')
      read(21) ngy
      read(21) (buff_pts(i),i=1,ngy)
      close(21)

      do i_point=1,nnhis
         if (mod(n_g,lp).eq.0) then
            i_buff = nnhis*nid+i_point
         else
            if (nid.lt.lp-npr) then
               i_buff = nnhis*nid+i_point
            else
               i_buff = nnhis*nid+lp-npr+i_point
            endif
         endif
         pts(2,i_point) = buff_pts(i_buff)
      enddo

!     Read interpolating mesh in z
      open(unit=22,form='unformatted',file='ZSTAT/z.fort')
      read(22) ngz
      read(22) (buff_pts(i),i=1,ngz)
      close(22)

      do i_point=1,nnhis
         if (mod(n_g,lp).eq.0) then
            i_buff = nnhis*nid+i_point
         else
            if (nid.lt.lp-npr) then
               i_buff = nnhis*nid+i_point
            else
               i_buff = nnhis*nid+lp-npr+i_point
            endif
         endif
         pts(3,i_point) = buff_pts(i_buff)
      enddo

      do sfn=1,nfiles
         write(pippa,'(i2.2)') sfn
         inputname1 = 'ZSTAT/a'//trim(pippa)//'per_hill0.f00001'
c         call read_hdr(inputname1,ftime) ! We read header to get the times
         if (nid.eq.0) write(*,*) '**FIELD,Time',sfn,ftime

         call load_field(inputname1)
c         call full_restart(inputname1,1)

!     Store quantities to interpolate
         do j=1,ntot
            struct(j,1) = vx(j,1,1,1)
            struct(j,2) = vz(j,1,1,1)
            struct(j,3) = T(j,1,1,1,1)
            struct(j,4) = vy(j,1,1,1)
         enddo

!     Do necessary checks before interpolation
         if(icalld.eq.0) then
            if(nid.eq.0) then
!     Compare number of points in interpolating mesh npoints
!     with maximum defined in SIZE file lhis
               if(nnhis.gt.lhis) then
                  write(6,*) 'Increase lhis to npoints',lhis,nnhis
                  call exitt
               endif
               write(6,*) 'found ', nnhis, ' points in interp. mesh'
            endif
!     Setup for interpolation tool. Use default tolerance of -1
            call intpts_setup(-1.0,inth_hpts)
        endif
      
!     Start interpolation
         if(icalld.eq.0) then
            
!     Conceptually, locate npt points. The data corresponding to each point
!     is whether it is inside an element, closesst to a border, not found.
!     Also identify the processor where the point was found, the element
!     where the point was found, parametric coordinates of the point and
!     the distance squared from found to sought point
      
            call fgslib_findpts(inth_hpts,rcode,1,
     &           proc,1,
     &           elid,1,
     &           rst,ndim,
     &           dist,1,
     &           pts(1,1),ndim,
     &           pts(2,1),ndim,
     &           pts(3,1),ndim,nnhis)
          
!     Check the return code
            do i=1,nnhis
!     Interpolating point is on boundary or outside the SEM mesh
               if(rcode(i).eq.1) then
                  if(dist(i).gt.1e-12) then
                     write(6,'(A,4E15.7)')
     &                    ' WARNING: point on boundary or outside
     &the mesh xy[z]d^2:',
     &                    (pts(k,i),k=1,ndim),dist(i)
                  endif
!     Interpolating ponit is not in the SEM mesh
               elseif(rcode(i).eq.2) then
                  nfail = nfail + 1
c               write(6,'(A,3E15.7)')
c     &              ' WARNING: point not within mesh xy[z]: !',
c     &              (pts(k,i),k=1,ndim)
               endif
            enddo
            icalld = 1
         endif                  ! icalld
       
         do i_field=1,n_out_fields
            call fgslib_findpts_eval(inth_hpts,
     &           buff_fieldout(1,i_field),1,
     &           rcode,1,
     &           proc,1,
     &           elid,1,
     &           rst,ndim,nnhis,
     &           struct(1,i_field)) ! write NEK calculated fields
         enddo                  ! each field to be interpolated
         
ccccccccccccccccccccccccccccccccccccccccc
         if (nid.eq.0) then
!     Open file to store data
            write(pippo,'(i2.2)') sfn
            inputnameu = 'ZSTAT/U'//trim(pippo)
            open(unit=37,form='unformatted',file=inputnameu)
!     Parameters to write on the header
            write(val1,'(1p15e17.9)') ftime ! Time
            write(val2,'(9i9)')      sfn ! Field number
!     Write header
            write(37) '(Time ='//trim(val1)
     &           //') (Field ='//trim(val2)
     &           //')'
!     Write values corresponding to the header
            write(37) ftime,sfn
            

!     Open file to store data
            write(pippo,'(i2.2)') sfn
            inputnamew = 'ZSTAT/W'//trim(pippo)
            open(unit=38,form='unformatted',file=inputnamew)
!     Parameters to write on the header
            write(val1,'(1p15e17.9)') ftime ! Time
            write(val2,'(9i9)')      sfn ! Field number
!     Write header
            write(38) '(Time ='//trim(val1)
     &           //') (Field ='//trim(val2)
     &           //')'
!     Write values corresponding to the header
            write(38) ftime,sfn
            
!     Open file to store data
            write(pippo,'(i2.2)') sfn
            inputnamel = 'ZSTAT/L'//trim(pippo)
            open(unit=39,form='unformatted',file=inputnamel)
!     Parameters to write on the header
            write(val1,'(1p15e17.9)') ftime ! Time
            write(val2,'(9i9)')      sfn ! Field number
!     Write header
            write(39) '(Time ='//trim(val1)
     &           //') (Field ='//trim(val2)
     &           //')'
!     Write values corresponding to the header
            write(39) ftime,sfn

!     Open file to store data
            write(pippo,'(i2.2)') sfn
            inputnamev = 'ZSTAT/V'//trim(pippo)
            open(unit=40,form='unformatted',file=inputnamev)
!     Parameters to write on the header
            write(val1,'(1p15e17.9)') ftime ! Time
            write(val2,'(9i9)')      sfn ! Field number
!     Write header
            write(40) '(Time ='//trim(val1)
     &           //') (Field ='//trim(val2)
     &           //')'
!     Write values corresponding to the header
            write(40) ftime,sfn
            
            
!     History file with number of points per core
            if(sfn.eq.1) then
               open(unit=33,access='append',file='history2.txt')
            endif
            
         endif ! Opening files

!     Store interpolated fields in arrays for communication
         do i_buff=1,nnhis
            b_u(i_buff)=buff_fieldout(i_buff,1)
            b_w(i_buff)=buff_fieldout(i_buff,2)
            b_l(i_buff)=buff_fieldout(i_buff,3)
            b_v(i_buff)=buff_fieldout(i_buff,4)
         enddo

! Communication among cores
         if (nid.ne.0) then
            dest=0
            tag0=0
            tag1=1
            tag2=2
            tag3=3
            tag4=4
            call mpi_send(nnhis,1,nekreal,
     &           dest,tag0,nekcomm,ierr)
            call mpi_send(b_u,nnhis,nekreal,
     &           dest,tag1,nekcomm,ierr)
            call mpi_send(b_w,nnhis,nekreal,
     &           dest,tag2,nekcomm,ierr)
            call mpi_send(b_l,nnhis,nekreal,
     &           dest,tag3,nekcomm,ierr)
            call mpi_send(b_v,nnhis,nekreal,
     &           dest,tag4,nekcomm,ierr)
         else
            
            write(37) (b_u(i),i=1,nnhis)
            write(38) (b_w(i),i=1,nnhis)
            write(39) (b_l(i),i=1,nnhis)
            write(40) (b_v(i),i=1,nnhis)

            if(sfn.eq.1) then
               write(33,'(9i9)') nid,nnhis
            endif
                        
            do rnk=1,lp-1
               tag0=0
               tag1=1
               tag2=2
               tag3=3
               tag4=4
               call mpi_recv(bnhis,1,nekreal,
     &              rnk,tag0,nekcomm,status,ierr)
               call mpi_recv(bb_u,lhis,nekreal,
     &              rnk,tag1,nekcomm,status,ierr)
               call mpi_recv(bb_w,lhis,nekreal,
     &              rnk,tag2,nekcomm,status,ierr)
               call mpi_recv(bb_l,lhis,nekreal,
     &              rnk,tag3,nekcomm,status,ierr)
               call mpi_recv(bb_v,lhis,nekreal,
     &              rnk,tag4,nekcomm,status,ierr)

               write(37) (bb_u(i),i=1,bnhis)
               write(38) (bb_w(i),i=1,bnhis)
               write(39) (bb_l(i),i=1,bnhis)
               write(40) (bb_v(i),i=1,bnhis)
               
               if(sfn.eq.1) then
                  write(33,'(9i9)') rnk,bnhis
               endif
                              
            enddo   
         endif

         close(33)
         close(37)
         close(38)
         close(39)
         close(40)

         
      enddo                     ! Loop over the fields

      if (nid.eq.0) write(*,*) 'AND WE ARE DONE!!'
      
      call exitt
      
      return
      end

c-----------------------------------------------------------------------
      
      subroutine load_field(field)

      implicit none

      include 'SIZE'
      include 'STATS'           ! Statistics specific variables
      include 'INPUT'           ! if3d
      include 'SOLN'
      include 'TSTEP'

      character*132 field
      
      call load_fld(field)
      
c      call outpost(vx,vy,vz,pr,t,'new')
      
      return
      end

c-----------------------------------------------------------------------
      subroutine read_hdr(field,mtimee)

      implicit none

      include 'SIZE'
      include 'STATS'           ! Statistics specific variables
      include 'INPUT'           ! if3d
      include 'SOLN'
      include 'TSTEP'

      character*132 field,hdr,fmt1
      character*10 tmpchar
      integer twdsize,mnelx,mnely,mnelz,nelo,isteps,fid0,nfileoo
      integer stat_gnum
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

      subroutine usrdat2
cc      include 'SIZE'
cc      include 'TOTAL'


cc      r0 = 1.
cc      rm = -r0
cc      call rescale_x(xm1,rm,r0) ! Rescale incoming pipe radius to be 1.0
cc      call rescale_x(ym1,rm,r0)

cc      z0 = 0.
cc      z1 = 1.
cc      call rescale_x(zm1,z0,z1) ! Make certain z_orig is on [0,1]


cc      n = nx1*ny1*nz1*nelv


c    rotating x into axial position
cc      do i=1,n
cc         x_original = xm1(i,1,1,1)
cc         y_original = ym1(i,1,1,1)
cc         z_original = zm1(i,1,1,1)

cc         xm1(i,1,1,1) =  z_original
cc         ym1(i,1,1,1) =  y_original
cc         zm1(i,1,1,1) = -x_original
cc      enddo

cc      phi_1 = 6.281407035
cc      do i=1,n
cc         prmtrc_t=xm1(i,1,1,1)*phi_1
cc         prmtrc_r=ym1(i,1,1,1)+ (1.0 / 0.3)
cc         xm1(i,1,1,1) = prmtrc_r*sin(prmtrc_t)
cc         ym1(i,1,1,1) = prmtrc_r*cos(prmtrc_t)
cc      enddo

cc      param(59) = 1.
cc      param(66) = 6.
cc      param(67) = 6.

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
