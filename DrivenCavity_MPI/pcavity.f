C
C     NAME: DRIVEN CAVITY PROBLEM BY MAC METHOD
C     AUTHOR: SOUVIK BISWAS
C     DESCRIPTION: SOLVES TWO DIMENSIONAL DRIVEN CAVITY PROBLEM IN
C                  RECTANGULAR CO-ORDINATE SYSTEM USING FINITE VOLUME
C                  APPROACH (MAC) AND SOR IETRATION FOR PRESSURE
C                  SOLUTION
C     DATE:07/03/04 11:35
C
      program main
C-----------------------------------------------------------------------C
C     Main Program
C-----------------------------------------------------------------------C
      implicit none
      include "mpif.h"
C     Process Grid
      integer comm2d
      integer myid,ierr,numprocs,stride
      integer nxproc,nyproc
      integer dims(2),coords(2),neighbor(8),bd(8)
      double precision stime,etime,ptime
      logical periods(2)
C     Grid variables
      integer nxp,nyp,nxp1,nyp1,nxp2,nyp2,maxn
      integer sx,ex,sy,ey,sxv,exv,syv,eyv,k
C     SOR variables
      double precision beta,eps
C     Convergence variables
      double precision unorm,vnorm

C     Initialize Variables
      parameter (maxn=200)
      parameter (nxproc=1,nyproc=1)
c     parameter (nxproc=1,nyproc=2)
c     parameter (nxproc=2,nyproc=2)
c     parameter (nxproc=2,nyproc=4)
C     parameter (nxproc=4,nyproc=4)
      parameter (beta=1.5d0)

C     Problem parameters
      double precision xlength,ylength,dt,h,time
      double precision uwall,rho,meu,endtime
      integer iPrint
C     Defining the Arrays
      double precision P(maxn,maxn),p0(maxn,maxn),a(maxn,maxn),
     1          u(maxn,maxn),v(maxn,maxn),ut(maxn,maxn),vt(maxn,maxn)
      double precision uu(maxn,maxn),vv(maxn,maxn),w(maxn,maxn),
     1                 pp(maxn,maxn)

C     Initialize MPI,get myrank, get mysize
      call MPI_INIT(ierr)
      call MPI_COMM_RANK(MPI_COMM_WORLD,myid,ierr)
      call MPI_COMM_SIZE(MPI_COMM_WORLD,numprocs,ierr)

C     Read the input data and broadcast it
      if (myid.eq.0) then
        call pgreadin(xlength,ylength,uwall,rho,meu,iPrint,endtime,
     1                nxp,nyp,nxp1,nyp1,nxp2,nyp2)
      endif

      call pgbcast(xlength,ylength,uwall,rho,meu,iPrint,endtime,
     1             nxp,nyp,nxp1,nyp1,nxp2,nyp2)
C     Initialize
      eps=1.0d0
      h=xlength/nxp
      time=0.0d0

C     Call routine to give processgrid,neighbors,domain decomposition
      call DECOMPOSTION(myid,numprocs,comm2d,dims,coords,periods,bd,
     1                  neighbor,nxproc,nyproc,nxp2,nyp2,
     1                  sx,ex,sy,ey,sxv,exv,syv,eyv)
      
C     Define mpi noncontigous datatype
      call MPI_TYPE_VECTOR(ey-sy+1,1,ex-sx+3,MPI_DOUBLE_PRECISION,
     1                     stride,ierr)
      call MPI_TYPE_COMMIT(stride,ierr)


C     Start Calculation
C     Initialize the local matrices
      call pginit(P,p0,a,u,v,ut,vt,bd,sx,ex,sy,ey)
      call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      stime=MPI_WTIME()
C     Time Loop
      do while (time.le.endtime)
C       Set Timestep Increment
        call setdt(comm2d,meu,rho,u,v,h,dt,sx,ex,sy,ey)
c       if(myid.eq.0) write(*,*)'time = ',time,' dt = ',dt

C       Velocity Boundary Conditions
        call pgbdry(u,v,bd,sx,ex,sy,ey,sxv,exv,syv,eyv,uwall)

C       Advect the Velocity
        call pgadvect(comm2d,neighbor,stride,u,v,ut,vt,meu,h,dt,
     1                sx,ex,sy,ey,sxv,exv,syv,eyv)

C       Solving the Pressure Equation
        call pgpressure(myid,comm2d,neighbor,stride,P,p0,a,ut,vt,
     1                  beta,eps,h,dt,
     1                  sx,ex,sy,ey,sxv,exv,syv,eyv)

C       Updating the Velocity
        call pgcorrect(comm2d,neighbor,stride,P,u,v,ut,vt,h,dt,
     1                 sx,ex,sy,ey,sxv,exv,syv,eyv)

        time=time+dt
      enddo
C     End of Calculation

      etime=MPI_WTIME()
      ptime=etime-stime

C     Calculating vorticity at the pressure points
      call pgpostproc(u,v,uu,vv,w,h,
     1                sx,ex,sy,ey,sxv,exv,syv,eyv)

C     Print Data
      if (iPrint.eq.1) then
        call pgprint(bd,P,h,time,'pres',coords,
     1               sx,ex,sy,ey)
        call pgprint(bd,uu,h,time,'uvel',coords,
     1               sx,ex,sy,ey)
        call pgprint(bd,vv,h,time,'vvel',coords,
     1               sx,ex,sy,ey)
        call pgprint(bd,w,h,time,'vort',coords,
     1               sx,ex,sy,ey)
      endif

      write(*,*) 'Process Time = ',ptime
      write(*,*) 'Done'

      call MPI_COMM_FREE( comm2d, ierr )
      call MPI_FINALIZE(ierr)

      end
C-----End-of-main-------------------------------------------------------C

      subroutine pgreadin(xlength,ylength,uwall,rho,meu,iPrint,endtime,
     1                   nxp,nyp,nxp1,nyp1,nxp2,nyp2)
C-----------------------------------------------------------------------C
C     Read input data(Done by root node only)
C-----------------------------------------------------------------------C
      implicit none

      double precision xlength,ylength
      double precision uwall,rho,meu,endtime
      integer nxp,nyp,nxp1,nyp1,nxp2,nyp2
      integer iPrint

      write(*,*)'Give endtime'
      read(*,*)endtime
      write(*,*)'endtime = ',endtime

      write(*,*)'Give x_length,y_length'
      read(*,*)xlength,ylength
      write(*,*)'x_length,y_length = ',xlength,ylength

      write(*,*)'Give Uwall,Density,Viscosity'
      read(*,*)uwall,rho,meu
      write(*,*)'Uwall,Density,Viscosity = ',uwall,rho,meu

      write(*,*)'Give nxp,nyp'
      read(*,*)nxp,nyp
      write(*,*)'nxp,nyp = ',nxp,nyp

      write(*,*)'Give iPrint'
      read(*,*)iPrint
      write(*,*)'iPrint = ',iPrint

      nxp1=nxp+1
      nyp1=nyp+1
      nxp2=nxp+2
      nyp2=nyp+2

      return
      end
C-----End-of-pgreadin---------------------------------------------------C

      subroutine pgbcast(xlength,ylength,uwall,rho,meu,iPrint,endtime,
     1                   nxp,nyp,nxp1,nyp1,nxp2,nyp2)
C-----------------------------------------------------------------------C
C     Broad cast input data read by root node
C-----------------------------------------------------------------------C
      implicit none
      include "mpif.h"

      double precision xlength,ylength
      double precision uwall,rho,meu,endtime
      integer nxp,nyp,nxp1,nyp1,nxp2,nyp2
      integer iPrint

      integer ierr

      call MPI_BCAST(xlength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     1               ierr)
      call MPI_BCAST(ylength,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     1               ierr)
      call MPI_BCAST(uwall,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rho,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(meu,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(endtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,
     1               ierr)

      call MPI_BCAST(nxp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nyp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nxp1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nyp1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nxp2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nyp2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iPrint,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      return
      end
C-----End-of-pgbcast----------------------------------------------------C

      subroutine DECOMPOSTION(myid,numprocs,comm2d,dims,coords,periods,
     1                        bd,neighbor,nxproc,nyproc,nxp2,nyp2,
     1                        sx,ex,sy,ey,sxv,exv,syv,eyv)
C-----------------------------------------------------------------------C
C     Routine to give processgrid,neighbors,domain decomposition
C-----------------------------------------------------------------------C
      implicit none
      include "mpif.h"

      integer myid,numprocs,comm2d,ierr
      integer nxproc,nyproc,nxp2,nyp2
      integer dims(2),coords(2),neighbor(8),bd(8)
      logical periods(2)

      integer sx,ex,sy,ey
      integer sxv,exv,syv,eyv

      integer ngb(3)
      integer i,j

      if(nxproc*nyproc.ne.numprocs) then
        write(*,*)'Aborting program: incorrect no. of processors'
        stop
      endif
      dims(1)=nxproc
      dims(2)=nyproc
      periods(1)=.false.
      periods(2)=.false.


C     Create a division of processors in a cartesian grid
      call MPI_DIMS_CREATE(numprocs,2,dims,ierr)

C     Create a new communicator with topology information
      call MPI_CART_CREATE(MPI_COMM_WORLD,2,dims,periods,.true.,
     1                     comm2d,ierr)

C     Get my position in this communicator
      call MPI_COMM_RANK(comm2d,myid,ierr)

C     Get neighbours
      do i=1,8
        bd(i)=0
      enddo
C     Get my co-ordinate in this communicator
      call MPI_CART_GET(comm2d,2,dims,periods,coords,ierr)
      
      if(coords(1).eq.0) then
        bd(4)=1
        bd(5)=1
        bd(6)=1
      endif
      if(coords(1).eq.dims(1)-1) then
        bd(8)=1
        bd(1)=1
        bd(2)=1
      endif
      if(coords(2).eq.0) then
        bd(2)=1
        bd(3)=1
        bd(4)=1
      endif
      if(coords(2).eq.dims(2)-1) then
        bd(6)=1
        bd(7)=1
        bd(8)=1
      endif

C     Line boundaries
      if(bd(1).eq.0) then
        ngb(1)=coords(1)+1
        ngb(2)=coords(2)
        call MPI_CART_RANK(comm2d,ngb,neighbor(1),ierr)
      else
        neighbor(1)=MPI_PROC_NULL
      endif
      if(bd(3).eq.0) then
        ngb(1)=coords(1)
        ngb(2)=coords(2)-1
        call MPI_CART_RANK(comm2d,ngb,neighbor(3),ierr)
      else
        neighbor(3)=MPI_PROC_NULL
      endif
      if(bd(5).eq.0) then
        ngb(1)=coords(1)-1
        ngb(2)=coords(2)
        call MPI_CART_RANK(comm2d,ngb,neighbor(5),ierr)
      else
        neighbor(5)=MPI_PROC_NULL
      endif
      if(bd(7).eq.0) then
        ngb(1)=coords(1)
        ngb(2)=coords(2)+1
        call MPI_CART_RANK(comm2d,ngb,neighbor(7),ierr)
      else
        neighbor(7)=MPI_PROC_NULL
      endif

C     Corner boundaries
      if(bd(2).eq.0) then
        ngb(1)=coords(1)+1
        ngb(2)=coords(2)-1
        call MPI_CART_RANK(comm2d,ngb,neighbor(2),ierr)
      else
        neighbor(2)=MPI_PROC_NULL
      endif
      if(bd(4).eq.0) then
        ngb(1)=coords(1)-1
        ngb(2)=coords(2)-1
        call MPI_CART_RANK(comm2d,ngb,neighbor(4),ierr)
      else
        neighbor(4)=MPI_PROC_NULL
      endif
      if(bd(6).eq.0) then
        ngb(1)=coords(1)-1
        ngb(2)=coords(2)+1
        call MPI_CART_RANK(comm2d,ngb,neighbor(6),ierr)
      else
        neighbor(6)=MPI_PROC_NULL
      endif
      if(bd(8).eq.0) then
        ngb(1)=coords(1)+1
        ngb(2)=coords(2)+1
        call MPI_CART_RANK(comm2d,ngb,neighbor(8),ierr)
      else
        neighbor(8)=MPI_PROC_NULL
      endif

C     Domain decomposition
      call MPE_DECOMP1D(nxp2-2,dims(1),coords(1),sx,ex)
      call MPE_DECOMP1D(nyp2-2,dims(2),coords(2),sy,ey)
      sx=sx+1
      sxv=sx
      ex=ex+1
      exv=ex
      if(bd(1).eq.1) exv=ex-1
      sy=sy+1
      syv=sy
      ey=ey+1
      eyv=ey
      if(bd(7).eq.1) eyv=ey-1

      return
      end
C-----End-of-DECOMPOSITI0N----------------------------------------------C

      subroutine MPE_DECOMP1D(n,numprocs,myid,s,e)
C-----------------------------------------------------------------------C
C     Routine to compute start and end indices of each processor
C-----------------------------------------------------------------------C
      implicit none
      integer n,numprocs,myid,s,e
      integer nlocal
      integer deficit
      nlocal=n/numprocs
      s=myid*nlocal + 1
      deficit=mod(n,numprocs)
      s=s+min(myid,deficit)
      if(myid .lt. deficit) then
        nlocal = nlocal + 1
      endif
      e=s+nlocal-1
      if(e.gt.n.or.myid.eq.numprocs-1) e=n

      return
      end
C-----End-of-MPE_DECOMP1D-----------------------------------------------C

      subroutine pginit(P,p0,a,u,v,ut,vt,bd,sx,ex,sy,ey,sxv,exv,syv,eyv)
C-----------------------------------------------------------------------C
C     Initialize arrays
C-----------------------------------------------------------------------C
      implicit none

      integer sx,ex,sy,ey,sxv,exv,syv,eyv
      integer bd(8)
      double precision P(sx-1:ex+1,sy-1:ey+1),p0(sx-1:ex+1,sy-1:ey+1),
     1     a(sx-1:ex+1,sy-1:ey+1),u(sx-1:ex+1,sy-1:ey+1),
     1     v(sx-1:ex+1,sy-1:ey+1),ut(sx-1:ex+1,sy-1:ey+1),
     1     vt(sx-1:ex+1,sy-1:ey+1)

      integer i,j

      do i=sx-1,ex+1
        do j=sy-1,ey+1
          P(i,j)=0.0d0
          p0(i,j)=0.0d0
          a(i,j)=1/4.0d0
          u(i,j)=0.0d0
          ut(i,j)=0.0d0
          v(i,j)=0.0d0
          vt(i,j)=0.0d0
        enddo
      enddo
      
      if(bd(1).eq.1) then
        do j=sy,ey
          a(ex,j)=1/3.0d0
        enddo
      endif
      if(bd(5).eq.1) then
        do j=sy,ey
          a(sx,j)=1/3.0d0
        enddo
      endif
      if(bd(3).eq.1) then
        do i=sx,ex
          a(i,sy)=1/3.0d0
        enddo
      endif
      if(bd(7).eq.1) then
        do i=sx,ex
          a(i,ey)=1/3.0d0
        enddo
      endif
      if(bd(1).eq.1.and.bd(3).eq.1) then
        a(ex,sy)=1/2.0d0
      endif
      if(bd(1).eq.1.and.bd(7).eq.1) then
        a(ex,ey)=1/2.0d0
      endif
      if(bd(5).eq.1.and.bd(3).eq.1) then
        a(sx,sy)=1/2.0d0
      endif
      if(bd(5).eq.1.and.bd(7).eq.1) then
        a(sx,ey)=1/2.0d0
      endif

      return
      end
C-----End-of-pginit-----------------------------------------------------C

      subroutine setdt(comm2d,meu,rho,u,v,h,dt,sx,ex,sy,ey)
C-----------------------------------------------------------------------C
C     set the time step increment
C-----------------------------------------------------------------------C
      implicit none
      include "mpif.h"

      integer comm2d
      integer sx,ex,sy,ey
      double precision meu,rho,h,dt
      double precision u(sx-1:ex+1,sy-1:ey+1),v(sx-1:ex+1,sy-1:ey+1)

      double precision dtvisc,dtadv,Velmax,Vel,gVelmax
      integer i,j,ierr

      dtvisc=0.25d0*h**2/meu

      Velmax=0.0d0
      gVelmax=0.0d0
      Vel=0.0d0
      do i=sx,ex
        do j=sy,ey
          Vel = sqrt((u(i,j)/2+u(i-1,j)/2)**2+(v(i,j)/2+v(i,j-1)/2)**2)
          if(Velmax.lt.Vel) then
            Velmax=Vel
          endif
        enddo
      enddo
      call MPI_ALLREDUCE(Velmax,gVelMax,1,MPI_DOUBLE_PRECISION,
     1                     MPI_SUM,comm2d,ierr)
      dtadv=0.5d0*2.0d0*meu/rho/gVelmax

      if (dtvisc.le.dtadv) then
        dt=dtvisc
      else
        dt =dtadv
      endif
      dt=dtvisc
      dt=0.5d0*dt

      if(dt.lt.1.0e-6) then
        write(*,*) 'dt is too small'
        write(*,*) 'stopping program'
        stop
      endif

      return
      end
C-----End-of-setdt------------------------------------------------------C

      subroutine pgbdry(u,v,bd,sx,ex,sy,ey,sxv,exv,syv,eyv,uwall)
C-----------------------------------------------------------------------C
C     Set Velocity boundary conditions
C-----------------------------------------------------------------------C
      implicit none

      integer sx,ex,sy,ey,sxv,exv,syv,eyv
      integer bd(8)
      double precision uwall
      double precision u(sx-1:ex+1,sy-1:ey+1),v(sx-1:ex+1,sy-1:ey+1)

      integer i,j

      if (bd(3).eq.1) then
        do i=sx-1,ex+1
          u(i,sy-1)=-u(i,sy)
        enddo
      endif
      if (bd(7).eq.1) then
        do i=sx-1,ex+1
          u(i,ey+1)=2*uwall-u(i,ey)
        enddo
      endif
      if (bd(5).eq.1) then
        do j=sy-1,ey+1
          v(sx-1,j)=-v(sx,j)
        enddo
      endif
      if (bd(1).eq.1) then
        do j=sy-1,ey+1
          v(ex+1,j)=-v(ex,j)
        enddo
      endif

      return
      end
C-----End-of-pgbdry-----------------------------------------------------C

      subroutine pgadvect(comm2d,neighbor,stride,u,v,ut,vt,meu,h,dt,
     1                    sx,ex,sy,ey,sxv,exv,syv,eyv)
C-----------------------------------------------------------------------C
C     Finding the projected velocity
C-----------------------------------------------------------------------C
      implicit none

      integer comm2d,neighbor(8)
      integer sx,ex,sy,ey,sxv,exv,syv,eyv,stride
      double precision dt,meu,h
      double precision u(sx-1:ex+1,sy-1:ey+1),v(sx-1:ex+1,sy-1:ey+1),
     1     ut(sx-1:ex+1,sy-1:ey+1),vt(sx-1:ex+1,sy-1:ey+1)

      integer i,j
      integer ierr

C     Finding Projected Velocity in X,
      do i=sxv,exv
        do j=sy,ey
          ut(i,j)=u(i,j)+dt*(-(0.25d0/h)*((u(i+1,j)+u(i,j))**2-(u(i,j)
     1        +u(i-1,j))**2+(u(i,j+1)+u(i,j))*(v(i+1,j)+v(i,j))-
     1        (u(i,j)+u(i,j-1))*(v(i+1,j-1)+v(i,j-1)))+(meu/h**2)*
     1        (u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j)))
        enddo
      enddo

C     Finding Projected Velocity in Y
      do i=sx,ex
        do j=syv,eyv
          vt(i,j)=v(i,j)+dt*(-(0.25d0/h)*((u(i,j+1)+u(i,j))*(v(i+1,j)+
     1        v(i,j))-(u(i-1,j+1)+u(i-1,j))*(v(i,j)+v(i-1,j))+
     1        (v(i,j+1)+v(i,j))**2-(v(i,j)+v(i,j-1))**2)+(meu/h**2)*
     1        (v(i+1,j)+v(i-1,j)+v(i,j+1)+v(i,j-1)-4*v(i,j)))
        enddo
      enddo
C     Exchange data
      call pgexchange(ut,comm2d,stride,neighbor,sx,ex,sy,ey)
      call pgexchange(vt,comm2d,stride,neighbor,sx,ex,sy,ey)
      call MPI_BARRIER(comm2d,ierr)

      return
      end
C-----End-of-pgadvect---------------------------------------------------C

      subroutine pgpressure(myid,comm2d,neighbor,stride,P,p0,a,ut,vt,
     1                      beta,eps,h,dt,
     1                      sx,ex,sy,ey,sxv,exv,syv,eyv)
C-----------------------------------------------------------------------C
C     Solves the Pressure Equation
C-----------------------------------------------------------------------C
      implicit none
      include "mpif.h"

      integer myid,comm2d,neighbor(8)
      integer sx,ex,sy,ey,sxv,exv,syv,eyv,stride
      double precision beta,eps,h,dt
      double precision P(sx-1:ex+1,sy-1:ey+1),p0(sx-1:ex+1,sy-1:ey+1),
     1     a(sx-1:ex+1,sy-1:ey+1),ut(sx-1:ex+1,sy-1:ey+1),
     1     vt(sx-1:ex+1,sy-1:ey+1)

      integer i,j,iter
      double precision epsglobal
      integer ierr

      iter=0
      epsglobal=1.0d0
      do while (epsglobal .gt. 1.0e-4)
c     do iter=1,80
        do i=sx-1,ex+1
          do j=sy-1,ey+1
            p0(i,j) = P(i,j)
          enddo
        enddo
        do i=sx,ex
          do j=sy,ey
            P(i,j)=beta*a(i,j)*(P(i+1,j)+P(i-1,j)+P(i,j+1)+P(i,j-1)-
     1         (h/dt)*(ut(i,j)-ut(i-1,j)+vt(i,j)-vt(i,j-1)))+
     1         (1-beta)*P(i,j)
          enddo
        enddo
        eps=0
        do i=sx,ex
          do j=sy,ey
            eps = eps+(P(i,j)-p0(i,j))**2
          enddo
        enddo
        eps=sqrt(eps/(ex-sx)/(ey-sy))
C       calculate global error
        call MPI_ALLREDUCE(eps,epsglobal,1,MPI_DOUBLE_PRECISION,
     1                     MPI_SUM,comm2d,ierr)
C       Exchange data
        call pgexchange(P,comm2d,stride,neighbor,sx,ex,sy,ey)
        call MPI_BARRIER(comm2d,ierr)
        iter=iter+1
      enddo
c     if(myid.eq.0) write(*,*) 'Iter = ',iter,' Error = ',epsglobal

      return
      end
C-----End-of-gpressure--------------------------------------------------C

      subroutine pgcorrect(comm2d,neighbor,stride,P,u,v,ut,vt,h,dt,
     1                     sx,ex,sy,ey,sxv,exv,syv,eyv)
C-----------------------------------------------------------------------C
C     Correct the velocities based on pressure
C-----------------------------------------------------------------------C
      implicit none

      integer comm2d,neighbor(8)
      integer sx,ex,sy,ey,sxv,exv,syv,eyv,stride
      double precision P(sx-1:ex+1,sy-1:ey+1),u(sx-1:ex+1,sy-1:ey+1),
     1     v(sx-1:ex+1,sy-1:ey+1),ut(sx-1:ex+1,sy-1:ey+1),
     1     vt(sx-1:ex+1,sy-1:ey+1)
      double precision h,dt

      integer i,j
      integer ierr
 
      do i=sxv,exv
        do j=sy,ey
          u(i,j) = ut(i,j)-(dt/h)*(P(i+1,j)-P(i,j))
        enddo
      enddo
      do i=sx,ex
        do j=syv,eyv
          v(i,j) = vt(i,j)-(dt/h)*(P(i,j+1)-P(i,j))
        enddo
      enddo
C     Exchange data
      call pgexchange(u,comm2d,stride,neighbor,sx,ex,sy,ey)
      call pgexchange(v,comm2d,stride,neighbor,sx,ex,sy,ey)
      call MPI_BARRIER(comm2d,ierr)

      return
      end
C-----End-of-pgcorrect--------------------------------------------------C

      subroutine pgpostproc(u,v,uu,vv,w,h,
     1                      sx,ex,sy,ey,sxv,exv,syv,eyv)
C-----------------------------------------------------------------------C
C     postprocess: calculate velocity,vorticity,stream function etc
C-----------------------------------------------------------------------C
      implicit none

      integer sx,ex,sy,ey,sxv,exv,syv,eyv
      double precision h
      double precision P(sx-1:ex+1,sy-1:ey+1),u(sx-1:ex+1,sy-1:ey+1),
     1          v(sx-1:ex+1,sy-1:ey+1)
      double precision uu(sx-1:ex+1,sy-1:ey+1),vv(sx-1:ex+1,sy-1:ey+1),
     1          w(sx-1:ex+1,sy-1:ey+1),pp(sx-1:ex+1,sy-1:ey+1)

      integer i,j

      do i=sxv,exv
        do j=syv,eyv
          uu(i,j)=0.5d0*(u(i,j+1)+u(i,j))
          vv(i,j)=0.5d0*(v(i+1,j)+v(i,j))
          w(i,j)=(1.0d0/2*h)*(u(i,j+1)-u(i,j)-v(i+1,j)+v(i,j))
        enddo
      enddo

      return
      end
C-----End-of-pgpostproc-------------------------------------------------C

      subroutine pgprint(bd,mat,h,time,fname,coords,
     1                   sx,ex,sy,ey)
C-----------------------------------------------------------------------C
C     Print the output data
C-----------------------------------------------------------------------C
      implicit none

      integer sx,ex,sy,ey
      integer coords(2),bd(8)
      double precision h,time
      double precision mat(sx-1:ex+1,sy-1:ey+1)
      character*6 fname

      character*1 num(10)
      integer i,j,sxp,exp,syp,eyp

      data (num(i),i=1,10)/'0','1','2','3','4','5','6','7','8','9'/

      sxp=sx
      exp=ex
      syp=sy
      eyp=ey
      if(bd(5).eq.1) sxp=sxp-1
      if(bd(1).eq.1) exp=exp+1
      if(bd(3).eq.1) syp=syp-1
      if(bd(7).eq.1) eyp=eyp+1

      fname(5:6)=num(coords(1)+1)
      fname(6:7)=num(coords(2)+1)
      open(8,file=fname,status='unknown')
      write(8,300)fname,time
      write(8,310)sx,ex,sy,ey
      write(8,320)((dble(i-1)*h),i=sx,ex)
      write(8,320)((dble(j-1)*h),j=sy,ey)
      do j=sy-1,ey+1
        write(8,320)(mat(i,j),i=sx-1,ex+1)
      enddo
300   format(a20,e12.5)
310   format(8i5)
320   format(8e13.5)
      close(8)

      return
      end
C-----End-of-pgprint----------------------------------------------------C

      subroutine pgexchange(a,comm2d,stride,neighbor,sx,ex,sy,ey)
C-----------------------------------------------------------------------C
C     Print the output data
C-----------------------------------------------------------------------C
      implicit none
      include "mpif.h"
      
      integer comm2d
      integer sx,ex,sy,ey,stride
      integer neighbor(8)
      integer status(MPI_STATUS_SIZE),nx,ierr
      double precision a(sx-1:ex+1,sy-1:ey+1)

      nx=ex-sx+1
C     SENDING LINES
C     Using Double Precision Datatype
	call MPI_SENDRECV(a(sx,ey),nx,MPI_DOUBLE_PRECISION,neighbor(7),0,
     &                  a(sx,sy-1),nx,MPI_DOUBLE_PRECISION,
     &                  neighbor(3),0,comm2d,status,ierr)
	call MPI_SENDRECV(a(sx,sy),nx,MPI_DOUBLE_PRECISION,neighbor(3),1,
     &                  a(sx,ey+1), nx, MPI_DOUBLE_PRECISION,
     &                  neighbor(7),1,comm2d,status,ierr)
C     Using strided Datatype
	call MPI_SENDRECV(a(ex,sy),1,stride,neighbor(1),0,
     &                  a(sx-1,sy),1,stride,neighbor(5),0,
     &                  comm2d,status,ierr)
	call MPI_SENDRECV(a(sx,sy),1,stride,neighbor(5),1,
     &                  a(ex+1,sy),1,stride,neighbor(1),1,
     &                  comm2d,status,ierr)

C     SENDING CORNERS
C     Using Double Precision Datatype
	call MPI_SENDRECV(a(sx,sy),1,MPI_DOUBLE_PRECISION,neighbor(8),0,
     &                  a(ex+1,ey+1),1,MPI_DOUBLE_PRECISION,
     &                  neighbor(4),0,comm2d,status,ierr)
	call MPI_SENDRECV(a(ex,ey),1,MPI_DOUBLE_PRECISION,neighbor(4),0,
     &                  a(sx-1,sy-1),1,MPI_DOUBLE_PRECISION,
     &                  neighbor(8),0,comm2d,status,ierr)
	call MPI_SENDRECV(a(sx,ey),1,MPI_DOUBLE_PRECISION,neighbor(6),0,
     &                  a(ex+1,sy-1),1,MPI_DOUBLE_PRECISION,
     &                  neighbor(6),0,comm2d,status,ierr)
	call MPI_SENDRECV(a(ex,sy),1,MPI_DOUBLE_PRECISION,neighbor(2),0,
     &                  a(sx-1,ey+1),1,MPI_DOUBLE_PRECISION,
     &                  neighbor(2),0,comm2d,status,ierr)

      return
      end
C-----End-of-pgexchline-------------------------------------------------C









