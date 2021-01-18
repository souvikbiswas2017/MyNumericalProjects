C----- THIS IS A FINITE DIFFERENCE SERIAL SOR SOLVER FOR 2D
C----- POISSON EQUATION WITH DIRICHLET BOUNDARY CONDITION. FOR
C----- SIMPLICITY UNIFORM GRID IS USED AND GRID SIZE IS SAME IN BOTH X &
C----- Y DIRECTION. 

C----- AUTHOR : SOUVIK BISWAS
C----- DATE : 04/23/05

C-----------------------------------------------------------------------C
      program sorsq
C-----------------------------------------------------------------------C
      implicit none
      
      integer max_level,max_num
      parameter(max_level=8)
      parameter(max_num=2**max_level+1)

      integer cnt,maxstep

      doubleprecision maxdiff
      real delapse,timearray(2)      

      doubleprecision u(max_num,max_num)
      doubleprecision f(max_num,max_num)
      
      write(*,*)'Running SOR'

      call read_fun(f,max_level,max_num)
      call read_guess(u,max_level,max_num)

      maxstep=500000

      delapse = etime(timearray)
      do cnt=1,maxstep
        call relaxation(u,f,max_num)
        call checkconv(u,f,max_num,maxdiff)
        write(*,300)cnt,maxdiff
      
        if(maxdiff .lt. 1e-9) goto 10
      enddo
  10  continue
      delapse = etime(timearray)

      write(*,*)delapse
      write(*,*)cnt,maxdiff

  300 format(i8,e13.5)

      call write_solution(u,max_level,max_num)

      stop
      end

C-----------------------------------------------------------------------C
      subroutine relaxation(u,f,n)
C-----------------------------------------------------------------------C
      implicit none

      integer n
      integer i,j

      doubleprecision h,length
      doubleprecision beta,pi

      doubleprecision u(n,n)
      doubleprecision f(n,n)

      length=1.0d0
      h=length/n
      pi=4*atan(1.0d0)
      beta=2.0d0/(1.0d0+sin(pi/n))

      do i=2,n-1
        do j=2,n-1
          u(i,j)=(1-beta)*u(i,j)
     1           +beta*((u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)
     1           -f(i,j)*h**2)/4.0d0)
        enddo
      enddo

      return
      end

C-----------------------------------------------------------------------C
      subroutine read_guess(u,ilevel,n)
C-----------------------------------------------------------------------C
      implicit none

      integer ilevel,n
      integer rlevel,rn
      integer i,j
      
      doubleprecision u(n,n)

      open(8,file='.//DATA//guess.dat',status='unknown')

      read(8,301)rlevel,rn
      if(rlevel.ne.ilevel .or. rn.ne.n) then
        write(*,*)'WARNING:guess file is wrong'
        stop
      endif

      do i=1,n
        read(8,302)(u(i,j),j=1,n)
      enddo

      close(8)

  301 format(2i5)
  302 format(8e13.5)

      return
      end

C-----------------------------------------------------------------------C
      subroutine read_fun(f,ilevel,n)
C-----------------------------------------------------------------------C
      implicit none

      integer ilevel,n
      integer rlevel,rn
      integer i,j
      
      doubleprecision f(n,n)

      open(8,file='.//DATA//fun.dat',status='unknown')

      read(8,301)rlevel,rn
      if(rlevel.ne.ilevel .or. rn.ne.n) then
        write(*,*)'WARNING:function file is wrong'
        stop
      endif

      do i=1,n
        read(8,302)(f(i,j),j=1,n)
      enddo

      close(8)

  301 format(2i5)
  302 format(8e13.5)

      return
      end

C-----------------------------------------------------------------------C
      subroutine write_solution(u,ilevel,n)
C-----------------------------------------------------------------------C
      implicit none

      integer ilevel,n
      integer i,j
      
      doubleprecision u(n,n)

      open(8,file='.//DATA//sor_sol.dat',status='unknown')
      write(8,301)ilevel,n
      do i=1,n
        write(8,302)(u(i,j),j=1,n)
      enddo
      close(8)

  301 format(2i5)
  302 format(8e13.5)


      return
      end

C-----------------------------------------------------------------------C
      subroutine checkconv(u,f,n,maxdiff)
C-----------------------------------------------------------------------C
      implicit none

      integer n
      integer i,j
      
      doubleprecision maxdiff,norm
      doubleprecision h,length

      doubleprecision r(n,n)
      doubleprecision u(n,n)
      doubleprecision f(n,n)

      length=1.0d0
      h=length/n

      do i=1,n
        do j=1,n
          r(i,j)=0.0d0
        enddo
      enddo
      norm=0.0d0
      do i=2,n-1
        do j=2,n-1
          r(i,j)=f(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))/
     1                  h**2        
          norm=norm+r(i,j)**2
        enddo
      enddo
      
      maxdiff=sqrt(norm)/n/n

      return
      end

