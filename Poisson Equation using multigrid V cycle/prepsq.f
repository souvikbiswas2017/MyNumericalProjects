C----- THIS IS A TEST DATA PREPARATION FILE FOR POISSON PROBLEM IN 
C----- IN SQUARE DOMAIN.

C----- AUTHOR : SOUVIK BISWAS
C----- DATE : 04/10/05

C-----------------------------------------------------------------------C
      program prepsq
C-----------------------------------------------------------------------C
      implicit none
      
      integer max_level,max_num
      parameter(max_level=8)
      parameter(max_num=2**max_level+1)

      doubleprecision length
      doubleprecision u(max_num,max_num),f(max_num,max_num)

      common/param/length

      length=1.0d0
      
      call write_guess(u,max_level,max_num)
      call write_fun(f,max_level,max_num)

      stop
      end

C-----------------------------------------------------------------------C
      subroutine write_guess(u,ilevel,n)
C-----------------------------------------------------------------------C
      implicit none

      integer ilevel,n
      integer i,j
      
      doubleprecision length
      doubleprecision u(n,n)

      do i=2,n-1
        do j=2,n-1
          u(i,j)=0.0d0
        enddo
      enddo

      do i=1,n
        u(i,1)=0.0d0
        u(i,n)=0.0d0
      enddo
      do j=1,n
        u(1,j)=0.0d0
        u(n,j)=0.0d0
      enddo

      open(8,file='.//DATA//guess.dat',status='unknown')
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
      subroutine write_fun(f,ilevel,n)
C-----------------------------------------------------------------------C
      implicit none

      integer ilevel,n
      integer i,j
      
      doubleprecision x,y,h
      doubleprecision length
      doubleprecision f(n,n)

      common/param/length

      h=length/(n-1)
      do i=1,n
        do j=1,n
          x=(i-1)*h
          y=(j-1)*h
          f(i,j)=0.0
          if(x.ge.length/4 .and. x.le.3*length/4) then
            if(y.ge.length/4 .and. y.le.3*length/4) then
              f(i,j)=1.0d0
            endif
          endif
        enddo
      enddo

      open(8,file='.//DATA//fun.dat',status='unknown')
      write(8,301)ilevel,n
      do i=1,n
        write(8,302)(f(i,j),j=1,n)
      enddo
      close(8)

  301 format(2i5)
  302 format(8e13.5)

      return
      end

