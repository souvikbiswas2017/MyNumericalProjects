C----- THIS IS A FINITE DIFFERENCE SERIAL MULTIGRID SOLVER FOR 2D
C----- POISSON EQUATION WITH DIRICHLET BOUNDARY CONDITION. FOR
C----- SIMPLICITY UNIFORM GRID IS USED AND GRID SIZE IS SAME IN BOTH X &
C----- Y DIRECTION. THE MAIN STORAGE IS WORK ARRAY WHICH IS ONLY
C----- ONE-DIMENSIONAL. LENGTH OF THIS ARRAY IS DETERMINED FROM THE
C----- NUMBER OF MULTIGRID LEVELS. FOR SIMPLICITY IN THIS PROGRAM WE ARE
C----- RESTRICTING THE NUMBER OF INTERVALS IN EACH DIRECTION AS A POWER
C----- OF 2 ONLY( i.e. 2^LVL ). THUS SPECIFYING THE GRIDSIZE IS NOT
C----- NECESSARY. IT CAN BE CALULATED FROM MAX. NO. OF LEVELS
C----- COPYVAL AND STOREVAL ROUTINES ARE WRITTEN TO ACCESS AND STORE
C----- VALUE FROM THE WORK ARRAY. THE WORK ARRAY ONLY CONTAINS THE
C----- INTERIOR POINTS AND STORE THE VALUE ROWWISE( i.e. IN BLOCKS OF
C----- ROW.)

C----- AUTHOR : SOUVIK BISWAS
C----- DATE : 04/10/05

C-----------------------------------------------------------------------C
      program mvsq
C-----------------------------------------------------------------------C
      implicit none
      
      integer max_level,memlen
      parameter(max_level=8)
      parameter(memlen=4*((4**max_level-1)/3+2**max_level-1)+max_level)

      integer ilevel
      integer cnt,maxstep

      doubleprecision length
      doubleprecision maxdiff
      real delapse,timearray(2)

C---- Assign global work arrays. These arrays are 1-D. Data form
C---- these arrays can be read into a 2-D temporary array by copyval
C---- subroutine. Data can be written into these arrays from a 2-D
C---- temporary array by storeval subroutine

      doubleprecision worku(memlen),workf(memlen)
      
      external mv

      common/param/length
      
      write(*,*)'Running MV cycle'

      length=1.0d0
      ilevel=max_level
      
      call init_1darray(worku,memlen)
      call init_1darray(workf,memlen)

      call read_fun(workf,max_level,memlen)
      call read_guess(worku,max_level,memlen)

      maxstep=500

      delapse = etime(timearray)
      do cnt=1,maxstep
        call cleanwork(workf,worku,max_level,memlen)
        call mv(ilevel,worku,workf,memlen,mv)
        call checkconv(worku,workf,max_level,memlen,maxdiff)
        write(*,300)cnt,maxdiff

        if(maxdiff .lt. 1e-9) goto 10
      enddo
  10  continue
      delapse = etime(timearray)

      write(*,*)delapse
      write(*,*)cnt,maxdiff

  300 format(i8,e13.5)

      call write_solution(worku,max_level,memlen)

      stop
      end

C-----------------------------------------------------------------------C
      subroutine mv(ilevel,worku,workf,memlen,dummv)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,iprevlevel
      integer n1,n2

      doubleprecision worku(memlen),workf(memlen)

      external dummv

      iprevlevel=ilevel-1
      n1=2
      n2=2

      call relaxation(worku,workf,ilevel,n1,memlen)

      if(ilevel.gt.1) then
        call fg_restriction(worku,workf,ilevel,iprevlevel,memlen)
        call dummv(iprevlevel,worku,workf,memlen,dummv)
        call cg_correction(worku,iprevlevel,ilevel,memlen)
      endif

      call relaxation(worku,workf,ilevel,n2,memlen)
        
      return
      end

C-----------------------------------------------------------------------C
      subroutine relaxation(worku,workf,ilevel,ncnt,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n
      integer ncnt
      integer i,j,k

      doubleprecision h,length

      doubleprecision worku(memlen),workf(memlen)


      doubleprecision u(2**ilevel+1,2**ilevel+1)
      doubleprecision f(2**ilevel+1,2**ilevel+1)

      common/param/length

      n=2**ilevel+1
      h=length/(n-1)

      call copyval(u,worku,ilevel,memlen)
      call copyval(f,workf,ilevel,memlen)

      do k=1,ncnt
        do i=2,n-1
          do j=2,n-1
            u(i,j)=(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)
     1             -f(i,j)*h**2)/4.0d0
          enddo
        enddo
      enddo

      call storeval(u,worku,ilevel,memlen)

      return
      end

C-----------------------------------------------------------------------C
      subroutine fg_restriction(worku,workf,ilevel,iprevlevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,iprevlevel,n,nprev
      integer i,j

      doubleprecision worku(memlen),workf(memlen)

      doubleprecision ufine(2**ilevel+1,2**ilevel+1)
      doubleprecision ffine(2**ilevel+1,2**ilevel+1)
      doubleprecision r(2**ilevel+1,2**ilevel+1)
      doubleprecision fcoarse(2**iprevlevel+1,2**iprevlevel+1)
      
      n=2**ilevel+1
      nprev=2**iprevlevel+1

      call init_2darray(r,n)
      call init_2darray(fcoarse,nprev)
      call copyval(ufine,worku,ilevel,memlen)
      call copyval(ffine,workf,ilevel,memlen)

      call residual(r,ufine,ffine,ilevel)

      do i=2,(n-1)/2
        do j=2,(n-1)/2
          fcoarse(i,j)=0.0625d0*(
     1      r(2*i-2,2*j-2)+r(2*i-2,2*j)+r(2*i,2*j-2)+r(2*i,2*j)+
     1   2*(r(2*i-2,2*j-1)+r(2*i-1,2*j-2)+r(2*i,2*j-1)+r(2*i-1,2*j))+
     1   4*r(2*i-1,2*j-1)
     1                        )
        enddo
      enddo

      call storeval(fcoarse,workf,iprevlevel,memlen)

      return
      end
      
C-----------------------------------------------------------------------C
      subroutine cg_correction(worku,ilevel,inextlevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,inextlevel,n,nnext
      integer i,j

      doubleprecision worku(memlen)

      doubleprecision ucoarse(2**ilevel+1,2**ilevel+1)
      doubleprecision ufine(2**inextlevel+1,2**inextlevel+1)
      doubleprecision corrfine(2**inextlevel+1,2**inextlevel+1)
      
      n=2**ilevel+1
      nnext=2**inextlevel+1

      call init_2darray(corrfine,nnext)
      call copyval(ucoarse,worku,ilevel,memlen)
      call copyval(ufine,worku,inextlevel,memlen)

      do i=2,n-1
        do j=2,n-1
          corrfine(2*i-1,2*j-1)=ucoarse(i,j)
        enddo
      enddo

      do i=1,n-1
        do j=2,n-1
          corrfine(2*i,2*j-1)=0.5d0*(ucoarse(i,j)+ucoarse(i+1,j))
        enddo
      enddo

      do i=2,n-1
        do j=1,n-1
          corrfine(2*i-1,2*j)=0.5d0*(ucoarse(i,j)+ucoarse(i,j+1))
        enddo
      enddo

      do i=1,n-1
        do j=1,n-1
          corrfine(2*i,2*j)=0.25d0*(ucoarse(i,j)+ucoarse(i,j+1)+
     1                         ucoarse(i+1,j)+ucoarse(i+1,j+1))
        enddo
      enddo

      do i=1,nnext
        do j=1,nnext
          ufine(i,j)=ufine(i,j)+corrfine(i,j)
        enddo
      enddo

      call storeval(ufine,worku,inextlevel,memlen)
      
      return
      end

C-----------------------------------------------------------------------C
      subroutine copyval(a,worka,ilevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n,iprevlevel,strtindex,endindex
      integer i,j,k

      doubleprecision worka(memlen)
      doubleprecision a(2**ilevel+1,2**ilevel+1)

      n=2**ilevel+1
      iprevlevel = ilevel-1

      strtindex = 4*((4**iprevlevel-1)/3+2**iprevlevel-1)+iprevlevel+1
      endindex = 4*((4**ilevel-1)/3+2**ilevel-1)+ilevel

      k=strtindex
      do j=1,n
        do i=1,n
          a(i,j)=worka(k)
          k=k+1
        enddo
      enddo

      if (k.ne.(endindex+1)) then
        write(*,*)'Warning:k and endindex not equal',k,endindex
        stop
      endif

      return
      end

C-----------------------------------------------------------------------C
      subroutine storeval(a,worka,ilevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n,iprevlevel,strtindex,endindex
      integer i,j,k

      doubleprecision worka(memlen)
      doubleprecision a(2**ilevel+1,2**ilevel+1)

      n=2**ilevel+1
      iprevlevel = ilevel-1
      
      strtindex = 4*((4**iprevlevel-1)/3+2**iprevlevel-1)+iprevlevel+1
      endindex = 4*((4**ilevel-1)/3+2**ilevel-1)+ilevel

      k=strtindex
      do j=1,n
        do i=1,n
          worka(k)=a(i,j)
          k=k+1
        enddo
      enddo

      if (k.ne.(endindex+1)) then
        write(*,*)'Warning:k and endindex not equal',k,endindex
        stop
      endif

      return
      end

C-----------------------------------------------------------------------C
      subroutine residual(r,u,f,ilevel)
C-----------------------------------------------------------------------C
      implicit none

      integer ilevel,n
      integer i,j

      doubleprecision h,length
      
      doubleprecision r(2**ilevel+1,2**ilevel+1)
      doubleprecision u(2**ilevel+1,2**ilevel+1)
      doubleprecision f(2**ilevel+1,2**ilevel+1)

      common/param/length

      n=2**ilevel+1
      h=length/(n-1)

      do i=2,n-1
        do j=2,n-1
          r(i,j)=f(i,j)-(u(i+1,j)+u(i-1,j)+u(i,j+1)+u(i,j-1)-4*u(i,j))
     1                  /h**2
        enddo
      enddo

      return
      end

C-----------------------------------------------------------------------C
      subroutine read_guess(worku,ilevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n
      integer rlevel,rn
      integer i,j
      
      doubleprecision worku(memlen)
      doubleprecision u(2**ilevel+1,2**ilevel+1)

      n=2**ilevel+1

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

      call storeval(u,worku,ilevel,memlen)

      return
      end

C-----------------------------------------------------------------------C
      subroutine read_fun(workf,ilevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n
      integer rlevel,rn
      integer i,j
      
      doubleprecision workf(memlen)
      doubleprecision f(2**ilevel+1,2**ilevel+1)

      n=2**ilevel+1

      open(8,file='.//DATA//fun.dat',status='unknown')

      read(8,301)rlevel,rn
      if(rlevel.ne.ilevel .or. rn.ne.n) then
        write(*,*)'WARNING:function file is wrong'
      endif

      do i=1,n
        read(8,302)(f(i,j),j=1,n)
      enddo

      close(8)

  301 format(2i5)
  302 format(8e13.5)

      call storeval(f,workf,ilevel,memlen)

      return
      end

C-----------------------------------------------------------------------C
      subroutine write_solution(worku,ilevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n
      integer i,j
      
      doubleprecision worku(memlen)
      doubleprecision u(2**ilevel+1,2**ilevel+1)

      n=2**ilevel+1

      call copyval(u,worku,ilevel,memlen)

      open(8,file='.//DATA//mg_sol.dat',status='unknown')
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
      subroutine cleanwork(workf,worku,ilevel,memlen)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer iprevlevel,ilevel
      integer i,endindex
      
      doubleprecision worku(memlen)
      doubleprecision workf(memlen)
      
      iprevlevel = ilevel-1
      endindex = 4*((4**iprevlevel-1)/3+2**iprevlevel-1)+iprevlevel

      do i=1,endindex
        worku(i)=0.0d0
        workf(i)=0.0d0
      enddo

      return
      end

C-----------------------------------------------------------------------C
      subroutine checkconv(worku,workf,ilevel,memlen,maxdiff)
C-----------------------------------------------------------------------C
      implicit none

      integer memlen
      integer ilevel,n
      integer i,j

      doubleprecision maxdiff,norm
      doubleprecision h,length

      doubleprecision worku(memlen),workf(memlen)

      doubleprecision r(2**ilevel+1,2**ilevel+1)
      doubleprecision u(2**ilevel+1,2**ilevel+1)
      doubleprecision f(2**ilevel+1,2**ilevel+1)

      common/param/length

      n=2**ilevel+1
      h=length/(n-1)

      call init_2darray(r,n)
      call copyval(u,worku,ilevel,memlen)
      call copyval(f,workf,ilevel,memlen)

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

C-----------------------------------------------------------------------C
      subroutine init_2darray(a,n)
C-----------------------------------------------------------------------C
      implicit none

      integer n
      integer i,j
      doubleprecision a(n,n)

      do i=1,n
        do j=1,n
          a(i,j)=0.0d0
        enddo
      enddo

      return
      end

C-----------------------------------------------------------------------C
      subroutine init_1darray(a,n)
C-----------------------------------------------------------------------C
      implicit none

      integer n
      integer i
      doubleprecision a(n)

      do i=1,n
        a(i)=0.0d0
      enddo

      return
      end

