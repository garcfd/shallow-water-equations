

      program main
      implicit none

      integer, parameter :: ni=101,nj=101

      integer live(-1:ni+2,-1:nj+2)

      real xn(ni+1),yn(nj+1)
      real u(ni,nj),v(ni,nj)
      real h(0:ni+1,0:nj+1)
      real b(0:ni+1,0:nj+1)
      real Q(0:ni+1,0:nj+1,3)
      real E(0:ni+1,0:nj+1,3)
      real F(0:ni+1,0:nj+1,3)
      real Qf(3)
      real Ef(ni+1,nj,3)
      real Ff(ni,nj+1,3)
      real sb(3)

      integer i,j,n,ii,jj,it,sa,nits,nsav,gc

      real dx,dy,dz,dt,xx,dtdx,dtdy,dbdx,dbdy
      real gg,hh,uu,vv,hmax,eps

      character(31) outfile
      character(3) string



C-----specify live dead regions (if needed)
      live = 0              ! dead
      live(1:ni,  1:nj) = 1 ! live
      live(40:60,40:60) = 0 ! dead

c      do i = 1, ni
c      do j = 1, nj
c      if ((i.ge.40).and.(i.le.60).and.(j.ge.40).and.(j.le.60)) then
c        live(i,j) = 0 ! dead region
c      endif
c      enddo
c      enddo


C-----specify constants
      gg = 9.8     ! gravity
      dx = 0.1     ! spacing
      dy = 0.1     ! spacing
      dz = 0.01    ! spacing
      dt = 0.005   ! timestep
      eps = 0.01   ! smoothing (0.05)

      dtdx = dt/dx
      dtdy = dt/dy

CCC   dbdx = 0.02  ! if constant gradient
CCC   dbdy = 0.0   ! if constant gradient

      hmax = 100.0

      nsav = 100
      nits = 200

C-----create uniform mesh nodes
      xn(1) = 0.0 ! start point
      yn(1) = 0.0 ! start point

      do i = 2, ni+1
        xn(i) = xn(i-1) + dx
      enddo

      do j = 2, nj+1
        yn(j) = yn(j-1) + dy
      enddo


C-----sea bed profile
      do i = 0, ni+1
      do j = 0, nj+1

CCC     xx = real(i-1)*dx
CCC     b(i,j) = 0.005*exp(-(xx-6.0)**2) + 0.0015*xx

        b(i,j) = 0.0

      enddo
      enddo


C-----initialise arrays with height pulse
      do i = 1, ni
      do j = 1, nj

        xx = real(i-1)*dx

        h(i,j) = 0.02 - b(i,j) + 0.005*exp(-(xx-0.0)**2.0)

        u(i,j) = 0.0
        v(i,j) = 0.0

      enddo
      enddo


C-----start iterations
      do sa = 0,nsav

      if (sa.gt.0) then

      do it = 1,nits

C-----step1 - vectors at cell centres
      do i = 1, ni
      do j = 1, nj
      if (live(i,j).eq.1) then

        hh = h(i,j)
        uu = u(i,j)
        vv = v(i,j)
        
        Q(i,j,1) =    hh
        Q(i,j,2) = uu*hh
        Q(i,j,3) = vv*hh

        E(i,j,1) =    uu*hh
        E(i,j,2) = uu*uu*hh + 0.5*gg*hh*hh
        E(i,j,3) = vv*uu*hh

        F(i,j,1) =    vv*hh
        F(i,j,2) = uu*vv*hh
        F(i,j,3) = vv*vv*hh + 0.5*gg*hh*hh

      endif
      enddo
      enddo
      

C-----step2 - vectors at ghost cell centres
      do i = 0,ni+1
      do j = 0,nj+1

        gc = 0 ! boundary marker

        if (live(i  ,j).eq.0) then  ! left wall
        if (live(i+1,j).eq.1) then  ! left wall
          hh = 2*h(i+1,j) - h(i+2,j)
          vv = 2*v(i+1,j) - v(i+2,j)
          uu = 0.0
          gc = 1
        endif
        endif

        if (live(i,  j).eq.0) then ! right wall
        if (live(i-1,j).eq.1) then ! right wall
          hh = 2*h(i-1,j) - h(i-2,j)
          vv = 2*v(i-1,j) - v(i-2,j)
          uu = 0.0
          gc = 1
        endif
        endif

        if (live(i  ,j).eq.0) then  ! bottom wall
        if (live(i,j+1).eq.1) then  ! bottom wall
          hh = 2*h(i,j+1) - h(i,j+2)
          uu = 2*u(i,j+1) - u(i,j+2)
          vv = 0.0
          gc = 1
        endif
        endif

        if (live(i,  j).eq.0) then ! top wall
        if (live(i,j-1).eq.1) then ! top wall
          hh = 2*h(i,j-1) - h(i,j-2)
          uu = 2*u(i,j-1) - u(i,j-2)
          vv = 0.0
          gc = 1
        endif
        endif

        if (gc.eq.1) then ! ghost cell
        Q(i,j,1) =       hh
        Q(i,j,2) =    uu*hh
        Q(i,j,3) =    vv*hh
        E(i,j,1) =    uu*hh
        E(i,j,2) = uu*uu*hh + 0.5*gg*hh*hh
        E(i,j,3) = vv*uu*hh
        F(i,j,1) =    vv*hh
        F(i,j,2) = uu*vv*hh
        F(i,j,3) = vv*vv*hh + 0.5*gg*hh*hh
        endif

      enddo
      enddo



C-----step3 i-face vectors (half step)
      do j = 1, nj   ! j-cells
      do i = 1, ni+1 ! i-faces 1,ni+1
      if ((live(i,j).eq.1).or.(live(i-1,j).eq.1)) then

        hh = 0.5*(h(i,j)+h(i-1,j)) ! at i-face

        dbdx = (b(i,j) - b(i-1,j))/dx
        dbdy = 0.0

        sb(1) = 0.0
        sb(2) = dbdx
        sb(3) = dbdy

        do n = 1,3
          Qf(n) = 0.5*(Q(i,j,n) + Q(i-1,j,n))
     &     - 0.5*dtdx*(E(i,j,n) - E(i-1,j,n))
     &     - 0.5*dt*gg*hh*sb(n)
        enddo

        hh = Qf(1)
        uu = Qf(2)/hh
        vv = Qf(3)/hh
        
        Ef(i,j,1) =    uu*hh
        Ef(i,j,2) = uu*uu*hh + 0.5*gg*hh*hh
        Ef(i,j,3) = vv*uu*hh

      endif
      enddo
      enddo


C-----step4 j-face vectors (half step)
      do i = 1, ni   ! i-cells 
      do j = 1, nj+1 ! j-faces 1,nj+1
      if ((live(i,j).eq.1).or.(live(i,j-1).eq.1)) then

        hh = 0.5*(h(i,j)+h(i,j-1)) ! at j-face

C-------bed gradient term (x-direction)
        dbdx = (b(i+1,j) - b(i-1,j) + b(i+1,j-1) - b(i-1,j-1))/(4.0*dx)
        dbdy = 0.0

        sb(1) = 0.0
        sb(2) = dbdx
        sb(3) = dbdy

        do n = 1,3
          Qf(n) = 0.5*(Q(i,j,n) + Q(i,j-1,n))
     &     - 0.5*dtdy*(F(i,j,n) - F(i,j-1,n))
     &     - 0.5*dt*gg*hh*sb(n)
        enddo

        hh = Qf(1)
        uu = Qf(2)/hh
        vv = Qf(3)/hh

        Ff(i,j,1) =    vv*hh
        Ff(i,j,2) = uu*vv*hh
        Ff(i,j,3) = vv*vv*hh + 0.5*gg*hh*hh

      endif
      enddo
      enddo


C-----step5 Q-vector (full step)
      hmax = 0.0

      do i = 1, ni
      do j = 1, nj
      if (live(i,j).eq.1) then

        hh = h(i,j)
        dbdx = (b(i+1,j) - b(i-1,j))/(2.0*dx)
        dbdy = 0.0

        sb(1) = 0.0
        sb(2) = dbdx 
     &        - eps*(q(i+1,j,2) - 2.0*q(i,j,2) + q(i-1,j,2))/(dx*dx)
        sb(3) = dbdy
     &        - eps*(q(i,j+1,3) - 2.0*q(i,j,3) + q(i,j-1,3))/(dy*dy)

        do n = 1,3
          Q(i,j,n) = Q(i,j,n)
     &             - dtdx*(Ef(i+1,j,n) - Ef(i,j,n))
     &             - dtdy*(Ff(i,j+1,n) - Ff(i,j,n))
     &             - dt*gg*hh*sb(n)
        enddo

        h(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2)/h(i,j)
        v(i,j) = Q(i,j,3)/h(i,j)

        if (h(i,j)+b(i,j).gt.hmax) hmax = h(i,j)+b(i,j)


      endif
      enddo
      enddo

      enddo ! its

      endif ! sa > 0

      write(6,*)"sa=",sa,"hmax=",hmax

C-----write vtk file
      write(unit=string, fmt='(I3.3)') sa
      outfile = './results/shallow-water-'//string//'.vtk'
      OPEN(UNIT=20, FILE=outfile)
      WRITE(20,10)'# vtk DataFile Version 2.0'
      WRITE(20,10)'# sample rectilinear grid'
      WRITE(20,10)'ASCII'
      WRITE(20,10)'DATASET RECTILINEAR_GRID'
      WRITE(20,20)'DIMENSIONS ',ni+1,nj+1,2
      WRITE(20,30)'X_COORDINATES ',ni+1,' float'
      WRITE(20,*) (xn(i),i=1,ni+1)
      WRITE(20,30)'Y_COORDINATES ',nj+1,' float'
      WRITE(20,*) (yn(j),j=1,nj+1)
      WRITE(20,30)'Z_COORDINATES ',2,' float'
      WRITE(20,*) '0.0',dz
      WRITE(20,40)'CELL_DATA ',ni*nj
C-----scalar 
      WRITE(20,10)'SCALARS live float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((live(i,j),i=1,ni),j=1,nj)
C-----scalar 
      WRITE(20,10)'SCALARS ground float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((b(i,j),i=1,ni),j=1,nj)
C-----scalar 
      WRITE(20,10)'SCALARS height float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((h(i,j)+b(i,j),i=1,ni),j=1,nj)
C-----vector
      WRITE(20,10)'VECTORS velocity float'
      WRITE(20,*)((u(i,j),v(i,j),0.0,i=1,ni),j=1,nj)
      CLOSE(20)
C-----end write vtk file

      enddo ! nsav
C-----iterations

   10 FORMAT(A)
   20 FORMAT(A,3I4)
   30 FORMAT(A,I4,A)
   40 FORMAT(A,I9)

      END 






