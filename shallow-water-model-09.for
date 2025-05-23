

      program main
      implicit none

      integer, parameter :: ni=200,nj=100

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

      integer i,j,n,it,sa,nits,nsav,gc
      integer iloc,jloc

      real dx,dy,dz,dt,dtdx,dtdy,dbdx,dbdy
      real gg,hh,uu,vv,hmax,eps,rr,uin,umax,vmax

      character(31) outfile
      character(3) string


C-----specify constants
      gg = 9.8     ! gravity
      dx = 0.1     ! spacing
      dy = 0.1     ! spacing
      dz = 0.01    ! spacing
      dt = 0.01    ! timestep
      eps = 0.1    ! smoothing (0.05)

      uin  = 0.5
      umax = 1.0
      vmax = 1.0

      hmax = 1.0
      nsav = 100
      nits = 100

      dtdx = dt/dx
      dtdy = dt/dy


C-----specify live/dead regions (if needed)
      live = 0            ! dead ghost cells needed for outer boundary
      live(1:ni,1:nj) = 1 ! live inner region

      iloc = 100 ! blockage centre
      jloc = 50  ! blockage centre
       
      do i = 1, ni
      do j = 1, nj
        rr = sqrt(real(i-iloc)**2 + real(j-jloc)**2)
        if (rr.le.10) then
          live(i,j) = 0 ! dead region
        endif
      enddo
      enddo



C-----create uniform mesh nodes
      xn(1) = 0.0 ! start point
      yn(1) = 0.0 ! start point

      do i = 2, ni+1
        xn(i) = xn(i-1) + dx
      enddo

      do j = 2, nj+1
        yn(j) = yn(j-1) + dy
      enddo


C-----initialise arrays with height pulse
      do i = 1, ni
      do j = 1, nj

        b(i,j) = 0.0 ! bed height
        h(i,j) = 0.1 ! water height
        u(i,j) = uin
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
          hh = 1*h(i+1,j) !- h(i+2,j)
          vv = 1*v(i+1,j) !- v(i+2,j)
          uu = 0.0
          gc = 1
        endif
        endif

        if (live(i,  j).eq.0) then ! right wall
        if (live(i-1,j).eq.1) then ! right wall
          hh = 1*h(i-1,j) !- h(i-2,j)
          vv = 1*v(i-1,j) !- v(i-2,j)
          uu = 0.0
          gc = 1
        endif
        endif

        if (live(i  ,j).eq.0) then  ! bottom wall
        if (live(i,j+1).eq.1) then  ! bottom wall
          hh = 1*h(i,j+1) !- h(i,j+2)
          uu = 1*u(i,j+1) !- u(i,j+2)
          vv = 0.0
          gc = 1
        endif
        endif

        if (live(i,  j).eq.0) then ! top wall
        if (live(i,j-1).eq.1) then ! top wall
          hh = 1*h(i,j-1) !- h(i,j-2)
          uu = 1*u(i,j-1) !- u(i,j-2)
          vv = 0.0
          gc = 1
        endif
        endif

        if ((i.eq.0).and.(live(1,j).eq.1)) then ! xmin inlet velocity
          gc = 1
          uu = uin
        endif

        if ((i.eq.ni+1).and.(live(ni,j).eq.1)) then ! xmax outlet boundary
          gc = 1
          uu = u(i-1,j)
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

        dbdx = 0.0 ! (b(i,j) - b(i-1,j))/dx
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
        dbdx = 0.0 ! (b(i+1,j) - b(i-1,j) + b(i+1,j-1) - b(i-1,j-1))/(4.0*dx)
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
     &        - eps*(q(i+1,j,1) - 2.0*q(i,j,1) + q(i-1,j,1))/(dx*dx)
     &        - eps*(q(i,j+1,1) - 2.0*q(i,j,1) + q(i,j-1,1))/(dy*dy)
        sb(2) = dbdx 
     &        - eps*(q(i+1,j,2) - 2.0*q(i,j,2) + q(i-1,j,2))/(dx*dx)
     &        - eps*(q(i,j+1,2) - 2.0*q(i,j,2) + q(i,j-1,2))/(dy*dy)
        sb(3) = dbdy
     &        - eps*(q(i+1,j,3) - 2.0*q(i,j,3) + q(i-1,j,3))/(dx*dx)
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

c-------velcoity limiter (if needed)
        if (u(i,j).gt.umax) u(i,j) = umax
        if (v(i,j).gt.vmax) v(i,j) = vmax

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






