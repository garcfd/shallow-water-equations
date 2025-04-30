
      program main
      implicit none

      integer, parameter :: ni=200,nj=100

      real xn(ni+1),yn(nj+1)
      real u(ni,nj),v(ni,nj),h(ni,nj)
      real Q(0:ni+1,0:nj+1,3)
      real E(0:ni+1,0:nj+1,3)
      real F(0:ni+1,0:nj+1,3)
      real Qf(3)
      real Ef(ni+1,nj,3)
      real Ff(ni,nj+1,3)

      integer im,jm,km,iloc,jloc
      integer i,j,ii,jj,it,sa,n

      real dx,dy,dt,rr,dtodx,dtody
      real gg,hh,uu,vv,hmax

      character(31) outfile
      character(3) string


C-----specify constants
      gg = 9.8  ! gravity
      dx = 1.0  ! spacing
      dy = 1.0  ! spacing
      dt = 0.01 ! timestep
      dtodx = dt/dx
      dtody = dt/dy
      hmax = 2.0

C-----create rectilinear mesh nodes
      xn(1) = 0.0 ! start point
      yn(1) = 0.0 ! start point
      do i = 2, ni+1
        xn(i) = xn(i-1) + dx
      enddo
      do j = 2, nj+1
        yn(j) = yn(j-1) + dy
      enddo


C-----initialise arrays with height pulse
      iloc =  0 ! i location of drop
      jloc = 50 ! j location of drop

      do i = 1, ni
      do j = 1, nj

        dx = real(i - iloc)
        dy = real(j - jloc)
        rr = sqrt(dx*dx + dy*dy)

        if (rr.lt.10) then        
          h(i,j) = 1.0 + 0.1*(10.0 - rr) ! 1.5m peak
        else
          h(i,j) = 1.0
        endif

        u(i,j) = 0.0 ! real(i)
        v(i,j) = 0.0 ! real(j)

      enddo
      enddo



C-----start iterations
      if (.true.) then
      do sa = 0,300

      if (sa.gt.0) then

      do it = 1,100

C-----step1 - vectors at cell centres
      do i = 1, ni
      do j = 1, nj

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

      enddo
      enddo
      

C-----step2 - vectors at ghost cells centres (imin/imax)
      do ii = 1,2
      do j = 1,nj

        if (ii.eq.1) then
          i  = 0
          hh = h(1,j) ! 2.0*h(1,j) - h(2,j)
          uu = 0.0    !    -u(1,j)
          vv = v(1,j) ! 2.0*v(1,j) - v(2,j)
        endif

        if (ii.eq.2) then
          i  = ni+1
          hh = h(ni,j) ! 2.0*h(ni,j) - h(ni-1,j)
          uu = 0.0     !    -u(ni,j)
          vv = v(ni,j) ! 2.0*v(ni,j) - v(ni-1,j)
        endif

        Q(i,j,1) =       hh
        Q(i,j,2) =    uu*hh
        Q(i,j,3) =    vv*hh

        E(i,j,1) =    uu*hh
        E(i,j,2) = uu*uu*hh + 0.5*gg*hh*hh
        E(i,j,3) = vv*uu*hh

        F(i,j,1) =    vv*hh
        F(i,j,2) = uu*vv*hh
        F(i,j,3) = vv*vv*hh + 0.5*gg*hh*hh

      enddo
      enddo
C----

C-----step3 - vectors at ghost cell centres (jmin/jmax)
      do jj = 1,2
      do i = 1,ni

        if (jj.eq.1) then ! jmin
          j  = 0
          hh = h(i,1) ! 2.0*h(i,1) - h(i,2)
          uu = u(i,1) ! 2.0*u(i,1) - u(i,2)
          vv = 0.0    !    -v(i,1)
        endif

        if (jj.eq.2) then ! jmax
          j  = nj+1
          hh = h(i,nj) ! 2.0*h(i,nj) - h(i,nj-1)
          uu = u(i,nj) ! 2.0*u(i,nj) - u(i,nj-1)
          vv = 0.0     !    -v(i,nj)
        endif

        Q(i,j,1) =       hh
        Q(i,j,2) =    uu*hh
        Q(i,j,3) =    vv*hh

        E(i,j,1) =    uu*hh
        E(i,j,2) = uu*uu*hh + 0.5*gg*hh*hh
        E(i,j,3) = vv*uu*hh

        F(i,j,1) =    vv*hh
        F(i,j,2) = uu*vv*hh
        F(i,j,3) = vv*vv*hh + 0.5*gg*hh*hh

      enddo
      enddo
C----


C-----step4 i-face vectors (half step)
      do j = 1, nj   ! j-cells
      do i = 1, ni+1 ! i-faces 1,ni+1

        do n = 1,3
          Qf(n) = 0.5*(Q(i,j,n) + Q(i-1,j,n))
     &    - 0.5*dtodx*(E(i,j,n) - E(i-1,j,n))
        enddo

        hh = Qf(1)
        uu = Qf(2)/hh
        vv = Qf(3)/hh

        Ef(i,j,1) =    uu*hh
        Ef(i,j,2) = uu*uu*hh + 0.5*gg*hh*hh
        Ef(i,j,3) = vv*uu*hh

      enddo
      enddo

C-----step5 j-face vectors (half step)
      do i = 1, ni   ! i-cells 
      do j = 1, nj+1 ! j-faces 1,nj+1

        do n = 1,3
          Qf(n) = 0.5*(Q(i,j,n) + Q(i,j-1,n))
     &    - 0.5*dtody*(F(i,j,n) - F(i,j-1,n))
        enddo

        hh = Qf(1)
        uu = Qf(2)/hh
        vv = Qf(3)/hh

        Ff(i,j,1) =    vv*hh
        Ff(i,j,2) = uu*vv*hh
        Ff(i,j,3) = vv*vv*hh + 0.5*gg*hh*hh

      enddo
      enddo

C-----step4 Q-vector (full step)
      hmax = 0.0
      do i = 1, ni
      do j = 1, nj

        do n = 1,3
          Q(i,j,n) = Q(i,j,n)
     &             - dtodx*(Ef(i+1,j,n) - Ef(i,j,n))
     &             - dtody*(Ff(i,j+1,n) - Ff(i,j,n))
        enddo

        h(i,j) = Q(i,j,1)
        u(i,j) = Q(i,j,2)/h(i,j)
        v(i,j) = Q(i,j,3)/h(i,j)

        if (h(i,j).gt.hmax) hmax = h(i,j)

      enddo
      enddo

      enddo ! its

      endif ! sa > 0

      write(6,*)"sa=",sa,"hmax=",hmax

C-----write vtk file
      write(unit=string, fmt='(I3.3)') sa
      outfile = './results/shallow-water-'//string//'.vtk'
CCC   write(6,*) "outfile = ", outfile
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
      WRITE(20,*) '0.0 1.0'
      WRITE(20,40)'CELL_DATA ',ni*nj
C-----scalar 
      WRITE(20,10)'SCALARS height float'
      WRITE(20,10)'LOOKUP_TABLE default'
      WRITE(20,*)((h(i,j),i=1,ni),j=1,nj)
C-----vector
      WRITE(20,10)'VECTORS velocity float'
      WRITE(20,*)((u(i,j),v(i,j),0.0,i=1,ni),j=1,nj)
      CLOSE(20)
C-----end write vtk file

      enddo ! nsav
      endif ! true/false
C-----iterations

   10 FORMAT(A)
   20 FORMAT(A,3I4)
   30 FORMAT(A,I3,A)
   40 FORMAT(A,I9)

      END 






