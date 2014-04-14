cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
c
c
c    $Date$
c    $Revision$
c
c
c     ROTATION VIA PROJECTION  (FORTRAN 77 AND 90 VERSIONS).
c
c     Requires FFT and Associated Legendre Function Libraries.
c
c     User-callable f77 routine is rotviaprojvar. 
c     User-callable f90 routine is rotviaprojvarf90.
c     The other routines are used internally.
c
c
c***********************************************************************
ccc      subroutine rotviaprojvar(alpha,nterms,m1,m2,mpole,lmp,marray2,
      subroutine rotviaproj(alpha,nterms,m1,m2,mpole,lmp,marray2,
     1           lmpn,w,lw,lused)
c***********************************************************************
c     Purpose:
c
c     Fast and stable algorithm for applying rotation operator about
c     the y-axis determined by angle alpha.
c
c     There is some loss in speed over using recurrence relations 
c     but it is stable to all orders whereas the recurrence schemes 
c     are not.
c     If the rotation operator is to be used multiple times, and
c     memory is available, one can precompute and store the 
c     multipliers used in evalall (see below). This has not yet been
c     implemented.
c
c     In the simplest version of the method (rotproj.f), projection
c     is carried out on the equator of the rotated sphere.
c
c     In the present version, a sequence of constant latitude circles
c     are used: 
c
c     On the equator, the dynamic range of Y_n^m(cos theta) is 
c     modest but the number of sampling nodes needed is 2*n+2 for 
c     terms of degree n. Thus, nterms**3 work is needed even if one
c     only wants orders of m up to m2 << nterms in the rotated expansion.
c     It is possible to do this with approx. (nterms**2)*m2 work, using
c     a more complicated scheme that requires some quadrature ideas
c     which we have tabulated in an approximate sense.
c
c     In essence, the idea is that on constant latitude circles closer 
c     to the poles, the dynamic range of Y_n^m(cos theta) increases over 
c     that on the equator, but high "m" modes decay rapidly. Thus,
c     the number of quadrature nodes can be less than 2*n+2 for any 
c     specified finite precision.
c     This is a rather subtle numerical issue, having to do with 
c     sampling of functions which are not strictly speaking
c     band-limited but for which the high frequency modes are 
c     significanty attenuated. We have chosen to pick values for which
c     13 digits or better are achieved (making no attempt to be more
c     efficient for lower precision levels).
c       
c     The optimal choice of circle is difficult to determine 
c     analyically, but we have experimentally found a reasonable
c     table of theta values (thetaprojs(i)) and nodes (nquads(i)) 
c     for varying values of n and m2. 
c     These are currently defined at the beginning of the code.
c     
c     For each of the fixed circles, the coefficients of the 
c     rotated expansion can be obtained by projection using the FFT. 
c
C---------------------------------------------------------------------
c     INPUT:
c
c     alpha:  the rotation angle about the y-axis.
c     nterms: order of multipole expansion
c
c     m1    : max m index for first expansion  
c     m2    : max m index for second expansion 
C     mpole   coefficients of original multiple expansion
C     lmp     leading dim for mpole (must exceed nterms)
C     lmpn    leading dim for marray2 (must exceed nterms)
c     w     :  work array 
c     lw    :  length of work array 
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     marray2  coefficients of rotated expansion.
c     lused    amount of workspace used.
c     ier      error return flag
c                0 successful execution
c                1 insufficient memory
c
c---------------------------------------------------------------------
c
      implicit real *8 (a-h,o-z)
      integer nquads(4),n2s(4)
      real *8 w(lw),pi,thetaprojs(4)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
c
c
c     if m2 is sufficiently large, just use equator.
c     This covers the case m1 < m2 ("small to big").
c
      pi = 4.0d0*datan(1.0d0)
      n2s(1) = nterms
      n2s(2) = nterms
      n2s(3) = nterms
      n2s(4) = nterms
      nquads(1) = 2*nterms+2
      nquads(2) = 2*nterms+2
      nquads(3) = 2*nterms+2
      nquads(4) = 2*nterms+2
      thetaprojs(1) = pi/2.0d0
      thetaprojs(2) = pi/2.0d0
      thetaprojs(3) = pi/2.0d0
      thetaprojs(4) = pi/2.0d0
c
c     now optimize for smaller and smaller values of m2
c     assuming m1 > m2.
c
      if ((m1.ge.m2).and.(m2.le.160)) then
         m3 = 96
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 12*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(500) )
         thetaprojs(3) = asin( 5*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 5*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.80)) then
         m3 = 48
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 12*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(250) )
         thetaprojs(3) = asin( 5*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 5*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.40)) then
         m3 = 36
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 8*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(200) )
         thetaprojs(3) = asin( 4*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 4*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.20)) then
         m3 = 24
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 8*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(200) )
         thetaprojs(3) = asin( 4*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 4*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.10)) then
         m3 = 12
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 8*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(200) )
         thetaprojs(3) = asin( 4*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 4*dfloat(m3)/dfloat(2400) )
      endif
c      
      do ii = 1,4
         if (ii.eq. 1) then
            n1 = 0
            n2 = n2s(ii)
            thetaproj = thetaprojs(ii)
            nquad = nquads(ii)
         else if (ii.ge.2) then
            n1 = n2+1
            n2 = n2s(ii)
            if (m1 .gt. 4*m2) then
               thetaproj = thetaprojs(ii)
               nquad = nquads(ii)
            else
               thetaproj = pi/2.0d0
               nquad = 2*nterms+2
            endif
         endif     
c
c     compute rotation for each set of degrees [n1,n2]
c     on successively greater latitude circles.
c
         ictheta = 1
         istheta = ictheta+nquad
         icphi = istheta+nquad
         isphi = icphi+nquad
         iynm = isphi+nquad
         iynmd = iynm + (nterms+1)**2
         irat1 = iynmd + (nterms+1)**2
         irat2 = irat1 + (nterms+1)**2
         iuval = irat2 + (nterms+1)**2
         iuder = iuval + 2*nquad*(nterms+1)
         iephi = iuder + 2*nquad*(nterms+1)
         iwsave = iephi + 2*(2*nterms+1)
         iavec = iwsave + 4*nquad+20
         ibvec = iavec + 2*nquad
         lused = ibvec + 2*nquad
         if (lused.gt.lw) stop
c
         call rotviaprojvar0(alpha,nquad,thetaproj,nterms,n1,n2,m1,m2,
     1           mpole,lmp,marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),w(iuval),w(iuder),
     1           w(iephi),w(iwsave),w(iavec),w(ibvec))
         if (n2.ge.nterms) return
      enddo
c
      return
      end
c
c
c
c
c
c***********************************************************************
      subroutine rotviaprojvar0(alpha,nquad,thetaproj,nterms,n1,n2,
     1       m1,m2,mpole,lmp,marray2,lmpn,cthetas,sthetas,cphis,sphis,
     1       ynm,ynmd,rat1,rat2,uval,uder,ephis,wsave,avec,bvec)
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad,nterms
      real *8 cthetas(nquad),cphis(nquad), thetaproj
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      real *8 rat1(0:nterms,0:nterms)
      real *8 rat2(0:nterms,0:nterms)
      complex *16 avec(nquad)
      complex *16 bvec(nquad)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 ephis(-nterms:nterms)
      real *8 wsave(4*nquad+20)
C
C
      call getmeridianvar(alpha,nquad,thetaproj,
     1     cthetas,sthetas,cphis,sphis)    
      call evalallvar(alpha,nquad,thetaproj,n1,n2,m1,cthetas,sthetas,
     1     cphis,sphis,mpole,lmp,nterms,uval,uder,ynm,ynmd,ephis,
     2     rat1,rat2)
      call projectonynmvar(nquad,thetaproj,uval,uder,ynm,ynmd,marray2,
     1     lmpn,nterms,n1,n2,m2,wsave,avec,bvec)
      return
      end
C
C
C***********************************************************************
      subroutine getmeridianvar(alpha,nquad,thetaproj,cthetas,sthetas,
     1           cphis,sphis)
C***********************************************************************
C     Purpose:
C
C     For a rotation of angle ALPHA about the y-axis, this
C     subroutine returns the NQUAD equispaced nodes on a circle 
C     of constant latitude (on or above the equator) in the rotated 
c     frame, but in the original coordinate system.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     alpha     = angle of rotation
C     nquad     = number of quadrature nodes in equator.
c     thetaproj = theta value that defines latitude in rotated frame
c                 on which evaluation is being carried out.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     cthetas = cos(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     sthetas = sin(theta) values in original coordinate system of 
C                    nquad equispaced nodes
C     cphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C     sphis =  cos(phi) values in original coordinate system of 
C                    nquad equispaced nodes
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad
      real *8 cthetas(nquad)
      real *8 sthetas(nquad)
      real *8 cphis(nquad)
      real *8 sphis(nquad)
C
      pi = 4.0d0*datan(1.0d0)
C
      ca = cos(alpha)
      sa = sin(alpha)
      do i = 1,nquad
	 im1 = i-1
         phi = 2*pi*im1/nquad
         xp = cos(phi)*sin(thetaproj)
         yp = sin(phi)*sin(thetaproj)
         zp = cos(thetaproj)
ccc	 write(19,*) xp,yp,zp
         x = ca*xp + sa*zp
         y = yp
         z = -sa*xp + ca*zp
ccc	 write(20,*) x,y,z
         proj = sqrt(x**2+y**2)
	 if (proj.le.1.0d-16) then
	    cphis(i) = 1.0d0
	    sphis(i) = 0.0d0
	 else
	    cphis(i) = x/proj
	    sphis(i) = y/proj
	 endif
	 cthetas(i) = z
	 sthetas(i) = proj
      enddo
      return
      end
C
C***********************************************************************
      subroutine evalallvar(alpha,nquad,thetaproj,n1,n2,m1,
     1           cthetas,sthetas,cphis,sphis,mpole,lmp,nterms,
     2           uval,uder,ynm,ynmd,ephis,rat1,rat2)
C***********************************************************************
C
C     This subroutine evaluates the multipole expansion for each
C     order at the nquad nodes on the rotated equator.
C
C---------------------------------------------------------------------
C     INPUT:
C
C     alpha    : angle of rotation about y-axis.
C     nquad    : number of target point son unit sphere
c     thetaproj : theta value that defines latitude on which evaluation
c                 is being carried out.
c     n1       : lowest degree to be rotated
c     n2       : highest degree to be rotated
c     m1       : maximum order (second index) in mpole
C     cthetas  : cos(theta) values of target points.
C     sthetas  : sin(theta) values of target points.
C     cphis    : cos(phi) values of target points.
C     sphis    : sin(phi) values of target points.
C     mpole    : original multipole expansion
C     nterms   : order of multipole expansion
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     ephis    : work array for exp(i m phi) values
C     rat1     : work array for accelerating ynm calculation.
C     rat2     : work array for accelerating ynm calculation.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C     uval(i,j) : contribution to potential 
C                 of multipole terms of order j at ith quad node.
C     uder(i,j) : contributions to theta derivative of potential
C                 of multipole terms of order j at ith quad node.
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer ndeg,morder, nquad
      real *8 cthetas(nquad),cphis(nquad)
      real *8 sthetas(nquad),sphis(nquad)
      real *8 ynm(0:n2,0:n2)
      real *8 ynmd(0:n2,0:n2)
      real *8 thetahat(3)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 ephi1,ephis(-nterms:nterms)
      complex *16 uder(nquad,0:nterms),uval(nquad,0:nterms)
      complex *16 uv,utheta,uphi,ztmp1,ztmp2,ztsum
      complex *16 ux,uy,uz,imag
      real *8 rat1(0:n2,0:n2)
      real *8 rat2(0:n2,0:n2)
C
      data imag/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      calpha = cos(alpha)
      salpha = -sin(alpha)
      call ylgndrini(n2,rat1,rat2)
      do jj=1,nquad
	 ctheta = cthetas(jj)
	 stheta = sthetas(jj)
	 cphi = cphis(jj)
	 sphi = sphis(jj)
         dir1 = -salpha
         dir2 = 0
         dir3 = calpha
ccc
	 phiold = 2*pi*(jj-1)/nquad
	 x0 = cos(phiold)*cos(thetaproj)
	 y0 = sin(phiold)*cos(thetaproj)
	 z0 = -sin(thetaproj)
         dir1 = calpha*x0 - salpha*z0 
         dir2 = y0
         dir3 = salpha*x0 + calpha*z0
         dir1 = -dir1
         dir2 = -dir2
         dir3 = -dir3
ccc
         tang1 = cphi*ctheta
         tang2 = sphi*ctheta
         tang3 = -stheta
         proj2 = tang1*dir1 + tang2*dir2 + tang3*dir3
         tang1 = -sphi
         tang2 = cphi
         tang3 = 0
         proj1 = tang1*dir1 + tang2*dir2 + tang3*dir3
         thetax = ctheta*cphi
         thetay = ctheta*sphi
         thetaz = -stheta
         phix = -sphi
         phiy = cphi
         phiz = 0.0d0
         thetahat(1) = -salpha
         thetahat(2) = 0
         thetahat(3) = -calpha
	 call ylgndr2sf_trunc(n2,m1,ctheta,ynm,ynmd,rat1,rat2)
         ephi1 = dcmplx(cphis(jj),sphis(jj))
	 ephis(1) = ephi1
	 ephis(-1) = dconjg(ephi1)
	 do i = 2,nterms
	    ephis(i) = ephis(i-1)*ephi1
	    ephis(-i) = ephis(-i+1)*dconjg(ephi1)
	 enddo
c
c	 
	 do ndeg = n1,n2
	    uv = ynm(ndeg,0)*mpole(ndeg,0)
	    utheta=ynmd(ndeg,0)*stheta*mpole(ndeg,0)
	    uphi = 0.0d0
	    do morder = 1,min(ndeg,m1)
               ztmp1 = ephis(morder)*mpole(ndeg,morder)
               ztmp2 = ephis(-morder)*mpole(ndeg,-morder)
	       ztsum = ztmp1+ztmp2
	       uv = uv + stheta*ynm(ndeg,morder)*ztsum
	       utheta = utheta + ynmd(ndeg,morder)*ztsum
	       uphi = uphi - ynm(ndeg,morder)*morder*(ztmp1-ztmp2)
	    enddo
	    ux = utheta*thetax + imag*uphi*phix
	    uy = utheta*thetay + imag*uphi*phiy
	    uz = utheta*thetaz + imag*uphi*phiz
	    ux = utheta*thetax 
	    uy = utheta*thetay 
	    uz = utheta*thetaz
            uval(jj,ndeg) = uv
            uder(jj,ndeg) = utheta*proj2+uphi*imag*proj1
	 enddo
      enddo
      return
      end
C
C
C
C
C
C***********************************************************************
      subroutine projectonynmvar(nquad,thetaproj,uval,uder,
     1           ynm,ynmd,marray,lmpn,nterms,n1,n2,m2,wsave,avec,bvec)
C***********************************************************************
C
C     This subroutine projects from values on equator for each multipole
C     order (uval, uder = dudthteta) 
C     onto spherical harmonics
C
C---------------------------------------------------------------------
C     INPUT:
C
C     nquad    : number of points on equator
c     thetaproj : theta value that defines latitude on which evaluation
c                 is being carried out.
C     uval     : F values on equator
C     uder     : dFdtheta values on equator
C     ynm      : work array for ynm values
C     ynmd     : work array for ynmd values
C     lmpn     : leading dim of marray (must exceed nterms)
C     nterms   : order of expansion
C     n1       : lowest degree to be projected
C     n2       : highest degree to be projected
C     m2       : maximum order (second index) computed in marray.
C     wsave    : work array for FFT (dimension at least 4*nquad+20)
C     avec     : work array of length nquad for FFT (complex)
C     bvec     : work array of length nquad for FFT (complex)
C---------------------------------------------------------------------
C     OUTPUT:
C
C     marray   : rotated expansion 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nquad, norder
      complex *16 ephi,ephi1,uval(nquad,0:1)
      complex *16 uder(nquad,0:1)
      complex *16 utheta,uphi,ztmp1,ztmp2
      complex *16 alpha,beta,imag
      complex *16 marray(0:lmpn,-lmpn:lmpn)
      real *8 ynm(0:n2,0:n2),ynms,ynmds
      real *8 ynmd(0:n2,0:n2)
      real *8 wsave(4*nquad+20)
      complex *16 avec(nquad)
      complex *16 bvec(nquad)
C
      data imag/(0.0d0,1.0d0)/
      pi = 4.0d0*datan(1.0d0)
C
      do norder=n1,n2
         do m=-norder,norder
            marray(norder,m) = 0.0d0
         enddo
      enddo
      ctheta = cos(thetaproj)
      stheta = sin(thetaproj)
      h = 1.0d0/nquad
      call ylgndr2s_trunc(n2,m2,ctheta,ynm,ynmd)
      call zffti(nquad,wsave)
      do norder=n1,n2
	 do ii = 1,nquad
	    avec(ii) = uval(ii,norder)
	    bvec(ii) = uder(ii,norder)
         enddo
	 call zfftf(nquad,avec,wsave)
	 call zfftf(nquad,bvec,wsave)
         do m = -min(m2,norder),min(m2,norder)
	    if (m.ge.0)  alpha = avec(m+1)*h
	    ynms = ynm(norder,abs(m))
	    ynmds = ynmd(norder,abs(m))
	    if (m.ne.0)  ynms = ynms*stheta
	    if (m.eq.0)  ynmds = ynmds*stheta
	    if (m.lt.0)  alpha = avec(nquad+m+1)*h
	    if (m.ge.0)  beta = bvec(m+1)*h
	    if (m.lt.0)  beta = bvec(nquad+m+1)*h
            beta = beta/(norder+0.5d0)
            ynmds = ynmds/(norder+0.5d0)
            marray(norder,m) = (alpha*ynms - beta*ynmds)/
     1        (ynms**2 + ynmds**2)
         enddo
      enddo
      return
      end

c
c
c***********************************************************************
      subroutine rotviaprojf90
     $     (alpha,nterms,m1,m2,mpole,lmp,marray2,lmpn)
c***********************************************************************
c       Purpose:
c
c	Fast and stable algorithm for applying rotation operator about
c	the y-axis determined by angle alpha.
c
c       The method is based on computing the induced potential and
c       its theta-derivative on the rotated equator
c       for each order (first index). The coefficients of  the rotated
c       expansion can then be obtained by FFT and projection.
c
c       There is some loss in speed over using recurrence relations 
c       but it is stable to all orders whereas the recurrence schemes 
c       are not.
c       If the rotation operator is to be used multiple times, and
c       memory is available, one can precompute and store the 
c       multipliers used in evalall (see below). This has not yet been
c       implemented.
c
C---------------------------------------------------------------------
c       INPUT:
c
c       alpha:  the rotation angle about the y-axis.
c       nterms: order of multipole expansion
C       mpole   coefficients of original multiple expansion
C       lmp     leading dim for mpole (must exceed nterms)
C       lmpn    leading dim for marray2 (must exceed nterms)
c
C---------------------------------------------------------------------
c       OUTPUT:
c
c       marray2  coefficients of rotated expansion.
c
C---------------------------------------------------------------------
c
      implicit real *8 (a-h,o-z)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
c
      call rotviaprojvarf90
     $   (alpha,nterms,nterms,nterms,mpole,lmp,marray2,lmpn)
c
      return
      end
c
c
c***********************************************************************
      subroutine rotviaprojvarf90(alpha,nterms,m1,m2,mpole,lmp,marray2,
     1           lmpn)
c***********************************************************************
c     Purpose:
c
c     Fast and stable algorithm for applying rotation operator about
c     the y-axis determined by angle alpha.
c
c     There is some loss in speed over using recurrence relations 
c     but it is stable to all orders whereas the recurrence schemes 
c     are not.
c     If the rotation operator is to be used multiple times, and
c     memory is available, one can precompute and store the 
c     multipliers used in evalall (see below). This has not yet been
c     implemented.
c
c     In the simplest version of the method (rotproj.f), projection
c     is carried out on the equator of the rotated sphere.
c
c     In the present version, a sequence of constant latitude circles
c     are used: 
c
c     On the equator, the dynamic range of Y_n^m(cos theta) is 
c     modest but the number of sampling nodes needed is 2*n+2 for 
c     terms of degree n. Thus, nterms**3 work is needed even if one
c     only wants orders of m up to m2 << nterms in the rotated expansion.
c     It is possible to do this with approx. (nterms**2)*m2 work, using
c     a more complicated scheme that requires some quadrature ideas
c     which we have tabulated in an approximate sense.
c
c     In essence, the idea is that on constant latitude circles closer 
c     to the poles, the dynamic range of Y_n^m(cos theta) increases over 
c     that on the equator, but high "m" modes decay rapidly. Thus,
c     the number of quadrature nodes can be less than 2*n+2 for any 
c     specified finite precision.
c     This is a rather subtle numerical issue, having to do with 
c     sampling of functions which are not strictly speaking
c     band-limited but for which the high frequency modes are 
c     significanty attenuated. We have chosen to pick values for which
c     13 digits or better are achieved (making no attempt to be more
c     efficient for lower precision levels).
c       
c     The optimal choice of circle is difficult to determine 
c     analyically, but we have experimentally found a reasonable
c     table of theta values (thetaprojs(i)) and nodes (nquads(i)) 
c     for varying values of n and m2. 
c     These are currently defined at the beginning of the code.
c     
c     For each of the fixed circles, the coefficients of the 
c     rotated expansion can be obtained by projection using the FFT. 
c
C---------------------------------------------------------------------
c     INPUT:
c
c     alpha:  the rotation angle about the y-axis.
c     nterms: order of multipole expansion
c
c     m1    : max m index for first expansion  
c     m2    : max m index for second expansion 
C     mpole   coefficients of original multiple expansion
C     lmp     leading dim for mpole (must exceed nterms)
C     lmpn    leading dim for marray2 (must exceed nterms)
c     w     :  work array 
c     lw    :  length of work array 
c
C---------------------------------------------------------------------
c     OUTPUT:
c
c     marray2  coefficients of rotated expansion.
c     lused    amount of workspace used.
c     ier      error return flag
c                0 successful execution
c                1 insufficient memory
c
c---------------------------------------------------------------------
c
      implicit real *8 (a-h,o-z)
      integer nquads(4),n2s(4),ier
      real *8 pi,thetaprojs(4)
      real *8 w
      allocatable w(:)
      complex *16 mpole(0:lmp,-lmp:lmp)
      complex *16 marray2(0:lmpn,-lmpn:lmpn)
c
c
c     if m2 is sufficiently large, just use equator.
c     This covers the case m1 < m2 ("small to big").
c
      pi = 4.0d0*datan(1.0d0)
      n2s(1) = nterms
      n2s(2) = nterms
      n2s(3) = nterms
      n2s(4) = nterms
      nquads(1) = 2*nterms+2
      nquads(2) = 2*nterms+2
      nquads(3) = 2*nterms+2
      nquads(4) = 2*nterms+2
      thetaprojs(1) = pi/2.0d0
      thetaprojs(2) = pi/2.0d0
      thetaprojs(3) = pi/2.0d0
      thetaprojs(4) = pi/2.0d0
c
c     now optimize for smaller and smaller values of m2
c     assuming m1 > m2.
c
      if ((m1.ge.m2).and.(m2.le.160)) then
         m3 = 96
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 12*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(500) )
         thetaprojs(3) = asin( 5*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 5*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.80)) then
         m3 = 48
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 12*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(250) )
         thetaprojs(3) = asin( 5*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 5*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.40)) then
         m3 = 36
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 8*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(200) )
         thetaprojs(3) = asin( 4*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 4*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.20)) then
         m3 = 24
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 8*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(200) )
         thetaprojs(3) = asin( 4*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 4*dfloat(m3)/dfloat(2400) )
      endif
      if ((m1.ge.m2).and.(m2.le.10)) then
         m3 = 12
         n2s(1) = 2*m3
	 n2s(1) = min(n2s(1),nterms)
         n2s(2) = 8*m3
	 n2s(2) = min(n2s(2),nterms)
         n2s(3) = 1200
	 n2s(3) = min(n2s(3),nterms)
         n2s(4) = nterms
         nquads(1) = 2*n2s(1)+2
         nquads(2) = 8*m3+2
         nquads(3) = 8*m3+2
         nquads(4) = 8*m3+2
         thetaprojs(1) = pi/2.0d0
         thetaprojs(2) = asin( 4*dfloat(m3)/dfloat(200) )
         thetaprojs(3) = asin( 4*dfloat(m3)/dfloat(1200) )
         thetaprojs(4) = asin( 4*dfloat(m3)/dfloat(2400) )
      endif
c      
      do ii = 1,4
         if (ii.eq. 1) then
            n1 = 0
            n2 = n2s(ii)
            thetaproj = thetaprojs(ii)
            nquad = nquads(ii)
         else if (ii.ge.2) then
            n1 = n2+1
            n2 = n2s(ii)
            if (m1 .gt. 4*m2) then
               thetaproj = thetaprojs(ii)
               nquad = nquads(ii)
            else
               thetaproj = pi/2.0d0
               nquad = 2*nterms+2
            endif
         endif     
c
c     compute rotation for each set of degrees [n1,n2]
c     on successively greater latitude circles.
c
         ictheta = 1
         istheta = ictheta+nquad
         icphi = istheta+nquad
         isphi = icphi+nquad
         iynm = isphi+nquad
         iynmd = iynm + (nterms+1)**2
         irat1 = iynmd + (nterms+1)**2
         irat2 = irat1 + (nterms+1)**2
         iuval = irat2 + (nterms+1)**2
         iuder = iuval + 2*nquad*(nterms+1)
         iephi = iuder + 2*nquad*(nterms+1)
         iwsave = iephi + 2*(2*nterms+1)
         iavec = iwsave + 4*nquad+20
         ibvec = iavec + 2*nquad
         lused = ibvec + 2*nquad
         allocate (w(lused), stat=ier)
         if (ier.ne.0) then
            write(6,*) ' alloc failure in rotviaprojvarf90'
            stop
         endif
c
         call rotviaprojvar0(alpha,nquad,thetaproj,nterms,n1,n2,m1,m2,
     1           mpole,lmp,marray2,lmpn,w(ictheta),w(istheta),
     1           w(icphi),w(isphi),w(iynm),w(iynmd),
     1           w(irat1),w(irat2),w(iuval),w(iuder),
     1           w(iephi),w(iwsave),w(iavec),w(ibvec))
         deallocate(w)
         if (n2.ge.nterms) return
      enddo
      return
      end
c
c
