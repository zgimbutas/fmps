cc Copyright (C) 2009-2010: Leslie Greengard and Zydrunas Gimbutas
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
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c     multipole shift routines, f95 version, using allocate
c
C***********************************************************************
      subroutine h3dmpmpquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for mpole expansion
C           x0y0z0 = center of original multiple expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           sc2     = scaling parameter for shifted expansion
C           xnynzn = center of shifted expansion
C           nterms2 = order of shifted expansion
C           radius  = radius of sphere on which mpole expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta direction
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted mpole expansion
C           ier   = error return flag
C
C                   CURRENTLY UNUSED.
C
C     Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           ynm      = work array 
C           phitemp  = work array 
C           fhs      = work array 
C           fhder    = work array 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
c
c     local allocated workspace array
c
      real *8, allocatable :: w(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ldc = max(nterms,nterms2)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2  
      iynmd = iynm+lynm
      iphitemp = iynmd+lynm
      lphitemp = nq*(2*ldc+1)*2 + 7  
      iphitemp2 = iphitemp+lphitemp
      ifhs = iphitemp2+lphitemp
      ifhder = ifhs + 2*(nterms+1) + 3
      lused = ifhder + 2*(nterms+1) +3 + 100
      allocate (w(lused))
c
ccc      call prinf(' ier before mpmpquad0 is *',ier,1)
      call h3dmpmpquad0(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,nterms2,w(imarray),w(imarray1),
     2           ldc,w(iephi),
     3           radius,xnodes,wts,nquad,nq,w(iynm),
     4           w(iynmd),w(iphitemp),w(iphitemp2),
     7           w(ifhs),w(ifhder),ier)
      return
      end
c
c
c
C***********************************************************************
      subroutine h3dmpmpquadu_add(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mpolen,ldc,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for mpole expansion
C           x0y0z0 = center of original multiple expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           sc2     = scaling parameter for shifted expansion
C           xnynzn = center of shifted expansion
C           nterms2 = order of shifted expansion
C           radius  = radius of sphere on which mpole expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta direction
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted mpole expansion
C           ier   = error return flag
C
C                   CURRENTLY UNUSED.
C
C     Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           ynm      = work array 
C           phitemp  = work array 
C           fhs      = work array 
C           fhder    = work array 
C
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:ldc,-ldc:ldc)
      complex *16 wavek,imag
c
c     local allocated workspace array
c
      complex *16, allocatable :: mptemp(:,:)
c
      data imag/(0.0d0,1.0d0)/
C
      allocate( mptemp(0:nterms2,-nterms2:nterms2) )

      call h3dmpmpquadu(wavek,sc1,x0y0z0,mpole,nterms,
     1           sc2,xnynzn,mptemp,nterms2,
     2           radius,xnodes,wts,nquad,ier)

      do l = 0,min(ldc,nterms2)
         do m=-l,l
            mpolen(l,m) = mpolen(l,m)+mptemp(l,m)
         enddo
      enddo
c
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmpmpquad0(wavek,sc1,x0y0z0,mpole,nterms,sc2,
     1           xnynzn,mpolen,nterms2,marray,marray1,ldc,ephi,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,
     3           phitemp,phitemp2,fhs,fhder,ier)
C***********************************************************************
C
C     Usage:
C
C           Shift multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of shifted expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for shifted expansion
C           radius  = radius of sphere on which shifted expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta
C           nq      = max(nquad, 2*ldc+2) = used to allocate
C                     work arrays for both z-shift and rotations.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted expansion
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           ldc      determines dimension of dc
c                    must exceed max(nterms,nterms2).
C           rd     = work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer  nterms, lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 radius, rshift
      real *8 xnodes(1),wts(1)
      real *8 d,theta,ctheta,phi,sc1,sc2,rvec(3)
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitemp2(nq,-ldc:ldc)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
c
      complex *16 ephi(-ldc-1:ldc+1),imag
ccc      complex *16 ephi2(-ldc-1:ldc+1)
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polar(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c----- a rotation of THETA radians about the Yprime axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      mpole and ephi inside the following loop. 
c
      do l=0,nterms
         do m=-l,l
            marray1(l,m)=mpole(l,m)*ephi(m)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms .ge. 30 ) then 
      call rotviaprojf90(theta,nterms,nterms,nterms,marray1,nterms,
     1        marray,ldc)
      else
      call rotviarecur3f90(theta,nterms,nterms,nterms,marray1,nterms,
     1        marray,ldc)
      endif
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmpmpzshift_fast
     $   (wavek,sc1,marray,ldc,nterms,sc2,mpolen,
     1           nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2           ynm,phitemp,fhs,fhder,ier)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      else
ccc      call rotviarecur3f90(-theta,nterms2,nterms,nterms2,mpolen,
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmpmpquadu_trunc(wavek,sc1,x0y0z0,mpole,nterms,
     1           nterms1,sc2,xnynzn,mpolen,nterms2,
     2           radius,xnodes,wts,nquad,ier)
C***********************************************************************
C
C     memory management wrapper for 
C     subroutine h3dmpmpquad0 (below).
C
C     Usage:
C
C           Shift center of multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then shifts
C           along the Z-axis, and then rotates back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           sc1     = scaling parameter for mpole expansion
C           x0y0z0 = center of original multiple expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           sc2     = scaling parameter for shifted expansion
C           xnynzn = center of shifted expansion
C           nterms2 = order of shifted expansion
C           radius  = radius of sphere on which mpole expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta direction
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted mpole expansion
C
C           ier   = error return flag
C
C                   CURRENTLY NOT USED
C
C     Work arrays carved out of w.
C
C           marray   = work array used to hold various intermediate 
c                      rotated expansions.
C           dc       = work array contain the square roots of 
C                      some binomial coefficients.
C           rd1,rd2  = work arrays used to compute rotation matrices
C                      about Y-axis recursively.
C           ephi     = work array 
C           ynm      = work array 
C           phitemp  = work array 
C           fhs      = work array 
C           fhder    = work array 
C
C
C     WORK ESTIMATE:
C
C     Let nmax = max(nterms,nterms2). Then
C     w should be at least 9*(nmax)^2+ 27*(nmax) 
C                          + 8*nquad*nmax + 2*nquad + 100
C
C***********************************************************************
      implicit real *8 (a-h,o-z)
      integer  nterms,ier,l,m,jnew,knew
      real *8 x0y0z0(3),xnynzn(3)
      real *8 xnodes(1),wts(1)
      real *8 sc1,sc2
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 wavek,imag
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: w(:)
c
      data imag/(0.0d0,1.0d0)/
C
      ier = 0
      ldc = max(nterms,nterms2)
      ldc = max(ldc,nterms1)
      nq = max(nquad,2*ldc+2)
      imarray = 1
      lmarray = 2*(ldc+1)*(2*ldc+1) + 3 
      imarray1 = imarray+lmarray
      lmarray1 = 2*(ldc+1)*(2*ldc+1) + 3 
      iephi = imarray1+lmarray1
      lephi = 2*(2*ldc+3) + 3 
      iynm = iephi+lephi
      lynm = (ldc+1)**2  
      iynmd = iynm+lynm
      iphitemp = iynmd+lynm
      lphitemp = nq*(2*ldc+1)*2 + 7  
      iphitemp2 = iphitemp+lphitemp
      ifhs = iphitemp2+lphitemp
      ifhder = ifhs + 2*(nterms+1) + 3
      lused = ifhder + 2*(nterms+1) +3
      allocate(w(lused))
c
ccc      call prinf(' ier before mpmpquad0 is *',ier,1)
      call h3dmpmpquad_trunc0(wavek,sc1,x0y0z0,mpole,nterms,nterms1,
     1           sc2,xnynzn,mpolen,nterms2,w(imarray),w(imarray1),
     2           ldc,w(iephi),
     3           radius,xnodes,wts,nquad,nq,w(iynm),
     4           w(iynmd),w(iphitemp),w(iphitemp2),
     7           w(ifhs),w(ifhder),ier)
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h3dmpmpquad_trunc0(wavek,sc1,x0y0z0,mpole,
     $           nterms,nterms1,sc2,
     1           xnynzn,mpolen,nterms2,marray,marray1,ldc,ephi,
     2           radius,xnodes,wts,nquad,nq,ynm,ynmd,
     3           phitemp,phitemp2,fhs,fhder,ier)
C***********************************************************************
C
C     Usage:
C
C           Shift multipole expansion.
C           This is a reasonably fast "point and shoot" version which
C           first rotates the coordinate system, then doing the shifting
C           along the Z-axis, and then rotating back to the original
C           coordinates.
C
C---------------------------------------------------------------------
C     INPUT:
C
C           wavek  = Helmholtz parameter
C           x0y0z0 = center of original multiple expansion
C           xnynzn = center of shifted expansion
C           mpole  = coefficients of original multiple expansion
C           nterms = order of multipole expansion
C           nterms2 = order of shifted expansion
C           sc1     = scaling parameter for mpole expansion
C           sc2     = scaling parameter for shifted expansion
C           radius  = radius of sphere on which shifted expansion is
C                     computed
C           xnodes  = Legendre nodes (precomputed)
C           wts     = Legendre weights (precomputed)
C           nquad   = number of quadrature nodes in theta
C           nq      = max(nquad, 2*ldc+2) = used to allocate
C                     work arrays for both z-shift and rotations.
C
C---------------------------------------------------------------------
C     OUTPUT:
C
C           mpolen = coefficients of shifted expansion
C
C     Work Arrays:
C
C           marray = work array used to hold various intermediate 
c                    expansions.
C           dc     = work array contain the square roots of 
C                    some binomial coefficients.
C           ldc      determines dimension of dc
c                    must exceed max(nterms,nterms2).
C           rd     = work arrays used to store rotation matrices
C                    about Y-axis.
C           ephi    = work array 
C
C***********************************************************************
C
      implicit none
      integer  nterms,nterms1
      integer lw, lused, ier, nq, nquad, nquse,ldc,nterms2
      real *8 x0y0z0(3),xnynzn(3)
      real *8 radius, rshift
      real *8 xnodes(1),wts(1)
      real *8 d,theta,ctheta,phi,sc1,sc2,rvec(3)
      real *8 ynm(0:ldc,0:ldc)
      real *8 ynmd(0:ldc,0:ldc)
      complex *16 phitemp(nq,-ldc:ldc)
      complex *16 phitemp2(nq,-ldc:ldc)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 marray1(0:nterms1,-nterms1:nterms1)
      complex *16 mpolen(0:nterms2,-nterms2:nterms2)
      complex *16 marray(0:ldc,-ldc:ldc)
      complex *16 wavek
c
      complex *16 ephi(-ldc-1:ldc+1),imag
      integer  l,m,jnew,knew
      data imag/(0.0d0,1.0d0)/
C
      rvec(1) = xnynzn(1) - x0y0z0(1)
      rvec(2) = xnynzn(2) - x0y0z0(2)
      rvec(3) = xnynzn(3) - x0y0z0(3)
      call cart2polar(rvec,d,theta,phi)
c
      ephi(1) = exp(imag*phi)
      ephi(0)=1.0d0
      ephi(-1)=dconjg(ephi(1))
c
c----- create array of powers e^(i*m*phi).
c
      do l = 1,ldc
         ephi(l+1) = ephi(l)*ephi(1)
         ephi(-1-l) = dconjg(ephi(l+1))
      enddo
c
c----- a rotation of THETA radians about the Yprime axis after PHI
c      radians about the z-axis.
c      The PHI rotation is carried out on the fly by multiplying 
c      mpole and ephi inside the following loop. 
c
      do l=0,nterms1
         do m=-l,l
            marray1(l,m)=mpole(l,m)*ephi(m)
         enddo
      enddo
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=0.0d0
         enddo
      enddo
c
      if( nterms .ge. 30 ) then 
      call rotviaprojf90(theta,nterms1,nterms1,nterms1,marray1,nterms1,
     1        marray,ldc)
      else
      call rotviarecur3f90(theta,nterms1,nterms1,nterms1,marray1,
     1        nterms1,marray,ldc)
      endif
c
c
c----- shift the mpole expansion from X0Y0Z0 to XNYNZN along
c      the Z-axis.
c
      rshift = d
      call h3dmpmpzshift_fast
     $   (wavek,sc1,marray,ldc,nterms1,sc2,mpolen,
     1           nterms2,nterms2,radius,rshift,xnodes,wts,nquad,
     2           ynm,phitemp,fhs,fhder,ier)
c
c
c     Reverse THETA rotation.
c     I.e. rotation of -THETA radians about Yprime axis.
c
      if( nterms2 .ge. 30 ) then
      call rotviaprojf90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      else
      call rotviarecur3f90(-theta,nterms2,nterms2,nterms2,mpolen,
     1        nterms2,marray,ldc)
      endif
c
c
c----- rotate back PHI radians about the Z-axis in the above system.
c
      do l=0,nterms2
         do m=-l,l
            mpolen(l,m)=ephi(-m)*marray(l,m)
         enddo
      enddo
      return
      end
c
c
c
c
c***********************************************************************
      subroutine h3dmpmpzshift_fast
     $     (zk,scale,mpole,lmp,nterms,scale2,mpolen,
     1      lmpn,nterms2,radius,zshift,xnodes,wts,nquad,ynm,
     2      phitemp,fhs,fhder,ier)
c***********************************************************************
c
c     This subroutine converts a multipole expansion centered at the 
c     origin to a multipole expansion centered at (0,0,zhift).
c     The expansion is rescaled to that of the shifted expansion.
c
c INPUT:
c
c     zk       : Helmholtz coefficient
c     scale    : scale parameter for mpole
c     mpole    : coefficients of original multipole exp.
c     lmp      : leading dim of mpole (may be a work array)
c     nterms   : number of terms in the orig. expansion
c
c     scale2   : scale parameter for new expansion (mpolen)
c     lmpn     : leading dim of shifted (may be work array)
c     nterms2  : number of terms in output expansion
c     radius   : radius of sphere on which mpole is evaluated
c                           in projeciton step
c     zshift   : shifting distance along z-axis
c                              (always assumed positive)
C     xnodes   : Legendre nodes (precomputed)
C     wts      : Legendre weights (precomputed)
C     nquad    : number of quadrature nodes in theta direction
c
c OUTPUT:
c
c     mpolen  (complex *16)  : coefficients of shifted exp.
c
c***********************************************************************
      implicit real *8 (a-h,o-z)
      integer nterms,nterms2,nquad,ier,lmp,lmpn,ldc,iynm,lynm
      real *8   zshift,scale,scale2,radius
      real *8   xnodes(1),wts(1)
      real *8   ynm(0:nterms,0:nterms)
      complex *16 phitemp(nquad,-nterms:nterms)
      complex *16 fhs(0:nterms)
      complex *16 fhder(0:nterms)
      complex *16 mpole(0:lmp,-lmp:lmp),zk
      complex *16 mpolen(0:lmpn,-lmpn:lmpn)
c
c     local allocated workspace arrays - no more passed workspace
c
      real *8, allocatable :: w(:)
c
      integer l,m,jnew,knew
C
C----- shift along z-axis by evaluating field on target sphere and
C     projecting onto spherical harmonics and scaling by j_n(kR).
C
C    OPTIMIZATION NOTES:
C
C    Suppose you are shifting from a very small sphere to the center
C    of a very large sphere (nterms2 >> nterms).
C    Then, ALONG THE Z-AXIS, the number of azimuthal modes that
C    need to be computed is only nterms (not nterms2). 
C    Subroutines h3dmpevalspherenm, h3dprojlocnmsep allow for this.
C    The final step of the point and shoot algorithm then distributes
C    these nterms (azimuthal) modes to nterms2 (azimuthal) in the
C    "laboratory frame".
C
C    cost is (nterms^2 x nterms2) rather than (nterms x nterms2^2)
C

        ldc = max(nterms,nterms2)
        irat1=1
        lrat1=(ldc+1)**2
        irat2=irat1+lrat1
        lrat2=(ldc+1)**2
        lused=irat2+lrat2
        allocate(w(lused))
c
ccc      call prinf(' allocated in shift fast *',lused,1)
      call h3dmpevalspherenm_fast(mpole,zk,scale,
     1     zshift,radius,nterms,lmp,ynm,
     2     phitemp,nquad,xnodes,fhs,fhder,w(irat1),w(irat2))
      call h3dprojlocnmsep_fast
     $   (nterms2,lmpn,nquad,nterms,xnodes,wts,
     1     phitemp,mpolen,ynm,w(irat1),w(irat2))
      call h3drescalemp(nterms2,lmpn,mpolen,radius,zk,
     1               scale2,fhs,fhder)
      return
      end
C
