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
c      This file contains the accelerated subroutines for 
c      forming and evaluating multipole expansions.
c
c
c**********************************************************************
      subroutine h3dmpevalall_trunc(zk,rscale,center,mpole,
     $     nterms,nterms1,ztarg,nt,
     1     ifpot,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c
c     pot =  sum sum  mpole(n,m) h_n(k r) Y_nm(theta,phi)
c             n   m
c     fld = -gradient(pot) if iffld = 0.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     nterms :    order of the multipole expansion
c     lw     :    length of workspace (greater than 
c                    (nterms+1)**2 + 8*(nterms+1) + 100
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter (see formmp1h3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     w      :    workspace
c     ztarg  :    target location
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     lused  :    amount of work space w that is actually used
c     ier    :    error return code
c		      ier=0  successful execution
c		      ier=8  insufficient work space 
c		      ier=16 ztarg too close to center
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 center(3),ztarg(3,1)
      complex *16 zk,pot(1),fld(3,1),pot0,fld0(3)
      complex *16 mpole(0:nterms,-nterms:nterms)
      real *8, allocatable :: w(:)
c
      ier=0
c
c     Carve up workspace:
c
c     for Ynm and Ynm'
c
      lpp=(nterms+1)**2+5
      ipp=1
      ippd = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
      if (iffld.eq.1) then
         iephi=ippd+lpp
      else
         iephi = ippd+2
      endif
      lephi=2*(2*nterms+3)+5
c
c     workspace for spherical Hankel functions
c
      ifhs=iephi+lephi
      lfhs=2*(nterms+1)+5
c
c     workspace for derivatives of spherical Hankel functions
c
      ifhder=ifhs+lfhs
      lfhder=2*(nterms+1)+5
c
      lused=ifhder+lfhder
      allocate(w(lused))
c
      do i=1,nt
      call h3dmpeval_trunc0(jer,zk,rscale,center,mpole,
     $   nterms,nterms1,ztarg(1,i),
     1   pot0,iffld,fld0,w(ipp),w(ippd),w(iephi),w(ifhs),w(ifhder),
     $   wlege,nlege)
      if( ifpot .eq. 1 ) pot(i)=pot(i)+pot0
      if( iffld .eq. 1 ) then
        fld(1,i)=fld(1,i)+fld0(1)
        fld(2,i)=fld(2,i)+fld0(2)
        fld(3,i)=fld(3,i)+fld0(3)
      endif
      enddo
c
ccc      if (jer.ne.0) ier=16
c
      return
      end
c
c
c**********************************************************************
      subroutine h3dmpeval_trunc(zk,rscale,center,mpole,
     $     nterms,nterms1,ztarg,
     1     pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates the potential and gradient of the 
c     potential due to an outgoing partial wave expansion.
c
c     pot =  sum sum  mpole(n,m) h_n(k r) Y_nm(theta,phi)
c             n   m
c     fld = -gradient(pot) if iffld = 0.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     iffld  :   flag controlling evaluation of gradient:
c                   iffld = 0, do not compute gradient.
c                   iffld = 1, compute gradient.
c     nterms :    order of the multipole expansion
c     lw     :    length of workspace (greater than 
c                    (nterms+1)**2 + 8*(nterms+1) + 100
c     zk     :    Helmholtz parameter
c     rscale :    scaling parameter (see formmp1h3d)
c     center :    expansion center
c     mpole  :    multipole expansion in 2d matrix format
c     w      :    workspace
c     ztarg  :    target location
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     lused  :    amount of work space w that is actually used
c     ier    :    error return code
c		      ier=0  successful execution
c		      ier=8  insufficient work space 
c		      ier=16 ztarg too close to center
c     pot    :    potential at ztarg
c     fld    :    gradient at ztarg (if requested)
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 center(3),ztarg(3)
      complex *16 zk,pot,fld(3)
      complex *16 mpole(0:nterms,-nterms:nterms)
      real *8, allocatable :: w(:)
c
      ier=0
c
c     Carve up workspace:
c
c     for Ynm and Ynm'
c
      lpp=(nterms+1)**2+5
      ipp=1
      ippd = ipp+lpp
c
c     workspace for azimuthal argument (ephi)
c
      if (iffld.eq.1) then
         iephi=ippd+lpp
      else
         iephi = ippd+2
      endif
      lephi=2*(2*nterms+3)+5
c
c     workspace for spherical Hankel functions
c
      ifhs=iephi+lephi
      lfhs=2*(nterms+1)+5
c
c     workspace for derivatives of spherical Hankel functions
c
      ifhder=ifhs+lfhs
      lfhder=2*(nterms+1)+5
c
      lused=ifhder+lfhder
      allocate(w(lused))
c
      call h3dmpeval_trunc0(jer,zk,rscale,center,mpole,
     $   nterms,nterms1,ztarg,
     1	   pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifhs),w(ifhder),
     $   wlege,nlege)
      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dmpeval_trunc0(ier,zk,rscale,center,mpole,
     $     nterms,nterms1,
     1		ztarg,pot,iffld,fld,ynm,ynmd,ephi,fhs,fhder,
     $     wlege,nlege)
c**********************************************************************
c
c     See h3dmpeval for comments.
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 center(3),ztarg(3),zdiff(3)
      real *8 ynm(0:nterms,0:nterms)
      real *8 ynmd(0:nterms,0:nterms)
      complex *16 zk,pot,fld(3),ephi1,ur,utheta,uphi,ux,uy,uz
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
      complex *16 fhs(0:nterms),fhder(0:nterms)
c
      complex *16 eye
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum,z
c
      data eye/(0.0d0,1.0d0)/
      data thresh/1.0d-15/
c
      ier=0
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
      if (abs(zk*r).lt.thresh) then
         ier=8
         return
      endif
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=done
      ephi(1)=ephi1
      cphi = dreal(ephi1)
      sphi = dimag(ephi1)
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
      if (iffld.eq.1) then
         rx = stheta*cphi
         thetax = ctheta*cphi/r
         phix = -sphi/r
         ry = stheta*sphi
         thetay = ctheta*sphi/r
         phiy = cphi/r
         rz = ctheta
         thetaz = -stheta/r
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
ccc         call ylgndr2s(nterms,ctheta,ynm,ynmd)
         call ylgndr2sfw(nterms1,ctheta,ynm,ynmd,wlege,nlege)
      else
ccc         call ylgndr(nterms,ctheta,ynm)
         call ylgndrfw(nterms1,ctheta,ynm,wlege,nlege)
      endif
c
c     get the Hankel functions:
c
      ifder=iffld
      z=zk*r
      call h3dall(nterms1,z,rscale,fhs,ifder,fhder)
c
c     initialize computed values and 
c     scale derivatives of Hankel functions so that they are
c     derivatives with respect to r.
c
      if (iffld.eq.1) then
         do i=0,nterms1
            fhder(i)=fhder(i)*zk
         enddo
         ur = mpole(0,0)*fhder(0)
         utheta = 0.0d0
         uphi = 0.0d0
      endif
      pot=mpole(0,0)*fhs(0)
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms1
	    pot=pot+mpole(n,0)*fhs(n)*ynm(n,0)
	    ur = ur + fhder(n)*ynm(n,0)*mpole(n,0)
	    utheta = utheta -mpole(n,0)*fhs(n)*ynmd(n,0)*stheta
	    do m=1,n
	       ztmp1=fhs(n)*ynm(n,m)*stheta
	       ztmp2 = mpole(n,m)*ephi(m) 
	       ztmp3 = mpole(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + fhder(n)*ynm(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fhs(n)*ynmd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fhs(n)*ynm(n,m)*ztmpsum
            enddo
         enddo
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms1
	    pot=pot+mpole(n,0)*fhs(n)*ynm(n,0)
	    do m=1,n
	       ztmp1=fhs(n)*ynm(n,m)
	       ztmp2 = mpole(n,m)*ephi(m) + mpole(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine h3dtaevalall_trunc(wavek,rscale,center,locexp,nterms,
     $     nterms1,ztarg,nt,
     $     ifpot,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates a j-expansion centered at CENTER
c     at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) j_n(k r) Y_nm(theta,phi)
c             n   m
c
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek      : the Helmholtz coefficient
c     rscale     : scaling parameter used in forming expansion
c                                   (see h3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the j-expansion
c     nterms     : order of the h-expansion
c     ztarg(3)   : target vector
c     iffld      : flag for gradient computation
c		                    iffld=0  - gradient is not computed
c		                    iffld=1  - gradient is computed
c     w          : workspace (real)
c     lw         : workspace length
c                            at least (nterms+1)**2 +
c                                     6*(nterms+1) + 7000
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier        : error return code
c		      ier=0	returned successfully
c		      ier=8 insuffficient workspace 
c		      ier=16 insufficient memory 
c                            in subroutine "jfuns3d"
c     pot        : potential at ztarg(3)
c     fld(3)     : gradient at ztarg (if requested)
c     lused      : amount of work space "w" used
c
c     NOTE: Parameter lwfjs is set to nterms+1000
c           Should be sufficient for any Helmholtz parameter
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer lwfjs
      real *8 center(3),ztarg(3,1)
      real *8, allocatable :: w(:)
      complex *16 wavek,pot(1),fld(3,1),pot0,fld0(3)
      complex *16 locexp(0:nterms,-nterms:nterms)
c
c ... Assigning work spaces for various temporary arrays:
c
      ier=0
c
      lwfjs=nterms+1000
      ipp=1
      lpp=(nterms+1)**2+3
      ippd  = ipp+lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      iiscale=iephi+lephi
      liscale=(lwfjs+1)+3
c
      ifjs=iiscale+liscale
      lfjs=2*(lwfjs+1)+3
c
      ifjder=ifjs+lfjs
      lfjder=2*(nterms+1)+3
c
      lused=ifjder+lfjder
      allocate(w(lused))
c
      do i=1,nt
      call h3dtaeval_trunc0(jer,wavek,rscale,center,locexp,
     $   nterms,nterms1,ztarg(1,i),
     1	     pot0,iffld,fld0,w(ipp),w(ippd),w(iephi),w(ifjs),
     2       w(ifjder),lwfjs,w(iiscale),wlege,nlege)
      if( ifpot .eq. 1 ) pot(i)=pot(i)+pot0
      if( iffld .eq. 1 ) then
        fld(1,i)=fld(1,i)+fld0(1)
        fld(2,i)=fld(2,i)+fld0(2)
        fld(3,i)=fld(3,i)+fld0(3)
      endif
      enddo
ccc      if (jer.ne.0) ier=16
c
      return
      end
c
c
c**********************************************************************
      subroutine h3dtaeval_trunc(wavek,rscale,center,locexp,nterms,
     $     nterms1,
     1		ztarg,pot,iffld,fld,wlege,nlege,ier)
c**********************************************************************
c
c     This subroutine evaluates a j-expansion centered at CENTER
c     at the target point ZTARG. 
c
c     pot =  sum sum  locexp(n,m) j_n(k r) Y_nm(theta,phi)
c             n   m
c
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek      : the Helmholtz coefficient
c     rscale     : scaling parameter used in forming expansion
c                                   (see h3dformmp1)
c     center     : coordinates of the expansion center
c     locexp     : coeffs of the j-expansion
c     nterms     : order of the h-expansion
c     ztarg(3)   : target vector
c     iffld      : flag for gradient computation
c		                    iffld=0  - gradient is not computed
c		                    iffld=1  - gradient is computed
c     w          : workspace (real)
c     lw         : workspace length
c                            at least (nterms+1)**2 +
c                                     6*(nterms+1) + 7000
c
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier        : error return code
c		      ier=0	returned successfully
c		      ier=8 insuffficient workspace 
c		      ier=16 insufficient memory 
c                            in subroutine "jfuns3d"
c     pot        : potential at ztarg(3)
c     fld(3)     : gradient at ztarg (if requested)
c     lused      : amount of work space "w" used
c
c     NOTE: Parameter lwfjs is set to nterms+1000
c           Should be sufficient for any Helmholtz parameter
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer lwfjs
      real *8 center(3),ztarg(3)
      real *8, allocatable :: w(:)
      complex *16 wavek,pot,fld(3)
      complex *16 locexp(0:nterms,-nterms:nterms)
c
c ... Assigning work spaces for various temporary arrays:
c
      ier=0
c
      lwfjs=nterms+1000
      ipp=1
      lpp=(nterms+1)**2+3
      ippd  = ipp+lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      iiscale=iephi+lephi
      liscale=(lwfjs+1)+3
c
      ifjs=iiscale+liscale
      lfjs=2*(lwfjs+1)+3
c
      ifjder=ifjs+lfjs
      lfjder=2*(nterms+1)+3
c
      lused=ifjder+lfjder
      allocate(w(lused))
c
      call h3dtaeval_trunc0(jer,wavek,rscale,center,locexp,
     $   nterms,nterms1,ztarg,
     1	     pot,iffld,fld,w(ipp),w(ippd),w(iephi),w(ifjs),
     2       w(ifjder),lwfjs,w(iiscale),wlege,nlege)
      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dtaeval_trunc0(ier,wavek,rscale,center,locexp,
     $     nterms,nterms1,ztarg,
     $     pot,iffld,fld,pp,ppd,ephi,fjs,fjder,lwfjs,iscale,
     $     wlege,nlege)
c**********************************************************************
c
c     See h3dtaeval for comments.
c     (pp and ppd are storage arrays for Ynm and Ynm')
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer iscale(0:1)
      real *8 center(3),ztarg(3),zdiff(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      complex *16 wavek,pot,fld(3),ephi1,ephi1inv
      complex *16 locexp(0:nterms,-nterms:nterms)
      complex *16 ephi(-nterms-1:nterms+1)
      complex *16 fjsuse,fjs(0:1),fjder(0:1)
c
      complex *16 eye,ur,utheta,uphi
      complex *16 ztmp,z
      complex *16 ztmp1,ztmp2,ztmp3,ztmpsum
      complex *16 ux,uy,uz
c
      data eye/(0.0d0,1.0d0)/
c
      ier=0
      done=1.0d0
c
      zdiff(1)=ztarg(1)-center(1)
      zdiff(2)=ztarg(2)-center(2)
      zdiff(3)=ztarg(3)-center(3)
c
c     Convert to spherical coordinates
c
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute e^{eye*m*phi} array.
c
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     For the n=0 mode, it is not relevant. For n>0 modes,
c     we use the recurrence relation 
c
c     (2n+1)fjs_n(kr)/(kr) = fjs(n+1)*rscale + fjs(n-1)/rscale
c
c     to avoid division by r. The variable fjsuse is set to fjs(n)/r:
c
c           fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
c	    fjsuse = wavek*fjsuse/(2*n+1.0d0)
c
c     
c
      if (iffld.eq.1) then
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
      endif
c
c     get the associated Legendre functions:
c
      if (iffld.eq.1) then
ccc         call ylgndr2s(nterms,ctheta,pp,ppd)
         call ylgndr2sfw(nterms1,ctheta,pp,ppd,wlege,nlege)
      else
ccc         call ylgndr(nterms,ctheta,pp)
         call ylgndrfw(nterms1,ctheta,pp,wlege,nlege)
      endif
c
c     get the spherical Bessel functions and their derivatives.
c
      ifder=iffld
      z=wavek*r
      call jfuns3d(jer,nterms1,z,rscale,fjs,ifder,fjder,
     1	      lwfjs,iscale,ntop)
      if (jer.ne.0) then
         ier=8
         return
      endif
c
c     scale derivatives of Bessel functions so that they are
c     derivatives with respect to r.
c
c
      pot=locexp(0,0)*fjs(0)
      if (iffld.eq.1) then
         do i=0,nterms1
            fjder(i)=fjder(i)*wavek
         enddo
         ur = locexp(0,0)*fjder(0)
         utheta = 0.0d0
         uphi = 0.0d0
      endif
c
c     compute the potential and the field:
c
      if (iffld.eq.1) then
         do n=1,nterms1
            pot=pot+locexp(n,0)*fjs(n)*pp(n,0)
	    ur = ur + fjder(n)*pp(n,0)*locexp(n,0)
	    fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
	    fjsuse = wavek*fjsuse/(2*n+1.0d0)
	    utheta = utheta -locexp(n,0)*fjsuse*ppd(n,0)*stheta
	    do m=1,n
	       ztmp1=fjs(n)*pp(n,m)*stheta
	       ztmp2 = locexp(n,m)*ephi(m) 
	       ztmp3 = locexp(n,-m)*ephi(-m)
	       ztmpsum = ztmp2+ztmp3
	       pot=pot+ztmp1*ztmpsum
	       ur = ur + fjder(n)*pp(n,m)*stheta*ztmpsum
	       utheta = utheta -ztmpsum*fjsuse*ppd(n,m)
	       ztmpsum = eye*m*(ztmp2 - ztmp3)
	       uphi = uphi + fjsuse*pp(n,m)*ztmpsum
            enddo
         enddo
ccc	 call prin2(' ur is *',ur,2)
ccc	 call prin2(' utheta is *',utheta,2)
ccc	 call prin2(' uphi is *',uphi,2)
	 ux = ur*rx + utheta*thetax + uphi*phix
	 uy = ur*ry + utheta*thetay + uphi*phiy
	 uz = ur*rz + utheta*thetaz + uphi*phiz
	 fld(1) = -ux
	 fld(2) = -uy
	 fld(3) = -uz
      else
         do n=1,nterms1
	    pot=pot+locexp(n,0)*fjs(n)*pp(n,0)
	    do m=1,n
	       ztmp1=fjs(n)*pp(n,m)
	       ztmp2 = locexp(n,m)*ephi(m)+locexp(n,-m)*ephi(-m)
	       pot=pot+ztmp1*ztmp2
            enddo
         enddo
      endif
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h3dformmp_trunc(ier,zk,scale,sources,charge,ns,center,
     1                  nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     zk              : Helmholtz parameter 
C     scale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     w               : workspace
C     lw              : length of  workspace
C                             at least (nterms+1)**2 +
C                                   8*(nterms+1) + 100
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		                   ier=0  returned successfully
c		                   ier=8  insufficient memory 
c		                   ier=16 insufficient memory 
c                                         in subroutine "jfuns3d"
c                                         called in h3dformmp1
c    
c
c     mpole           : coeffs of the h-expansion
c                  
c     lused           : amount of work space "w" used
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ns,i,l,m, ier, ier1, lused
      real *8 center(1),sources(3,ns)
      real *8 scale
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,zk,charge(1)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      ier = 0
      do i = 1, ns
         call h3dformmp_trunc1
     $   (ier1,zk,scale,sources(1,i),charge(i),center,
     1        nterms,nterms1,mpole,wlege,nlege)
      enddo
      if (ier1.ne.0) ier = ier1
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = mpole(l,m)*eye*zk
         enddo
      enddo
C
      return
      end
C
C***********************************************************************
      subroutine h3dformmp_add_trunc
     $     (ier,zk,scale,sources,charge,ns,center,
     1                  nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS sources 
C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     zk              : Helmholtz parameter 
C     scale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     charge(ns)      : source strengths
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
C     w               : workspace
C     lw              : length of  workspace
C                             at least (nterms+1)**2 +
C                                   8*(nterms+1) + 100
C
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		                   ier=0  returned successfully
c		                   ier=8  insufficient memory 
c		                   ier=16 insufficient memory 
c                                         in subroutine "jfuns3d"
c                                         called in h3dformmp1
c    
c
c     mpole           : coeffs of the h-expansion
c                  
c     lused           : amount of work space "w" used
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ns,i,l,m, ier, ier1, lused
      real *8 center(1),sources(3,ns)
      real *8 scale
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,zk,charge(1)
      complex *16, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
        allocate( mptemp(0:nterms,-nterms:nterms) )

        do l = 0,nterms
          do m=-l,l
             mptemp(l,m) = 0
          enddo
        enddo

        call h3dformmp_trunc
     $     (ier,zk,scale,sources,charge,ns,center,
     1     nterms,nterms1,mptemp,wlege,nlege)

      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
         enddo
      enddo
C
      return
      end
C
c**********************************************************************
      subroutine h3dformmp_trunc1(ier,zk,rscale,source,charge,center,
     1		nterms,nterms1,mpole,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the h-expansion about CENTER
c     due to a charge located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk      : the Helmholtz coefficient
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     charge  : complex charge strength
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c     w       : workspace
c     lw      : workspace length
c                         at least (nterms+1)**2 +
c                            6*(nterms+1) + 7000
c
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      ier=8	insufficient memory 
c		      ier=16 insufficient memory 
c                            in subroutine "jfuns3d"
c                            called in h3dformmp0
c                            
c     mpole   : coeffs of the h-expansion
c            
c     lused   : amount of work space "w" used
c
c     NOTE: Parameter lwfjs is set to nterms+1000
c           Should be sufficient for any Helmholtz parameter
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer lwfjs
      real *8 source(3),center(3)
      real *8, allocatable :: w(:)
      complex *16 zk,mpole(0:nterms,-nterms:nterms)
      complex *16 charge
c
c ... Assign work spaces:
c
      ier=0
c
      lwfjs=nterms+1000
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifjder=iephi+lephi
      lfjder=2*(nterms+1)+7
c
      ifjs=ifjder+lfjder
      lfjs=2*(lwfjs+1)+7
c
      iiscale=ifjs+lfjs
      liscale=(lwfjs+1)+7
c
      lused=iiscale+liscale
      allocate(w(lused))
c
      call h3dformmp_trunc0(jer,zk,rscale,source,charge,center,
     $   nterms,nterms1,
     1		mpole,w(ipp),w(ippd),w(iephi),w(ifjs),lwfjs,
     2          w(iiscale),w(ifjder),wlege,nlege)
      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dformmp_trunc0(ier,zk,rscale,source,charge,center,
     1		nterms,nterms1,
     $     mpole,pp,ppd,ephi,fjs,lwfjs,iscale,fjder,wlege,nlege)
c**********************************************************************
c
c     See h3dformmp1 for comments.
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer iscale(0:1)
      real *8 source(3),center(3),zdiff(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      complex *16 zk,mpole(0:nterms,-nterms:nterms)
      complex *16 charge
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fjs(0:1),ztmp,fjder(0:1),z
      data thresh/1.0d-15/
c
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     get the associated Legendre functions:
c
ccc      call ylgndr(nterms,ctheta,pp)
      call ylgndrfw(nterms1,ctheta,pp,wlege,nlege)
ccc      call ylgndr2s(nterms,ctheta,pp,ppd)
ccc      call prinf(' after ylgndr with nterms = *',nterms,1)
ccc      call prinm2(pp,nterms)
c
c     get Bessel functions:
c
      ifder=0
      z=zk*r
      call jfuns3d(jer,nterms1,z,rscale,fjs,ifder,fjder,
     1	      lwfjs,iscale,ntop)
c
c
c     multiply all jn by charge strength.
c
      do n = 0,nterms1
         fjs(n) = fjs(n)*charge
      enddo
      if (jer.ne.0) then
	 ier=16
	 return
      endif
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|S|) Ylm*(S) h_n(k|T|)Ylm(T)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of after all source contributions
c     have been included in the calling subroutine h3dformmp.
c
      mpole(0,0)= mpole(0,0) + fjs(0)
      do n=1,nterms1
         dtmp=pp(n,0)
         mpole(n,0)= mpole(n,0) + dtmp*fjs(n)
         do m=1,n
cc            ztmp=stheta*pp(n,m)*fjs(n)
            ztmp=pp(n,m)*fjs(n)
            mpole(n, m)= mpole(n, m) + ztmp*dconjg(ephi(m))
            mpole(n,-m)= mpole(n,-m) + ztmp*dconjg(ephi(-m))
         enddo
      enddo
c

c
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine h3dformta_trunc
     $     (ier,zk,rscale,sources,charge,ns,center,
     1     nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta1/h3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     zk       : Helmholtz coefficient
c     rscale   : scaling parameter
c                     should be less than one in magnitude.
c                     Needed for low frequency regime only
c                     with rsclale abs(wavek) recommended.
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c     w         : workspace
c     lw        : workspace length
c                           at least (nterms+1)**2 +
c                              6*(nterms+1) + 7000
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c	 	  ier=4  d is out of range in h3dall
c
c     locexp    : coeffs for the j-expansion
c     lused     : amount of work space "w" used
c
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(3,ns),center(3)
      complex *16 zk,locexp(0:nterms,-nterms:nterms), charge(ns)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call h3dformta_trunc1(ier,zk,rscale,sources(1,i),charge(i),
     1		center,nterms,nterms1,locexp,wlege,nlege)
      enddo
c
c     scale by (i*k)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m)*eye*zk
         enddo
      enddo
C
      return
      end
c
c
c**********************************************************************
      subroutine h3dformta_add_trunc
     $     (ier,zk,rscale,sources,charge,ns,center,
     1     nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS sources at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta1/h3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     zk       : Helmholtz coefficient
c     rscale   : scaling parameter
c                     should be less than one in magnitude.
c                     Needed for low frequency regime only
c                     with rsclale abs(wavek) recommended.
c     sources   : coordinates of the sources
c     charge    : charge strengths
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c     w         : workspace
c     lw        : workspace length
c                           at least (nterms+1)**2 +
c                              6*(nterms+1) + 7000
c
c
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c	 	  ier=4  d is out of range in h3dall
c
c     locexp    : coeffs for the j-expansion
c     lused     : amount of work space "w" used
c
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(3,ns),center(3)
      complex *16 zk,locexp(0:nterms,-nterms:nterms), charge(ns)
      complex *16 eye
      complex *16, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c

        allocate( mptemp(0:nterms,-nterms:nterms) )
c
        do l = 0,nterms
          do m=-l,l
             mptemp(l,m) = 0
          enddo
        enddo
c
        call h3dformta_trunc
     $     (ier,zk,rscale,sources,charge,ns,center,
     1     nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m)+mptemp(l,m)
         enddo
      enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine h3dformta_trunc1(ier,wavek,rscale,source,charge,center,
     &		nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single charge located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta0 below.
c
c---------------------------------------------------------------------
c INPUT:
c
c     wavek     : the Helmholtz coefficient
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     charge    : coordinates of the source
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c     w         : workspace
c     lw        : workspace length
c                           at least (nterms+1)**2 +
c                                    6*(nterms+1) + 1000
c
c---------------------------------------------------------------------
c OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   ier=2 insufficient memory in workspace w
c	 	   ier=4 d is out of range in h3dall
c     locexp : coefficients of the local expansion
c     lused  : amount of work space "w" used
c
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(3),center(3)
      real *8, allocatable :: w(:)
      complex *16 wavek,locexp(0:nterms,-nterms:nterms), charge
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      iephi=ipp+lpp
      lephi=2*(2*nterms+1)+7
c
      ifhs=iephi+lephi
      lfhs=2*(nterms+1)+7
c
      lused=ifhs+lfhs
      allocate(w(lused))
c
      call h3dformta_trunc0(jer,wavek,rscale,source,charge,center,
     &		nterms,nterms1,locexp,w(ipp),w(iephi),w(ifhs),
     $   wlege,nlege)
      if (jer.ne.0) ier=4
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dformta_trunc0(ier,wavek,rscale,source,charge,
     &		center,nterms,nterms1,locexp,pp,ephi,fhs,wlege,nlege)
c**********************************************************************
c
c     See h3dformta/h3dformta1 for comments
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(3),center(3),zdiff(3)
      real *8 pp(0:nterms,0:nterms)
      complex *16 wavek,locexp(0:nterms,-nterms:nterms), charge
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fhs(0:nterms),ztmp,fhder(0:1),z
      data thresh/1.0d-15/
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      do i=2,nterms
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
c
c     get the Ynm
c
ccc      call ylgndr(nterms,ctheta,pp)
      call ylgndrfw(nterms,ctheta,pp,wlege,nlege)
c
c     compute Hankel functions and scale them by charge strength.
c
      ifder=0
      z=wavek*r
      if (abs(z).lt.thresh) then
         ier = 4
         return
      endif
      call h3dall(nterms1,z,rscale,fhs,ifder,fhder)
      do n = 0, nterms1
         fhs(n) = fhs(n)*charge
      enddo
c
c     Compute contributions to locexp
c
      locexp(0,0)=locexp(0,0) + fhs(0)
      do n=1,nterms1
         locexp(n,0)=locexp(n,0) + pp(n,0)*fhs(n)
         do m=1,n
            ztmp=pp(n,m)*fhs(n)
	    locexp(n,m)=locexp(n,m) + ztmp*ephi(-m)
	    locexp(n,-m)=locexp(n,-m) + ztmp*ephi(m)
         enddo
      enddo
      return
      end
c
c
c
C***********************************************************************
      subroutine h3dformmp_dp_trunc
     $     (ier,zk,scale,sources,dipstr,dipvec,ns,
     1                  center,nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS 
c     dipole sources C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     zk              : Helmholtz parameter 
C     scale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipstr(ns)      : source strengths
C     dipvec(3,ns)    : dipole vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		          ier=0  returned successfully
c		          ier=8  insufficient memory 
c                         ier=16 insufficient memory 
c                            in subroutine "jfuns3d" called
c                            in h3dformmp1_dp
c
c     mpole           : coeffs of the h-expansion
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ns,i,l,m, ier, lused
      real *8 center(1),sources(3,ns)
      real *8 dipvec(3,ns)
      real *8 scale
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,zk,dipstr(1)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1, ns
         call h3dformmp1_dp_trunc(ier,zk,scale,sources(1,i),dipstr(i),
     1        dipvec(1,i),center,nterms,nterms1,mpole,wlege,nlege)
      enddo
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = mpole(l,m)*eye*zk
         enddo
      enddo
C
      return
      end
C
C***********************************************************************
      subroutine h3dformmp_dp_add_trunc
     $     (ier,zk,scale,sources,dipstr,dipvec,ns,
     1                  center,nterms,nterms1,mpole,wlege,nlege)
C***********************************************************************
C
C     Constructs multipole (h) expansion about CENTER due to NS 
c     dipole sources C     located at SOURCES(3,*).
C
c-----------------------------------------------------------------------
C     INPUT:
c
C     zk              : Helmholtz parameter 
C     scale           : the scaling factor.
C     sources(3,ns)   : coordinates of sources
C     dipstr(ns)      : source strengths
C     dipvec(3,ns)    : dipole vector direction 
C     ns              : number of sources
C     center(3)       : epxansion center
C     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
C     OUTPUT:
C
c     ier             : error return code
c		          ier=0  returned successfully
c		          ier=8  insufficient memory 
c                         ier=16 insufficient memory 
c                            in subroutine "jfuns3d" called
c                            in h3dformmp1_dp
c
c     mpole           : coeffs of the h-expansion
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms,ns,i,l,m, ier, lused
      real *8 center(1),sources(3,ns)
      real *8 dipvec(3,ns)
      real *8 scale
      complex *16 mpole(0:nterms,-nterms:nterms)
      complex *16 eye,zk,dipstr(1)
      complex *16, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
C
C----- set mpole to zero
C
        allocate( mptemp(0:nterms,-nterms:nterms) )
c
      do l = 0,nterms
         do m=-l,l
            mptemp(l,m) = 0.0d0
         enddo
      enddo
c
      call h3dformmp_dp_trunc
     $     (ier,zk,scale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mptemp,wlege,nlege)
c
      do l = 0,nterms
         do m=-l,l
            mpole(l,m) = mpole(l,m)+mptemp(l,m)
         enddo
      enddo
C
      return
      end
C
c**********************************************************************
      subroutine h3dformmp1_dp_trunc
     $     (ier,zk,rscale,source,dipstr,dipvec,
     1		center,nterms,nterms1,mpole,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the h-expansion about CENTER
c     due to a dipole located at the point SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformmp0 below.
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     zk      : the Helmholtz coefficient
c     rscale  : scaling parameter
c     source  : coordinates of the charge
c     dipstr  : complex dipole strength
c     dipvec  : dipole direction vector
c     center  : coordinates of the expansion center
c     nterms  : order of the h-expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     ier     : error return code
c		      ier=0 returned successfully
c		      ier=8	insufficient memory 
c                     ier=16 insufficient memory 
c                         in subroutine "jfuns3d" called
c                         in h3dformmp1_dp
c     mpole   : coeffs of the h-expansion
c            
c     lused   : amount of work space "w" used
c
c
c     NOTE: Parameter lwfjs is set to nterms+1000
c           Should be sufficient for any Helmholtz parameter
c
c-----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer lwfjs
      real *8 source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 dipvec(3)
      complex *16 zk,mpole(0:nterms,-nterms:nterms)
      complex *16 dipstr
c
c ... Assign work spaces:
c
      ier=0
c
      lwfjs=nterms+1000
      ipp=1
      lpp=(nterms+1)**2+7
      ippd = ipp + lpp
c
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifjder=iephi+lephi
      lfjder=2*(nterms+1)+7
c
      ifjs=ifjder+lfjder
      lfjs=2*(lwfjs+1)+7
c
      iiscale=ifjs+lfjs
      liscale=(lwfjs+1)+7
c
      lused=iiscale+liscale
      allocate(w(lused))
c
      call h3dformmp0_dp_trunc(jer,zk,rscale,source,dipstr,dipvec,
     1		center,nterms,nterms1,mpole,w(ipp),w(ippd),w(iephi),
     2          w(ifjs),lwfjs,w(iiscale),w(ifjder),wlege,nlege)
      if (jer.ne.0) ier=16
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dformmp0_dp_trunc
     $     (ier,wavek,rscale,source,dipstr,dipvec,
     1     center,nterms,nterms1,
     $     mpole,pp,ppd,ephi,fjs,lwfjs,iscale,fjder,wlege,nlege)
c**********************************************************************
c
c     See h3dformmp1_dp for comments.
c
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer iscale(0:1)
      real *8 source(3),center(3),zdiff(3)
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      complex *16 wavek,mpole(0:nterms,-nterms:nterms)
      complex *16 dipstr
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fjs(0:1),ztmp,fjder(0:1),z
      complex *16 fjsuse,ux,uy,uz,ur,utheta,uphi,zzz
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(1.0d0-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     compute exp(eye*m*phi) array
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=dconjg(ephi1)
      do i=2,nterms+1
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi(-1)
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
c     In thetax, thetaty, phix, phiy we leave out the 1/r factors in the 
c     change of variables to avoid blow-up at the origin.
c     For the n=0 mode, it is not relevant. For n>0 modes,
c     we use the recurrence relation 
c
c     (2n+1)fjs_n(kr)/(kr) = fjs(n+1)*rscale + fjs(n-1)/rscale
c
c     to avoid division by r. The variable fjsuse is set to fjs(n)/r:
c
c           fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
c	    fjsuse = wavek*fjsuse/(2*n+1.0d0)
c
c     
c
         rx = stheta*cphi
ccc         thetax = ctheta*cphi/r
ccc         phix = -sphi/r
         thetax = ctheta*cphi
         phix = -sphi
         ry = stheta*sphi
ccc         thetay = ctheta*sphi/r
ccc         phiy = cphi/r
         thetay = ctheta*sphi
         phiy = cphi
         rz = ctheta
ccc         thetaz = -stheta/r
         thetaz = -stheta
         phiz = 0.0d0
c
c     get the associated Legendre functions:
c
ccc      call ylgndr2s(nterms,ctheta,pp,ppd)
      call ylgndr2sfw(nterms,ctheta,pp,ppd,wlege,nlege)
c
c     get the spherical Bessel functions and their derivatives.
c
 
      ifder=1
      z=wavek*r
      call jfuns3d(jer,nterms,z,rscale,fjs,ifder,fjder,
     1	      lwfjs,iscale,ntop)
      if (jer.ne.0) then
         ier=8
         return
      endif
c
c     scale derivatives of Bessel functions so that they are
c     derivatives with respect to r.
c
      do i=0,nterms
         fjder(i)=fjder(i)*wavek
      enddo
c
c
c     Compute contribution to mpole coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|S|) Ylm*(S) h_n(k|T|)Ylm(T)
c
c     so contribution is j_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of after all source contributions
c     have been included in the calling subroutine h3dformmp.
c
      ur = pp(0,0)*fjder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      ux = ur*rx + utheta*thetax + uphi*phix
      uy = ur*ry + utheta*thetay + uphi*phiy
      uz = ur*rz + utheta*thetaz + uphi*phiz
      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      mpole(0,0)= mpole(0,0) + zzz*dipstr
      do n=1,nterms
         fjsuse = fjs(n+1)*rscale + fjs(n-1)/rscale
         fjsuse = wavek*fjsuse/(2*n+1.0d0)
         ur = pp(n,0)*fjder(n)
         utheta = -fjsuse*ppd(n,0)*stheta
         uphi = 0.0d0
         ux = ur*rx + utheta*thetax + uphi*phix
         uy = ur*ry + utheta*thetay + uphi*phiy
         uz = ur*rz + utheta*thetaz + uphi*phiz
         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         mpole(n,0)= mpole(n,0) + zzz*dipstr
         do m=1,n
            ur = fjder(n)*pp(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fjsuse*ppd(n,m)
            uphi = -eye*m*ephi(-m)*fjsuse*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            mpole(n,m)= mpole(n,m) + zzz*dipstr
c
            ur = fjder(n)*pp(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fjsuse*ppd(n,m)
            uphi = eye*m*ephi(m)*fjsuse*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            mpole(n,-m)= mpole(n,-m) + zzz*dipstr
         enddo
      enddo
c
c
      return
      end
c
c
c
c
c**********************************************************************
      subroutine h3dformta_dp_trunc
     $     (ier,zk,rscale,sources,dipstr,dipvec,ns,
     1		           center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta1/h3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     zk       : Helmholtz coefficient
c     rscale   : scaling parameter
c                     should be less than one in magnitude.
c                     Needed for low frequency regime only
c                     with rsclale abs(wavek) recommended.
c     sources   : coordinates of the sources
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c	 	  ier=4  d is out of range in h3dall
c
c     locexp    : coeffs for the j-expansion
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(3,ns),center(3)
      real *8 dipvec(3,ns)
      complex *16 zk,locexp(0:nterms,-nterms:nterms), dipstr(ns)
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
      do l = 0,nterms
         do m = -l,l
            locexp(l,m) = 0.0d0
         enddo
      enddo
c
      do i = 1,ns
         call h3dformta1_dp_trunc(ier,zk,rscale,sources(1,i),dipstr(i),
     1		dipvec(1,i),center,nterms,nterms1,locexp,wlege,nlege)
      enddo
c
c     scale by (i*k)
c
      do l = 0,nterms
         do m=-l,l
            locexp(l,m) = locexp(l,m)*eye*zk
         enddo
      enddo
C
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dformta_dp_add_trunc
     $     (ier,zk,rscale,sources,dipstr,dipvec,ns,
     1		           center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates a local (j) expansion about the point
c     CENTER due to the NS dipoles at the locations SOURCES(3,*).
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta1/h3dformta0 below.
c
c----------------------------------------------------------------------
c     INPUT:
c
c     zk       : Helmholtz coefficient
c     rscale   : scaling parameter
c                     should be less than one in magnitude.
c                     Needed for low frequency regime only
c                     with rsclale abs(wavek) recommended.
c     sources   : coordinates of the sources
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     ns        : number of sources
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c----------------------------------------------------------------------
c     OUTPUT:
c
c     ier       : error return code
c		  ier=0	returned successfully;
c		  ier=2	insufficient memory in workspace w
c	 	  ier=4  d is out of range in h3dall
c
c     locexp    : coeffs for the j-expansion
c----------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(3,ns),center(3)
      real *8 dipvec(3,ns)
      complex *16 zk,locexp(0:nterms,-nterms:nterms), dipstr(ns)
      complex *16 eye
      complex *16, allocatable :: mptemp(:,:)
      data eye/(0.0d0,1.0d0)/
c
c     initialize local exp
c
        allocate( mptemp(0:nterms,-nterms:nterms) )

        do l = 0,nterms
          do m=-l,l
             mptemp(l,m) = 0
          enddo
        enddo

        call h3dformta_dp_trunc
     $     (ier,zk,rscale,sources,dipstr,dipvec,ns,
     1     center,nterms,nterms1,mptemp,wlege,nlege)
c       
        do l = 0,nterms
          do m=-l,l
            locexp(l,m) = locexp(l,m)+mptemp(l,m)
          enddo
        enddo
C
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine h3dformta1_dp_trunc
     $     (ier,wavek,rscale,source,dipstr,dipvec,
     &		center,nterms,nterms1,locexp,wlege,nlege)
c**********************************************************************
c
c     This subroutine creates the local expansion about CENTER
c     due to a single dipole located at SOURCE.
c     This is the memory management routine. Work is done in the
c     secondary call to h3dformta0 below.
c
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek     : the Helmholtz coefficient
c     rscale    : scaling parameter
c                         should be less than one in magnitude.
c                         Needed for low frequency regime only
c                         with rsclale abs(wavek) recommended.
c     source    : coordinates of the source
c     dipstr    : dipole strengths
c     dipvec    : dipole direction
c     center    : coordinates of the expansion center
c     nterms    : order of the j-expansion
c---------------------------------------------------------------------
c     OUTPUT:
c
c     ier    : error return code
c	           ier=0 successful execution
c		   ier=2 insufficient memory in workspace w
c	 	   ier=4 d is out of range in h3dall
c     locexp : coefficients of the local expansion
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(3),center(3)
      real *8, allocatable :: w(:)
      real *8 dipvec(3)
      complex *16 wavek,locexp(0:nterms,-nterms:nterms), dipstr
c
c     Carve up workspace
c
      ier=0
c
      ipp=1
      lpp=(nterms+1)**2+7
c
      ippd = ipp+lpp
      iephi=ippd+lpp
      lephi=2*(2*nterms+1)+7
c
      ifhs=iephi+lephi
      lfhs=2*(nterms+1)+7
c
      ifhder=ifhs+lfhs
      lfhder=2*(nterms+1)+5
c
      lused=ifhder+lfhder
      allocate(w(lused))
c
      call h3dformta0_dp_trunc(jer,wavek,rscale,source,dipstr,dipvec,
     &   center,nterms,nterms1,locexp,
     $   w(ipp),w(ippd),w(iephi),w(ifhs),w(ifhder),wlege,nlege)
      if (jer.ne.0) ier=4
c
      return
      end
c
c
c
c**********************************************************************
      subroutine h3dformta0_dp_trunc
     $     (ier,wavek,rscale,source,dipstr,dipvec,
     &     center,nterms,nterms1,
     $     locexp,pp,ppd,ephi,fhs,fhder,wlege,nlege)
c**********************************************************************
c
c     See h3dformta_dp/h3dformta1_dp for comments
c
c---------------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      real *8 source(3),center(3),zdiff(3)
      real *8 dipvec(3)
      real *8 pp(0:nterms,0:nterms)
      real *8 ppd(0:nterms,0:nterms)
      complex *16 wavek,locexp(0:nterms,-nterms:nterms), dipstr
      complex *16 ephi(-nterms:nterms),ephi1,ephi1inv
      complex *16 fhs(0:nterms),ztmp,fhder(0:nterms),z
      complex *16 ux,uy,uz,ur,utheta,uphi,zzz
      complex *16 eye
      data eye/(0.0d0,1.0d0)/
c
      ier=0
c
      zdiff(1)=source(1)-center(1)
      zdiff(2)=source(2)-center(2)
      zdiff(3)=source(3)-center(3)
c
      done=1
      call cart2polar(zdiff,r,theta,phi)
      ctheta = dcos(theta)
      stheta=sqrt(done-ctheta*ctheta)
      cphi = dcos(phi)
      sphi = dsin(phi)
      ephi1 = dcmplx(cphi,sphi)
c
c     Compute the e^{eye*m*phi} array
c
      ephi1inv=1.0d0/ephi1
c
      ephi(0)=1.0d0
      ephi(1)=ephi1
      ephi(-1)=ephi1inv
      do i=2,nterms
         ephi(i)=ephi(i-1)*ephi1
         ephi(-i)=ephi(-i+1)*ephi1inv
      enddo
c
c     compute coefficients in change of variables from spherical
c     to Cartesian gradients. In phix, phiy, we leave out the 
c     1/sin(theta) contribution, since we use values of Ynm (which
c     multiplies phix and phiy) that are scaled by 
c     1/sin(theta).
c
        rx = stheta*cphi
        thetax = ctheta*cphi/r
        phix = -sphi/r
        ry = stheta*sphi
        thetay = ctheta*sphi/r
        phiy = cphi/r
        rz = ctheta
        thetaz = -stheta/r
        phiz = 0.0d0
c
c
c     get the associated Legendre functions:
c
ccc      call ylgndr2s(nterms,ctheta,pp,ppd)
      call ylgndr2sfw(nterms,ctheta,pp,ppd,wlege,nlege)
c
c     get the spherical Hankel functions and their derivatives.
c
 
      ifder=1
      z=wavek*r
      call h3dall(nterms,z,rscale,fhs,ifder,fhder)
c
c     scale derivatives of Bessel functions so that they are
c     derivatives with respect to r.
c
      do i=0,nterms
         fhder(i)=fhder(i)*wavek
      enddo
c
c
c     Compute contribution to local coefficients.
c
c     Recall that there are multiple definitions of scaling for
c     Ylm. Using our standard definition, 
c     the addition theorem takes the simple form 
c
c        e^( i k r}/r = 
c         (ik) \sum_n \sum_m  j_n(k|T|) Ylm(T) h_n(k|S|)Ylm*(S)
c
c     so contribution is h_n(k|S|) times
c   
c       Ylm*(S)  = P_l,m * dconjg(ephi(m))               for m > 0   
c       Yl,m*(S)  = P_l,|m| * dconjg(ephi(m))            for m < 0
c                   
c       where P_l,m is the scaled associated Legendre function.
c
c     The factor (i*k) is taken care of after all source contributions
c     have been included in the calling subroutine h3dformmp.
c
      ur = pp(0,0)*fhder(0)
      utheta = 0.0d0
      uphi = 0.0d0
      ux = ur*rx + utheta*thetax + uphi*phix
      uy = ur*ry + utheta*thetay + uphi*phiy
      uz = ur*rz + utheta*thetaz + uphi*phiz
      zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
      locexp(0,0)= locexp(0,0) + zzz*dipstr
      do n=1,nterms
         ur = pp(n,0)*fhder(n)
         utheta = -fhs(n)*ppd(n,0)*stheta
         uphi = 0.0d0
         ux = ur*rx + utheta*thetax + uphi*phix
         uy = ur*ry + utheta*thetay + uphi*phiy
         uz = ur*rz + utheta*thetaz + uphi*phiz
         zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
         locexp(n,0)= locexp(n,0) + zzz*dipstr
         do m=1,n
            ur = fhder(n)*pp(n,m)*stheta*ephi(-m)
            utheta = -ephi(-m)*fhs(n)*ppd(n,m)
            uphi = -eye*m*ephi(-m)*fhs(n)*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            locexp(n,m)= locexp(n,m) + zzz*dipstr
c
            ur = fhder(n)*pp(n,m)*stheta*ephi(m)
            utheta = -ephi(m)*fhs(n)*ppd(n,m)
            uphi = eye*m*ephi(m)*fhs(n)*pp(n,m)
            ux = ur*rx + utheta*thetax + uphi*phix
            uy = ur*ry + utheta*thetay + uphi*phiy
            uz = ur*rz + utheta*thetaz + uphi*phiz
            zzz = dipvec(1)*ux + dipvec(2)*uy + dipvec(3)*uz
            locexp(n,-m)= locexp(n,-m) + zzz*dipstr
         enddo
      enddo
c
      return
      end
c
c
