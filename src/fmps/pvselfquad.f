
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and beginning of the
c       code proper.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This file contains code for generating quadratures for the
c       evaluation of integrals of the form
c
c           \int    K(x,y) \sigma(y) dS(y)                                (1)
c               S
c       where:
c
c           S is a surface element and dS is the corresponding
c           surface measure,
c
c           \sigma is a smooth function.
c           
c           K(x,y) is the double or single layer potential on S, and
c 
c           x is a specified point in S.
c
c       It is assumed that the surface element S is the image under
c       a parameterization
c
c          p: T \to \mathbb{R}^3
c
c       given over a triangle T.  
c
c       The behavior of the Jacobian dp of the parameterization
c       p at the point x has a large influence on the form of the 
c       integrand of (1).  This code proceeds by first composing the 
c       given parameterization p with an appropriate linear mapping 
c
c                 A: \mathbb{R}^2 \to \mathbb{R}^2
c
c       in order to  form a new parameterization p' = p \ocirc A such 
c       that the Jacobian of p' at the point x is conformal.  Then
c       the quadrature rules from radial.f are used to evaluate
c       the resulting integral.
c
c       The quadratures returned by raddiag are formed using precomputed
c       quadrature tables stored on the disk.  These precomputed
c       tables determine the possible orders for the polynomials p
c       and q in (2).  See radial_init for a list of possible orders.
c
c       The following subroutines are user-callable:
c
c   selfquad_init - this routine must be called before 
c       self_quadrature
c
c   selfquad - return a quadrature for evaluating an integral of
c       the form (1)
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine self_quadrature
     $     (ier,verts0,x0,y0,dr,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),a(2,2),ainv(2,2),verts0(2,3)
        dimension xs(1),ys(1),whts(1),b(3,2)
        double precision, allocatable :: rad(:)
        save rad
c
c       Return a quadrature formula for evaluating integrals of the
c       form (1).
c
c                            Input Parameters:
c
c    verts - a (2,3) matrix which specifys the vertices of the triangle
c       over which the surface element is parameterized
c    (x0,y0) - the coordinates (in the parameterization variables)
c       of the target point x
c    dr - the (3,2) Jacobian of the parameterization at (x0,y0)
c
c                           Output Parameters:
c
c    ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 128   means that the quadrature formula could not be
c                   constructed 
c
c    nquad - the number of nodes in the resulting quadrature formula
c    xs - this user-supplied array will contain the x coordinates of 
c       the quadrature nodes upon return
c    ys - this user-supplied array will contain the y coordinates of 
c       the quadrature nodes upone return
c    whts - this user-supplied array will contain the quadrature weights
c       upon return
c
c
        ier = 0
c
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        ainv=0
        call self_findmap(dr,a,ainv)
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        do 2000 i=1,3
        x = verts0(1,i)-x0
        y = verts0(2,i)-y0
c
        xx = ainv(1,1)*x+ainv(1,2)*y
        yy = ainv(2,1)*x+ainv(2,2)*y
c
        verts(1,i) = xx
        verts(2,i) = yy
 2000 continue
c
c       Fetch a quadrature on T_0 for radially singular functions.
c
        call raddiag(jer,rad,verts,nquad,xs,ys,whts)
c
cccccccccccccccccccccccccccccccccccc
c        iplot = 1
c        call plot_quadrature(iplot,verts,nquad,xs,ys)
ccccccccccccccccccccccccccccccccccccc
c
        if (jer .ne. 0) then
        ier = 128
        return
        endif
c
c       Apply the mapping A to the quadrature formula.
c
        det = abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
        sum1=0
        sum2=0
c
        do 3000 i=1,nquad
        s   = xs(i)
        t   = ys(i)
        wht = whts(i)
c
        sum1=sum1+wht
c       
        u = a(1,1)*s + a(1,2)*t
        v = a(2,1)*s + a(2,2)*t
        wht = wht*det
c
        sum2=sum2+wht
c
        xs(i)   = u+x0
        ys(i)   = v+y0
        whts(i) = wht
 3000 continue
c
cccccccccccccccccccccccccccccccccccc
c        iplot = 2
c        call plot_quadrature(iplot,verts0,nquad,xs,ys)
ccccccccccccccccccccccccccccccccccccc

        return

        entry self_quadrature_init(ier0,norder)
c
c       Read the precomputed quadrature formulas from the disc.
c
c                          Input Parameters:
c
c   norder - order for the self-interaction quadratures;
c       right now, the possible orders are 4, 8, 12, 16 
c   lrad - the length of the user-supplied rad array which will store
c       the quadrature tables
c
c                          Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful executionc
c       ier = 4     means that the rad array was of insufficient length
c       ier = 128   means that the quadrature table for the specified
c                   order could not be found
c       ier = 1024  means that their was an I/O error reading the
c                   table from the disk; this usually means that the
c                   data file containing it was improperly formatted
c
        ier0 = 0
        lrad = 2 000 000
        if(allocated(rad)) deallocate(rad)
        allocate(rad(lrad))
c
        call radial_init(jer,norder,rad,lrad,lkeep)
        if (jer .ne. 0) then
        ier0 = 1024
        return
        endif
c
        return
        end




        subroutine selfquad0(ier,n,nlege,verts0,x0,y0,
     -    nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts0(2,3),verts(2,3),xs(1),ys(1),whts(1)
c
        do 1000 i=1,3
        verts(1,i) = verts0(1,i)-x0
        verts(2,i) = verts0(2,i)-y0
 1000 continue
c
        call radquad(ier,verts,n,nlege,nquad,xs,ys,whts)
c
        do 2000 i=1,nquad
        xs(i) = xs(i)+x0
        ys(i) = ys(i)+y0
 2000 continue
c
        end


        subroutine selfquad(ier,rad,verts0,x0,y0,dr,nquad,xs,ys,whts)
        implicit double precision (a-h,o-z)
        dimension verts(2,3),dr(3,2),a(2,2),ainv(2,2),verts0(2,3)
        dimension xs(1),ys(1),whts(1),b(3,2)
c
c       Return a quadrature formula for evaluating integrals of the
c       form (1).
c
c                            Input Parameters:
c
c    verts - a (2,3) matrix which specifys the vertices of the triangle
c       over which the surface element is parameterized
c    (x0,y0) - the coordinates (in the parameterization variables)
c       of the target point x
c    dr - the (3,2) Jacobian of the parameterization at (x0,y0)
c
c                           Output Parameters:
c
c    ier - an error return code;
c       ier = 0     indicates successful execution
c       ier = 128   means that the quadrature formula could not be
c                   constructed 
c
c    nquad - the number of nodes in the resulting quadrature formula
c    xs - this user-supplied array will contain the x coordinates of 
c       the quadrature nodes upon return
c    ys - this user-supplied array will contain the y coordinates of 
c       the quadrature nodes upone return
c    whts - this user-supplied array will contain the quadrature weights
c       upon return
c
c
        ier = 0
c
c
c       Find a linear map A which makes the Jacobian at the target node
c       conformal.  Also, compute the inverse of A.
c
        ainv=0
        call self_findmap(dr,a,ainv)
c
c       Apply the inverse of A to the vertices of the triangle T in
c       order to find the vertices of the triangle T_0.
c
        do 2000 i=1,3
        x = verts0(1,i)-x0
        y = verts0(2,i)-y0
c
        xx = ainv(1,1)*x+ainv(1,2)*y
        yy = ainv(2,1)*x+ainv(2,2)*y
c
        verts(1,i) = xx
        verts(2,i) = yy
 2000 continue
c
c       Fetch a quadrature on T_0 for radially singular functions.
c
        call raddiag(jer,rad,verts,nquad,xs,ys,whts)
c
cccccccccccccccccccccccccccccccccccc
c        iplot = 1
c        call plot_quadrature(iplot,verts,nquad,xs,ys)
ccccccccccccccccccccccccccccccccccccc
c
        if (jer .ne. 0) then
        ier = 128
        return
        endif
c
c       Apply the mapping A to the quadrature formula.
c
        det = abs(a(1,1)*a(2,2)-a(1,2)*a(2,1))
c
        sum1=0
        sum2=0
c
        do 3000 i=1,nquad
        s   = xs(i)
        t   = ys(i)
        wht = whts(i)
c
        sum1=sum1+wht
c       
        u = a(1,1)*s + a(1,2)*t
        v = a(2,1)*s + a(2,2)*t
        wht = wht*det
c
        sum2=sum2+wht
c
        xs(i)   = u+x0
        ys(i)   = v+y0
        whts(i) = wht
 3000 continue
c
cccccccccccccccccccccccccccccccccccc
c        iplot = 2
c        call plot_quadrature(iplot,verts0,nquad,xs,ys)
ccccccccccccccccccccccccccccccccccccc

        return
        end



        subroutine selfquad_init(ier,norder,rad,lrad,lkeep)
        implicit double precision (a-h,o-z)
        dimension rad(1)
c
c       Read the precomputed quadrature formulas from the disc.
c
c                          Input Parameters:
c
c   norder - order for the self-interaction quadratures;
c       right now, the possible orders are 4, 8, 12, 16 
c   lrad - the length of the user-supplied rad array which will store
c       the quadrature tables
c
c                          Output Parameters:
c
c   ier - an error return code;
c       ier = 0     indicates successful executionc
c       ier = 4     means that the rad array was of insufficient length
c       ier = 128   means that the quadrature table for the specified
c                   order could not be found
c       ier = 1024  means that their was an I/O error reading the
c                   table from the disk; this usually means that the
c                   data file containing it was improperly formatted
c
        ier = 0
        call radial_init(ier,norder,rad,lrad,lkeep)
        return
        end


        subroutine self_findmap(a,b,binv)
        implicit double precision (a-h,o-z)
        dimension a(3,2),b(2,2),binv(2,2),q(3,2)
c
c       Given a (3,2) matrix A, return a (2,2) matrix B such that
c       A*B has orthonormal columns.  Also, return the inverse
c       BINV of the mapping B.
c
c
c                             --- WARNING ---
c       THIS ROUTINE IS A BIT ROUGH --- IT IS EXPECTED TO FAIL
c       IN THE EVENT THAT THE MATRIX A IS NEARLY SINGULAR.  BUT
c       WHAT ARE YOU DOING USING PARAMETERIZATION WITH JACOBIANS
c       WHICH ARE NEARLY SINGULAR ANYWAY?
c                              --------------
c
c
c                          Input Parameters:
c
c   a - the (3,2) input matrix 
c
c                         Output Parameters:
c
c   b - a (2,2) matrix such that a*b has orthonormal columns
c   binv - the (2,2) inverse of b
c
c
c       Orthonormalize the columns of the matrix a.
c
        sum1 = sqrt(a(1,1)**2 + a(2,1)** 2 + a(3,1)**2)
        q(1,1) = a(1,1) / sum1
        q(2,1) = a(2,1) / sum1
        q(3,1) = a(3,1) / sum1
c
        dip = a(1,2)*q(1,1)+a(2,2)*q(2,1)+a(3,2)*q(3,1)
c
        q(1,2) = a(1,2) - dip * q(1,1)
        q(2,2) = a(2,2) - dip * q(2,1)
        q(3,2) = a(3,2) - dip * q(3,1)
c
        sum2 = sqrt(q(1,2)**2 + q(2,2)** 2 + q(3,2)**2)
        q(1,2) = q(1,2) / sum2
        q(2,2) = q(2,2) / sum2
        q(3,2) = q(3,2) / sum2
c
c       Compute BINV = Q'*A.
c     
        binv(1,1) = q(1,1)*a(1,1) + q(2,1)*a(2,1) + q(3,1)*a(3,1)
        binv(1,2) = q(1,1)*a(1,2) + q(2,1)*a(2,2) + q(3,1)*a(3,2)
        binv(2,1) = q(1,2)*a(1,1) + q(2,2)*a(2,1) + q(3,2)*a(3,1)
        binv(2,2) = q(1,2)*a(1,2) + q(2,2)*a(2,2) + q(3,2)*a(3,2)
c
c       Compute B = (BINV)^(-1).
c
        det = binv(1,1)*binv(2,2) - binv(1,2)*binv(2,1)
        b(1,1) = binv(2,2)/det
        b(2,2) = binv(1,1)/det
        b(1,2) = -binv(1,2)/det
        b(2,1) = -binv(2,1)/det
c
        end

