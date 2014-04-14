        subroutine em3sphlin(ampole,nterms,cvec)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 cvec(1)
c
c       ... compress multipole expansion storage into linear array
c
        kk=0
        do n=0,nterms
        do m=-n,n
        kk=kk+1
        cvec(kk)=ampole(n,m)
        enddo
        enddo
c
ccc        call prinf('kk=*',kk,1)
c
        return
        end
c
c
c
c
c
        subroutine em3linsph(ampole,nterms,cvec)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 cvec(1)
c
c       ... unroll linear storage array into multipole expansion 
c
        kk=0
        do n=0,nterms
        do m=-n,n
        kk=kk+1
        ampole(n,m)=cvec(kk)
        enddo
        enddo
c
        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine empol(wavelength,rk,
     $     ampole,bmpole,nterms,radius,pvec,mvec)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 ima,rk,cd
        complex *16 pvec(3),mvec(3)
        data ima/(0.0d0,1.0d0)/

        done=1
        pi=4*atan(done)

c
c       ... far field signature factor
c
        cd = -ima / rk
c
c       ... size normalization factor, optional
c
ccc        cd = cd / radius**3 
ccc        cd = cd / (4.0d0/3.0d0) 
c
c
c       ... evaluate p vector
c
        pvec(1)=
     $     (bmpole(1,1)+bmpole(1,-1))/2*sqrt(3.0d0)/rk**2*(-ima)*cd
        pvec(2)=
     $     (bmpole(1,1)-bmpole(1,-1))/2*sqrt(3.0d0)/rk**2*cd
        pvec(3)=
     $     (bmpole(1,0))/sqrt(2.0d0)*sqrt(3.0d0)/rk**2*ima*cd
c
ccc        write(54,1000) wavelength, pvec(1),pvec(2),pvec(3)
c
c
c       ... evaluate m vector
c
        mvec(1)=
     $     (ampole(1,1)+ampole(1,-1))/2*sqrt(3.0d0)/rk**2*(-ima)*cd
        mvec(2)=
     $     (ampole(1,1)-ampole(1,-1))/2*sqrt(3.0d0)/rk**2*cd
        mvec(3)=
     $     (ampole(1,0))/sqrt(2.0d0)*sqrt(3.0d0)/rk**2*ima*cd
c
ccc        write(55,1000) wavelength, mvec(1),mvec(2),mvec(3)

 1000   format(7(1x,e11.5))

        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
        subroutine em3dmpta(nspheres,nterms,ncoefs,omega,eps0,cmu0,
     $     sphere_xyz,sphere_r,aompole,bompole,asmpole,bsmpole,
     $     rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
        dimension sphere_xyz(3,nspheres),sphere_r(nspheres)
        complex *16 aompole(ncoefs,nspheres)
        complex *16 bompole(ncoefs,nspheres)
        complex *16 asmpole(ncoefs,nspheres)
        complex *16 bsmpole(ncoefs,nspheres)
        complex *16 eps0,cmu0,rk
c
c
c       ... prepare to evaluate the incoming scattered field
c
        do 1190 kk=1,nspheres
        do 1180 ii=1,ncoefs
        asmpole(ii,kk)=0
        bsmpole(ii,kk)=0
 1180   continue
 1190   continue
c
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
c
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(nn)
        do 2200 kk=1,nspheres
c
        do 2100 nn=1,nspheres
        if( nn .eq. kk ) goto 2100
c
        call em3mpta3_add(rk,
     $     sphere_xyz(1,nn),aompole(1,nn),bompole(1,nn),
     $     nterms,
     $     sphere_xyz(1,kk),asmpole(1,kk),bsmpole(1,kk),
     $     nterms,
     $     sphere_r(kk),rnodes,weights,nphi,ntheta)
c       
 2100   continue
c       
 2200   continue
C$OMP END PARALLEL DO
c
c
        return
        end
c
c
c
c
c
        subroutine em3dtatab(nspheres,omega,eps0,cmu0,
     $     nterms,ncoefs,sphere_xyz,sphere_r,aompole,bompole,
     $     nterms0,ncoefs0,sphere_xyz0,sphere_r0,aompole0,bompole0,
     $     rnodes,weights,nphi,ntheta)
        implicit real *8 (a-h,o-z)
        dimension sphere_xyz(3,nspheres),sphere_r(nspheres)
        complex *16 aompole(ncoefs,nspheres)
        complex *16 bompole(ncoefs,nspheres)
        complex *16 aompole0(ncoefs0)
        complex *16 bompole0(ncoefs0)
        complex *16 eps0,cmu0,rk
        complex *16, allocatable :: atemp(:)
        complex *16, allocatable :: btemp(:)
c
c
c       ... prepare to evaluate the incoming scattered field
c
        do i=1,ncoefs0
        aompole0(ii)=0
        bompole0(ii)=0
        enddo
c
        allocate( atemp(ncoefs0) )
        allocate( btemp(ncoefs0) )
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
c
c
cccC$OMP PARALLEL DO DEFAULT(SHARED)
        do 2100 nn=1,nspheres
c
        call em3tata3trunc(rk,
     $     sphere_xyz(1,nn),aompole(1,nn),bompole(1,nn),
     $     nterms,
     $     sphere_xyz0,atemp,btemp,
     $     nterms0,
     $     sphere_r0,rnodes,weights,nphi,ntheta)
c       
        do i=1,ncoefs0
        aompole0(i)=aompole0(i)+atemp(i)
        bompole0(i)=bompole0(i)+btemp(i)
        enddo
c
 2100   continue
cccC$OMP END PARALLEL DO
c
c
        return
        end
c
c
c
c
c
        subroutine em3dmpoletargeval
     $     (nspheres,nterms,ncoefs,omega,eps0,cmu0,
     $     sphere_xyz,sphere_r,aompole,bompole,
     $     ntargets,targets,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension sphere_xyz(3,nspheres),sphere_r(3,nspheres)
        complex *16 aompole(ncoefs,nspheres)
        complex *16 bompole(ncoefs,nspheres)
        complex *16 eps0,cmu0,rk,evec0(3),hvec0(3)
        dimension targets(3,ntargets)
        complex *16 evec(3,ntargets)
        complex *16 hvec(3,ntargets)
c      
c       ... outgoing scattered fields 
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(kk,evec0,hvec0)
        do 2200 nn=1,ntargets
c
        evec(1,nn)=0
        evec(2,nn)=0
        evec(3,nn)=0
c
        hvec(1,nn)=0
        hvec(2,nn)=0
        hvec(3,nn)=0
c
        do 2100 kk=1,nspheres
c
c       ... scattered field due to each sphere
c
        call em3mpeval
     $     (rk,sphere_xyz(1,kk),aompole(1,kk),bompole(1,kk),nterms,
     $     targets(1,nn),evec0,hvec0)
c
        evec(1,nn)=evec(1,nn)+evec0(1)
        evec(2,nn)=evec(2,nn)+evec0(2)
        evec(3,nn)=evec(3,nn)+evec0(3)
c
        hvec(1,nn)=hvec(1,nn)+hvec0(1)
        hvec(2,nn)=hvec(2,nn)+hvec0(2)
        hvec(3,nn)=hvec(3,nn)+hvec0(3)
c
 2100   continue
c
 2200   continue
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine em3dlocaltargeval
     $     (nspheres,nterms,ncoefs,omega,eps0,cmu0,
     $     sphere_xyz,sphere_r,aimpole,bimpole,
     $     ntargets,targets,evec,hvec)
        implicit real *8 (a-h,o-z)
        dimension sphere_xyz(3,nspheres),sphere_r(3,nspheres)
        complex *16 aimpole(ncoefs,nspheres)
        complex *16 bimpole(ncoefs,nspheres)
        complex *16 eps0,cmu0,rk,evec0(3),hvec0(3)
        dimension targets(3,ntargets)
        complex *16 evec(3,ntargets)
        complex *16 hvec(3,ntargets)
c      
c       ... incoming fields 
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
c
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(kk,evec0,hvec0)
        do 2200 nn=1,ntargets
c
        evec(1,nn)=0
        evec(2,nn)=0
        evec(3,nn)=0
c
        hvec(1,nn)=0
        hvec(2,nn)=0
        hvec(3,nn)=0
c
        do 2100 kk=1,nspheres
c
c       ... incoming field due to each sphere
c
        call em3taeval
     $     (rk,sphere_xyz(1,kk),aimpole(1,kk),bimpole(1,kk),nterms,
     $     targets(1,nn),evec0,hvec0)
c
        evec(1,nn)=evec(1,nn)+evec0(1)
        evec(2,nn)=evec(2,nn)+evec0(2)
        evec(3,nn)=evec(3,nn)+evec0(3)
c
        hvec(1,nn)=hvec(1,nn)+hvec0(1)
        hvec(2,nn)=hvec(2,nn)+hvec0(2)
        hvec(3,nn)=hvec(3,nn)+hvec0(3)
c
 2100   continue
c
 2200   continue
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine emplanewtargeval(rk,cjvec,cmvec,
     $     ntargets,targets,evec,hvec)
        implicit real *8 (a-h,o-z)
        complex *16 cjvec(3),cmvec(3),rk,evec0(3),hvec0(3)
        dimension targets(3,ntargets)
        complex *16 evec(3,ntargets)
        complex *16 hvec(3,ntargets)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(kk,evec0,hvec0)
        do 2200 nn=1,ntargets
c
        evec(1,nn)=0
        evec(2,nn)=0
        evec(3,nn)=0
c
        hvec(1,nn)=0
        hvec(2,nn)=0
        hvec(3,nn)=0
c
c       ... scattered field due to each plane wave
c
        call emplanew(rk,targets(1,nn),cjvec,cmvec,
     $     evec0,hvec0)
c
        evec(1,nn)=evec(1,nn)+evec0(1)
        evec(2,nn)=evec(2,nn)+evec0(2)
        evec(3,nn)=evec(3,nn)+evec0(3)
c
        hvec(1,nn)=hvec(1,nn)+hvec0(1)
        hvec(2,nn)=hvec(2,nn)+hvec0(2)
        hvec(3,nn)=hvec(3,nn)+hvec0(3)
c
 2200   continue
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine emplanearbtargeval(rkvec,epol,
     $     ntargets,targets,evec,hvec)
        implicit real *8 (a-h,o-z)
        complex *16 rkvec(3),epol(3),rk,evec0(3),hvec0(3)
        dimension targets(3,ntargets)
        complex *16 evec(3,ntargets)
        complex *16 hvec(3,ntargets)

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(kk,evec0,hvec0)
        do 2200 nn=1,ntargets
c
        evec(1,nn)=0
        evec(2,nn)=0
        evec(3,nn)=0
c
        hvec(1,nn)=0
        hvec(2,nn)=0
        hvec(3,nn)=0
c
c       ... scattered field due to each plane wave
c
        call emplanearb(rkvec,epol,targets(1,nn),
     $     evec0,hvec0)
c
        evec(1,nn)=evec(1,nn)+evec0(1)
        evec(2,nn)=evec(2,nn)+evec0(2)
        evec(3,nn)=evec(3,nn)+evec0(3)
c
        hvec(1,nn)=hvec(1,nn)+hvec0(1)
        hvec(2,nn)=hvec(2,nn)+hvec0(2)
        hvec(3,nn)=hvec(3,nn)+hvec0(3)
c
 2200   continue
C$OMP END PARALLEL DO
c
        return
        end
c
c
c
c
c
        subroutine em3orient
     $     (sphere_rot,nspheres,frames)
        implicit real *8 (a-h,o-z)
        real *8 sphere_rot(3,1),frames(3,3,1)
        dimension vert1(3),vert2(3),vert3(3)
c
c       ... orientation of spheres 
c
        do kk=1,nspheres
        call framestd(vert1,vert2,vert3)
        call framerotf(vert1,vert2,vert3,
     $     sphere_rot(1,kk),sphere_rot(2,kk),sphere_rot(3,kk))
        frames(1,1,kk)=vert1(1)
        frames(2,1,kk)=vert1(2)
        frames(3,1,kk)=vert1(3)
        frames(1,2,kk)=vert2(1)
        frames(2,2,kk)=vert2(2)
        frames(3,2,kk)=vert2(3)
        frames(1,3,kk)=vert3(1)
        frames(2,3,kk)=vert3(2)
        frames(3,3,kk)=vert3(3)
        enddo
c
        return
        end
c
c
c
c
c
        subroutine empol_one_sphere(wavelength,rk,
     $     ampole,bmpole,nterms,radius,pvec,mvec)
        implicit real *8 (a-h,o-z)
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 ima,rk,cd
        complex *16 pvec(3),mvec(3)
        data ima/(0.0d0,1.0d0)/

        done=1
        pi=4*atan(done)

c
c       ... far field signature factor
c
        cd = -ima / rk
c
c       ... size normalization factor, optional
c
ccc        cd = cd / radius**3 
ccc        cd = cd / (4.0d0/3.0d0) 
c
c
c       ... evaluate p vector
c
        pvec(1)=
     $     (bmpole(1,1)+bmpole(1,-1))/2*sqrt(3.0d0)/rk**2*(-ima)*cd
        pvec(2)=
     $     (bmpole(1,1)-bmpole(1,-1))/2*sqrt(3.0d0)/rk**2*cd
        pvec(3)=
     $     (bmpole(1,0))/sqrt(2.0d0)*sqrt(3.0d0)/rk**2*ima*cd
c
        write(54,1000) wavelength, pvec(1),pvec(2),pvec(3)
c
c
c       ... evaluate m vector
c
        mvec(1)=
     $     (ampole(1,1)+ampole(1,-1))/2*sqrt(3.0d0)/rk**2*(-ima)*cd
        mvec(2)=
     $     (ampole(1,1)-ampole(1,-1))/2*sqrt(3.0d0)/rk**2*cd
        mvec(3)=
     $     (ampole(1,0))/sqrt(2.0d0)*sqrt(3.0d0)/rk**2*ima*cd
c
        write(55,1000) wavelength, mvec(1),mvec(2),mvec(3)

 1000   format(7(1x,e11.5))

        return
        end
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
        subroutine qsca_one_sphere
     $     (omega,eps0,cmu0,ampole,bmpole,nterms,radius)
        implicit real *8 (a-h,o-z)
        complex *16 eps0,cmu0,rk,cd
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
c
        complex *16 z0
c
        complex *16, allocatable :: jvals0(:),hvals0(:)
        complex *16, allocatable :: jders0(:),hders0(:)
        complex *16, allocatable :: emtjvals0(:),emrjvals0(:)
        complex *16, allocatable :: emthvals0(:),emrhvals0(:)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( jvals0(nterms+1) )
        allocate( jders0(nterms+1) )
        allocate( hvals0(nterms+1) )
        allocate( hders0(nterms+1) )
        allocate( emtjvals0(nterms+1) ) 
        allocate( emrjvals0(nterms+1) ) 
        allocate( emthvals0(nterms+1) ) 
        allocate( emrhvals0(nterms+1) ) 
c
c       ... test 
c
        done=1
        pi=4*atan(done)
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
        z0=rk*radius
cc        call prin2('rk=*',rk,2)
cc        call prin2('z0=*',z0,2)
cc        call prin2('radius=*',radius,1)

        call emjevalrt(nterms,z0,jvals0,jders0,emtjvals0,emrjvals0)
        call emhevalrt(nterms,z0,hvals0,hders0,emthvals0,emrhvals0)

        cd=0
        do n=1,nterms
        do m=-n,n

c        cd = cd + (2*n+1)*abs(ampole(n,m)*emthvals0(n+1))**2
c        cd = cd + (2*n+1)*abs(bmpole(n,m)*emthvals0(n+1))**2

c        cd = cd + dble(2*n+1)*abs(ampole(n,m))**2 
c        cd = cd + dble(2*n+1)*abs(bmpole(n,m))**2 

        cd = cd + abs(ampole(n,m))**2 
        cd = cd + abs(bmpole(n,m))**2 

        enddo
        enddo
c
ccc        cd =cd * 0.40528E+00
c        cd = cd * (2/pi)**2
c        cd = cd/rk**2

        wavelength = 2*pi/omega
        alpha = 2*pi*radius/wavelength
        cd = cd * 4/alpha**2  
        call prin2('directly, qsca, cd=*',cd,1)
c
c octave:17> 1/sqrt(0.40528E+00)
c ans =  1.5708
c octave:18> (pi/2)
c ans =  1.5708
c
c octave:21> (2/pi)**2
c ans =  0.40528

        return
        end
c
c
c
c
c
        subroutine qext_one_sphere
     $     (omega,eps0,cmu0,ampole,bmpole,aimpole,bimpole,nterms,radius)
        implicit real *8 (a-h,o-z)
        complex *16 eps0,cmu0,rk,cd
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 aimpole(0:nterms,-nterms:nterms)
        complex *16 bimpole(0:nterms,-nterms:nterms)
c
        complex *16 z0
        complex *16, allocatable :: jvals0(:),hvals0(:)
        complex *16, allocatable :: jders0(:),hders0(:)
        complex *16, allocatable :: emtjvals0(:),emrjvals0(:)
        complex *16, allocatable :: emthvals0(:),emrhvals0(:)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
        allocate( jvals0(nterms+1) )
        allocate( jders0(nterms+1) )
        allocate( hvals0(nterms+1) )
        allocate( hders0(nterms+1) )
        allocate( emtjvals0(nterms+1) ) 
        allocate( emrjvals0(nterms+1) ) 
        allocate( emthvals0(nterms+1) ) 
        allocate( emrhvals0(nterms+1) ) 
c
c       ... test 
c
        done=1
        pi=4*atan(done)
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
        z0=rk*radius
cc        call prin2('rk=*',rk,2)
cc        call prin2('z0=*',z0,2)

        call emjevalrt(nterms,z0,jvals0,jders0,emtjvals0,emrjvals0)
        call emhevalrt(nterms,z0,hvals0,hders0,emthvals0,emrhvals0)

        cd=0
        do n=1,nterms
        do m=-n,n

c        cd = cd + (2*n+1)*abs(ampole(n,m)*emthvals0(n+1))**2
c        cd = cd + (2*n+1)*abs(bmpole(n,m)*emthvals0(n+1))**2

c        cd = cd + dble(2*n+1)*abs(ampole(n,m))**2 
c        cd = cd + dble(2*n+1)*abs(bmpole(n,m))**2 

        cd = cd + dble(ampole(n,m)*conjg(aimpole(n,m)))
        cd = cd + dble(bmpole(n,m)*conjg(bimpole(n,m)))

        enddo
        enddo
c
        cd = -cd
c
ccc        cd =cd * 0.40528E+00
c        cd = cd * (2/pi)**2
c        cd = cd/rk**2

        wavelength = 2*pi/omega
        alpha = 2*pi*radius/wavelength
        cd = cd * 4/alpha**2  
        call prin2('directly, qext, cd=*',cd,1)
c
c octave:17> 1/sqrt(0.40528E+00)
c ans =  1.5708
c octave:18> (pi/2)
c ans =  1.5708
c
c octave:21> (2/pi)**2
c ans =  0.40528

        return
        end
c
c
c
c
c
        subroutine qscat_diff_one_sphere
     $     (omega,eps0,cmu0,ampole,bmpole,aimpole,bimpole,nterms,radius)
        implicit real *8 (a-h,o-z)
        complex *16 eps0,cmu0,rk,cd
        complex *16 ampole(0:nterms,-nterms:nterms)
        complex *16 bmpole(0:nterms,-nterms:nterms)
        complex *16 aimpole(0:nterms,-nterms:nterms)
        complex *16 bimpole(0:nterms,-nterms:nterms)
c
        complex *16 ima
        data ima/(0.0d0,1.0d0)/
c
c       ... test 
c
        done=1
        pi=4*atan(done)
c
        rk=omega*sqrt(eps0)*sqrt(cmu0)
        z0=rk*radius

        cd=0

        do n=1,nterms
        do m=-n,n

        cd = cd + dble(
     $     (aimpole(n,m))*conjg(ampole(n,m)))
        cd = cd + dble(
     $     (bimpole(n,m))*conjg(bmpole(n,m)))

        enddo
        enddo
c
c
ccc        cd =cd * 0.40528E+00
c        cd = cd * (2/pi)**2
c        cd = cd/rk**2

        wavelength = 2*pi/omega
        alpha = 2*pi*radius/wavelength
        cd = cd * 4/alpha**2  
        call prin2('directly, qscat_diff, cd=*',cd,1)
c
c octave:17> 1/sqrt(0.40528E+00)
c ans =  1.5708
c octave:18> (pi/2)
c ans =  1.5708
c
c octave:21> (2/pi)**2
c ans =  0.40528

        return
        end
c
c
c
c
c
