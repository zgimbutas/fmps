c
c
        subroutine patchmatcanaf90(itype_geotriainfo,npatches,
     $     norder,noversamp,cmatr,
     $     xyzs,xyznorms,xyztang1s,xyztang2s,whts,
     $     ikernel,rk)
        implicit real *8 (a-h,o-z)
c
        external fpatchpnt,qpatchpnt,cpatchpnt
        external rpatchpnt,spatchpnt
c
        external lfinter1,lfinter2,lfinter3,lfinter4
        external hfinter1,hfinter2,hfinter3,hfinter4
        external hfinter1m,hfinter2m,hfinter3m,hfinter4m
c
        dimension usout(100),vsout(100),wsout(100)
        dimension umatr(10000),vmatr(10000)
c
        dimension verts(3,100000),ifaces(3,100000)

        dimension triainfo(3,3,100000)
        complex *16 cmatr(*)
c
        allocatable :: ipatchinfo(:),refineinfo(:,:)
c
        allocatable :: ixyzs(:,:)
        real *8 xyzs(3,*),xyznorms(3,*)
        real *8 xyztang1s(3,*),xyztang2s(3,*)
        real *8 whts(*)
c
        complex *16 rk
c
        real *8, allocatable :: w(:)
c
c
c
        lw = 50 000 000
        allocate(w(lw))

c
c       ... retrieve the interpolation nodes
c
        itype=1
ccc        norder=5
        call ortho2siexps(itype,norder,npols,usout,vsout,
     1     umatr,vmatr,wsout)
c
        call prinf('norder=*',norder,1)
        call prinf('npols=*',npols,1)
        call prin2('usout=*',usout,npols)
        call prin2('vsout=*',vsout,npols)
        call prin2('wsout=*',wsout,npols)
c
        d=0
        do 1100 i=1,npols
        d=d+wsout(i)
 1100   continue
c
        call prin2('sum of weights=*',d,1)
c

        call prinf('npatches=*',npatches,1)
c
c       ... retrieve a triangulation from a library
c
        itype=3
        call rsolid(itype,verts,nverts,ifaces,nfaces)
        call gentriainfo(verts,nverts,ifaces,nfaces,triainfo)
        npatches=nfaces
c
        call prinf('npatches=*',npatches,1)
        call prinf('nverts=*',nverts,1)
        call prin2('verts=*',verts,3*nverts)
        call prinf('nfaces=*',nfaces,1)
        call prinf('ifaces=*',ifaces,3*nfaces)
        call prin2('triainfo=*',triainfo,3*3*npatches)


c       ... refine the triangulation
c
ccc        noversamp=2
c
        max_ref_tri = npatches*noversamp**2
        allocate( ipatchinfo(max_ref_tri), refineinfo(4,max_ref_tri) )
c
        call genrefineinfo(noversamp,npatches,
     $     npatchesout,ipatchinfo,refineinfo)
c
ccc        npatches=npatchesout
c

c
c       ... map the interpolation points on the patch into R^3 
c        
        npts=npols*npatchesout
        call prinf('npts=*',npts,1)
        call prinf('npatchesout=*',npatchesout,1)
        call prinf('npols=*',npols,1)

c
c       ... allocate work arrays for discretization 
c
        allocate( ixyzs(2,npts) )
c
c
c       ... call patchmatc
c
        call patchallpnts(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     npols,usout,vsout,ixyzs,xyzs,xyznorms,
     $     xyztang1s,xyztang2s,npts)
        call patchallwhts(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     npols,usout,vsout,wsout,whts,npts)
c
        do i=1,npts
        d=d+whts(i)
        enddo
c
        call prin2('after patchallwhts, sum of weights=*',d,1)
c
        if( ikernel .eq. 1 ) then
        call prinf('patchmatc: S_0*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     lfinter1,par5,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 2 ) then
        call prinf('patchmatc: D_0*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     lfinter2,par5,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 3 ) then
        call prinf('patchmatc: S_0 prime*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     lfinter3,par5,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 4 ) then
        call prinf('patchmatc: D_0 prime*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     lfinter4,par5,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
c
        if( ikernel .eq. 11 ) then
        call prinf('patchmatc: S_k*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter1,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 12 ) then
        call prinf('patchmatc: D_k*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter2,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 13 ) then
        call prinf('patchmatc: S_k prime*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter3,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 14 ) then
        call prinf('patchmatc: D_k prime*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter4,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
c
        if( ikernel .eq. 15 ) then
        call prinf('patchmatc: S_k-S_0*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter1m,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 16 ) then
        call prinf('patchmatc: D_k-D_0*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter2m,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 17 ) then
        call prinf('patchmatc: S_k prime - S_0 prime*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter3m,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
        if( ikernel .eq. 18 ) then
        call prinf('patchmatc: D_k prime - D_0 prime*',id,0)
        call patchmatc(npatchesout,rpatchpnt,
     $     triainfo,spatchpnt,ipatchinfo,refineinfo,
     $     norder,npols,usout,vsout,umatr,vmatr,
     $     ixyzs,xyzs,xyznorms,xyztang1s,xyztang2s,npts,
     $     hfinter4m,rk,par6,par7,par8,
     $     cmatr,w,lw,lused,ier)
        endif
c
c
        return
        end
c
c
c
c
c
