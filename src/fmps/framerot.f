c
c
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c       This is the end of the debugging code and the beginning 
c       of the frame rotation routines 
c       
c       Euler angle routines
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
c
c
c
        subroutine framestd(vert1,vert2,vert3)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine constructs (or rather retrieves) the
c       on the standard simplex in R^3 with the vertices
c
c       (1,0,0), (0,1,0), (0,0,1)
c
c       Output parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of simplex
c
c
        vert1(1)=1
        vert1(2)=0
        vert1(3)=0
c        
        vert2(1)=0
        vert2(2)=1
        vert2(3)=0
c        
        vert3(1)=0
        vert3(2)=0
        vert3(3)=1
c        
        return
        end
c
c
c
c
c
        subroutine framerotf(vert1,vert2,vert3,alpha,beta,gamma)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine performs the forward euler rotation of the
c       simplex frame
c
c       (1,0,0), (0,1,0), (0,0,1)
c
c       The rotations are performed in the following order (all in the
c       original coordinate system!)
c
c       First, perform a rotation around the z axis by the angle gamma
c       Then, perform a rotation around the y axis by the angle beta 
c       Finally, perform a rotation around the z axis by the angle alpha
c       
c
c
c       Input parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of simplex
c                     (these are changed on output)
c
c       Euler angles:
c       alpha   - the angle of the first rotation about z axis
c       beta    - the angle of the rotation about y axis
c       gamma   - the angle of the second rotation about z axis
c
c       Output parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of
c                                     the rotated frame
c
c
        call framerotz(vert1,vert2,vert3,gamma)
c        call prin2('vert1=*',vert1,3)
c        call prin2('vert2=*',vert2,3)
c        call prin2('vert3=*',vert3,3)
c        write(*,*)
        call frameroty(vert1,vert2,vert3,beta)
c        call prin2('vert1=*',vert1,3)
c        call prin2('vert2=*',vert2,3)
c        call prin2('vert3=*',vert3,3)
c        write(*,*)
        call framerotz(vert1,vert2,vert3,alpha)
c        call prin2('vert1=*',vert1,3)
c        call prin2('vert2=*',vert2,3)
c        call prin2('vert3=*',vert3,3)
c        write(*,*)
c
        return
        end
c
c
c
c
c
        subroutine framerotb(vert1,vert2,vert3,alpha,beta,gamma)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine performs the backward euler rotation of the
c       simplex frame
c
c       (1,0,0), (0,1,0), (0,0,1)
c
c       The rotations are performed in the following order (all in the
c       original coordinate system!)
c
c       First, perform a rotation around the z axis by the angle -alpha 
c       Then, perform a rotation around the y axis by the angle -beta 
c       Finally, perform a rotation around the z axis by the angle -gamma 
c
c       Input parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of simplex
c                     (these are changed on output)
c
c       Euler angles:
c       alpha   - the angle of the first rotation about z axis
c       beta    - the angle of the rotation about y axis
c       gamma   - the angle of the second rotation about z axis
c
c       Output parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of the 
c                                     rotated frame
c
c
        call framerotz(vert1,vert2,vert3,-alpha)
c        call prin2('vert1=*',vert1,3)
c        call prin2('vert2=*',vert2,3)
c        call prin2('vert3=*',vert3,3)
c        write(*,*)
        call frameroty(vert1,vert2,vert3,-beta)
c        call prin2('vert1=*',vert1,3)
c        call prin2('vert2=*',vert2,3)
c        call prin2('vert3=*',vert3,3)
c        write(*,*)
        call framerotz(vert1,vert2,vert3,-gamma)
c        call prin2('vert1=*',vert1,3)
c        call prin2('vert2=*',vert2,3)
c        call prin2('vert3=*',vert3,3)
c        write(*,*)
c
        return
        end
c
c
c
c
c
        subroutine framerotx(vert1,vert2,vert3,alpha)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine rotates vertices vert1, vert2, and vert3 in R^3
c       about the x axis by the angle alpha
c
c       Input parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of the vertices
c       alpha   - the angle of the rotation 
c
c       Output parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of the rotated
c                                     vertices
c
        call framepntr(vert1(3),vert1(2),alpha)
        call framepntr(vert2(3),vert2(2),alpha)
        call framepntr(vert3(3),vert3(2),alpha)
c
        return
        end
c
c
c
c
c
        subroutine frameroty(vert1,vert2,vert3,alpha)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine rotates vertices vert1, vert2, and vert3 in R^3
c       about the y axis by the angle alpha
c
        call framepntr(vert1(1),vert1(3),alpha)
        call framepntr(vert2(1),vert2(3),alpha)
        call framepntr(vert3(1),vert3(3),alpha)
c
        return
        end
c
c
c
c
c
        subroutine framerotz(vert1,vert2,vert3,alpha)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine rotates vertices vert1, vert2, and vert3 in R^3
c       about the z axis by the angle alpha
c
c       Input parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of the vertices
c       alpha   - the angle of the rotation 
c
c       Output parameters:
c
c       vert1,vert2,vert3 (real *8) - the x,y,z coordinates of the rotated
c                                     vertices
c
        call framepntr(vert1(1),vert1(2),alpha)
        call framepntr(vert2(1),vert2(2),alpha)
        call framepntr(vert3(1),vert3(2),alpha)
c
        return
        end
c
c
c
c
c
        subroutine framepntr(x,y,alpha)
        implicit real *8 (a-h,o-z)
c
        cosa=cos(alpha)
        sina=sin(alpha)
c
        u=+x*cosa+y*sina
        v=-x*sina+y*cosa
c
        x=u
        y=v
c
        return
        end
c
c
c
c
c
        subroutine pointrotf(vert1,alpha,beta,gamma)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine performs the forward euler rotation of the
c       point in R^3
c
c       The rotations are performed in the following order (all in the
c       original coordinate system!)
c
c       First, perform a rotation around the z axis by the angle gamma
c       Then, perform a rotation around the y axis by the angle beta 
c       Finally, perform a rotation around the z axis by the angle alpha
c       
c
c       Input parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the point
c                     (these are changed on output)
c
c       Euler angles:
c       alpha   - the angle of the first rotation about z axis
c       beta    - the angle of the rotation about y axis
c       gamma   - the angle of the second rotation about z axis
c
c       Output parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the rotated point
c
c
        call pointrotz(vert1,gamma)
c        call prin2('vert1=*',vert1,3)
c        write(*,*)
        call pointroty(vert1,beta)
c        call prin2('vert1=*',vert1,3)
c        write(*,*)
        call pointrotz(vert1,alpha)
c        call prin2('vert1=*',vert1,3)
c        write(*,*)
c
        return
        end
c
c
c
c
c
        subroutine pointrotb(vert1,alpha,beta,gamma)
        implicit real *8 (a-h,o-z)
        dimension vert1(3)
c
c       This subroutine performs the backward euler rotation of the
c       point in R^3
c
c       The rotations are performed in the following order (all in the
c       original coordinate system!)
c
c       First, perform a rotation around the z axis by the angle -alpha 
c       Then, perform a rotation around the y axis by the angle -beta 
c       Finally, perform a rotation around the z axis by the angle -gamma 
c
c       Input parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of point
c                     (these are changed on output)
c
c       Euler angles:
c       alpha   - the angle of the first rotation about z axis
c       beta    - the angle of the rotation about y axis
c       gamma   - the angle of the second rotation about z axis
c
c       Output parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the rotated point
c
c
        call pointrotz(vert1,-alpha)
c        call prin2('vert1=*',vert1,3)
c        write(*,*)
        call pointroty(vert1,-beta)
c        call prin2('vert1=*',vert1,3)
c        write(*,*)
        call pointrotz(vert1,-gamma)
c        call prin2('vert1=*',vert1,3)
c        write(*,*)
c
        return
        end
c
c
c
c
c
        subroutine pointrotx(vert1,vert2,vert3,alpha)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine rotates point vert1 in R^3
c       about the x axis by the angle alpha
c
c       Input parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the point
c       alpha   - the angle of the rotation 
c
c       Output parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the rotated point
c                                     
c
        call pointpntr(vert1(3),vert1(2),alpha)
c
        return
        end
c
c
c
c
c
        subroutine pointroty(vert1,alpha)
        implicit real *8 (a-h,o-z)
        dimension vert1(3),vert2(3),vert3(3)
c
c       This subroutine rotates point vert1  in R^3
c       about the y axis by the angle alpha
c
        call pointpntr(vert1(1),vert1(3),alpha)
c
        return
        end
c
c
c
c
c
        subroutine pointrotz(vert1,alpha)
        implicit real *8 (a-h,o-z)
        dimension vert1(3)
c
c       This subroutine rotates point vert1 in R^3
c       about the z axis by the angle alpha
c
c       Input parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the point
c       alpha   - the angle of the rotation 
c
c       Output parameters:
c
c       vert1 (real *8) - the x,y,z coordinates of the rotated
c                                     point
c
        call pointpntr(vert1(1),vert1(2),alpha)
c
        return
        end
c
c
c
c
c
        subroutine pointpntr(x,y,alpha)
        implicit real *8 (a-h,o-z)
c
        cosa=cos(alpha)
        sina=sin(alpha)
c
        u=+x*cosa+y*sina
        v=-x*sina+y*cosa
c
        x=u
        y=v
c
        return
        end
c
c
c
c
c
