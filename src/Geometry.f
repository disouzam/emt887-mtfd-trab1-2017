*******************************************************
* FILE: Geometry.f
*******************************************************
      module Geometry
            implicit none

      contains
* BasicGrid
      subroutine Grid1DUni(Nx,Lx,NPx,NIx,DelX,DeltaX)
* Grid1DUni calculates the node positions in a one-dimensional setting
* Uni means mesh is uniformly spaced
* Interface is located at half distance from both neighbor nodes, except for
* boundary nodes, where node and interface coincides

* Nx: number of nodes in direction X
* Lx: length of direction X
* NPx: array containing node positions in X axis (to be calculated)
* NIx: array containing interface positions in X axis (to be calculated)
* DelX: array containing distance between neighboring nodes (to be calculated)
* DeltaX: array containing distance between parallel interfaces of a volume element in X axis (to be calculated)



* Nodes: index ranges from 0 to Nx-1. There are Nx nodes
* Interfaces: index ranges from 0 to Nx. There is one interface more than the number of nodes
* Distances between nodes (DelX): index ranges from 0 to Nx-2. One less than the number of nodes
* Distance between interfaces (DeltaX): index ranges from 0 to Nx-1. The same number of nodes (but one less than the number of interfaces)


          implicit none

          logical :: debugmode
          common /dbgMode/ debugmode

          integer :: Nx, I
          double precision :: NPx(0:999),NIx(0:999), Lx, DX
          double precision :: DelX(0:999), DeltaX(0:999)

          double precision :: NPxt(0:999)

*         Mesh for X direction
          if (Nx .GT. 1) then
            DX = Lx / float(Nx - 1)
          else
            goto 100
          end if

          NPx(0) = 0.0D0
          NPxt(0) = NPx(0)

          NIx(0) = 0.0D0
          DelX(0) = DX

*         Loop trough all nodes, except the last one at Nx - 1 index
          do I = 1, Nx - 2, 1
            DelX(I) = DX
            NPx(I) = NPx(I-1) + DX
            NPxt(I)= NPx(I)
            NIx(I) = (NPx(I) + NPx(I-1)) / 2.0D0
          end do

          NPx(Nx-1) = Lx
          NIx(Nx-1) = (NPx(Nx-1) + NPx(Nx-2)) / 2.0D0
          NIx(Nx) = Lx

          do I = 0, Nx - 1, 1
            DeltaX(I) = NIx(I+1) - NIx(I)
          end do

          if (debugmode) then
            call ShGeom1D(Nx,NPx,NIx,DelX,DeltaX)
          end if

100     continue

      end subroutine Grid1DUni




      subroutine Grid2DUni(N,Lth,NPx,NIx,NPy,NIy,DelX,DeltaX,DelY,DeltaY)
* Grid2DUni calculates the node positions in a bi-dimensional setting
* Uni means mesh is uniformly spaced

*     N: number of nodes in each direction
*        1,2,3 indices refer to x,y,z respectively

*     Lth: Lenght along each direction
*          1,2,3 indices refer to x,y,z respectively

* NPx: array containing node positions in X axis (to be calculated)
* NIx: array containing interface positions in X axis (to be calculated)

* NPy: array containing node positions in Y axis (to be calculated)
* NIy: array containing interface positions in Y axis (to be calculated)

* DelX: array containing distance between neighboring nodes in X axis(to be calculated)
* DeltaX: array containing distance between parallel interfaces of a volume element in X axis (to be calculated)

* DelY: array containing distance between neighboring nodes in Y axis(to be calculated)
* DeltaY: array containing distance between parallel interfaces of a volume element in Y axis (to be calculated)


* Nodes: index ranges from 0 to N-1. There are Nx nodes
* Interfaces: index ranges from 0 to N. There is one interface more than the number of nodes
* Distances between nodes (DelX,DelY): index ranges from 0 to N-2. One less than the number of nodes
* Distance between interfaces (DeltaY,DeltaY): index ranges from 0 to N-1. The same number of nodes (but one less than the number of interfaces)

        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: N(1:3), I
        double precision :: Lth(1:3)
        double precision :: NPx(0:999), NIx(0:999),DX
        double precision :: DelX(0:999), DeltaX(0:999)

        double precision :: NPy(0:999), NIy(0:999),DY
        double precision :: DelY(0:999), DeltaY(0:999)

*       Mesh for X direction
        if (N(1) .GT. 1) then
          DX = Lth(1) / float(N(1) - 1)
        else
          goto 100
        end if

        NPx(0) = 0
        NIx(0) = 0
        DelX(0) = DX

        do I = 1, N(1) - 2, 1
          DelX(I) = DX
          NPx(I) = NPx(I-1) + DX
          NIx(I) = (NPx(I) + NPx(I-1)) / 2.0D0
        end do

        NPx(N(1)-1) = Lth(1)
        NIx(N(1)-1) = (NPx(N(1)-1) + NPx(N(1)-2)) / 2.0D0
        NIx(N(1)) = Lth(1)

        do I=0, N(1) - 1, 1
          DeltaX(I) = NIx(I+1) - NIx(I)
        end do

*       Mesh for Y direction
        if (N(2) .GT. 1) then
          DY = Lth(2) / float(N(2) - 1)
        else
          goto 100
        end if

        NPy(0) = 0
        NIy(0) = 0
        DelY(0) = DY

        do I = 1, N(2) - 2, 1
          DelY(I) = DY
          NPy(I) = NPy(I-1) + DY
          NIy(I) = (NPy(I) + NPy(I-1)) / 2.0D0
        end do

        NPy(N(2)-1) = Lth(2)
        NIy(N(2)-1) = (NPy(N(2)-1) + NPy(N(2)-2)) / 2.0D0
        NIy(N(2)) = Lth(2)

        do I = 0, N(2) - 1, 1
          DeltaY(I) = NIy(I+1) - NIy(I)
        end do

100   continue
      end subroutine Grid2DUni




      subroutine Grid3DUni(N,Lth,NPx,NIx,NPy,NIy,NPz,NIz,DelX,DeltaX,DelY,DeltaY,DelZ,DeltaZ)

* Grid3DUni calculates the node positions in a tri-dimensional setting
* Uni means mesh is uniformly spaced

*     N: number of nodes in each direction
*        1,2,3 indices refer to x,y,z respectively

*     Lth: Lenght along each direction
*          1,2,3 indices refer to x,y,z respectively
* NPx: array containing node positions in X axis (to be calculated)
* NIx: array containing interface positions in X axis (to be calculated)

* NPy: array containing node positions in Y axis (to be calculated)
* NIy: array containing interface positions in Y axis (to be calculated)

* NPz: array containing node positions in Z axis (to be calculated)
* NIz: array containing interface positions in Z axis (to be calculated)

* DelX: array containing distance between neighboring nodes in X axis(to be calculated)
* DeltaX: array containing distance between parallel interfaces of a volume element in X axis (to be calculated)

* DelY: array containing distance between neighboring nodes in Y axis(to be calculated)
* DeltaY: array containing distance between parallel interfaces of a volume element in Y axis (to be calculated)

* DelZ: array containing distance between neighboring nodes in Y axis(to be calculated)
* DeltaZ: array containing distance between parallel interfaces of a volume element in Y axis (to be calculated)

* Nodes: index ranges from 0 to N-1. There are Nx nodes
* Interfaces: index ranges from 0 to N. There is one interface more than the number of nodes
* Distances between nodes (DelX,DelY,DelZ): index ranges from 0 to N-2. One less than the number of nodes
* Distance between interfaces (DeltaY,DeltaY,DeltaZ): index ranges from 0 to N-1. The same number of nodes (but one less than the number of interfaces)

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: N(1:3), I, J, K
        double precision :: Lth(1:3)
        double precision :: NPx(0:999), NIx(0:999),DX
        double precision :: DelX(0:999), DeltaX(0:999)

        double precision :: NPy(0:999), NIy(0:999),DY
        double precision :: DelY(0:999), DeltaY(0:999)

        double precision :: NPz(0:999), NIz(0:999),DZ
        double precision :: DelZ(0:999), DeltaZ(0:999)

        double precision :: X,Y,Z

*       Mesh for X direction
        if (N(1) .GT. 1) then
          DX = Lth(1) / float(N(1) - 1)
        else
          goto 100
        end if

        NPx(0) = 0
        NIx(0) = 0
        DelX(0) = DX

        do I = 1, N(1) - 2, 1
          DelX(I) = DX
          NPx(I) = NPx(I-1) + DX
          NIx(I) = (NPx(I) + NPx(I-1)) / 2.0D0
        end do

        NPx(N(1)-1) = Lth(1)
        NIx(N(1)-1) = (NPx(N(1)-1) + NPx(N(1)-2)) / 2.0D0
        NIx(N(1)) = Lth(1)

        do I = 0, N(1) - 1, 1
          DeltaX(I) = NIx(I+1) - NIx(I)
        end do

*       Mesh for Y direction
        if (N(2) .GT. 1) then
          DY = Lth(2) / float(N(2) - 1)
        else
          goto 100
        end if

        NPy(0) = 0
        NIy(0) = 0
        DelY(0) = DY

        do I = 1, N(2) - 2, 1
          DelY(I) = DY
          NPy(I) = NPy(I-1) + DY
          NIy(I) = (NPy(I) + NPy(I-1)) / 2.0D0
        end do

        NPy(N(2)-1) = Lth(2)
        NIy(N(2)-1) = (NPy(N(2)-1) + NPy(N(2)-2)) / 2.0D0
        NIy(N(2)) = Lth(2)

        do I = 0, N(2) - 1, 1
          DeltaY(I) = NIy(I+1) - NIy(I)
        end do

*       Mesh for Z direction
        if (N(3) .GT. 1) then
          DZ = Lth(3) / float(N(3) - 1)
        else
          goto 100
        end if

        NPz(0) = 0
        NIz(0) = 0
        DelZ(0) = DZ

        do I = 1, N(3) - 2, 1
          DelZ(I) = DZ
          NPz(I) = NPz(I-1) + DZ
          NIz(I) = (NPz(I) + NPz(I-1)) / 2.0D0
        end do

        NPz(N(2)-1) = Lth(2)
        NIz(N(2)-1) = (NPz(N(2)-1) + NPz(N(2)-2)) / 2.0D0
        NIz(N(2)) = Lth(2)

        do I = 0, N(2) - 1, 1
          DeltaZ(I) = NIz(I+1) - NIz(I)
        end do

500     FORMAT(' ', A,I3,A,I3,A,I3,A,F5.3,A,F5.3,A,F5.3,A)
        if (debugmode) then
          do I = 0, N(1)-1,1
            X = NPx(I)
            do J = 0, N(2)-1,1
              Y = NPy(J)
              do K=0,N(3)-1,1
                Z = NPz(K)
      write(*,500) "N(",I,",",J,",",K,")=(",X,",",Y,",",Z,")"
              end do
              write(*,*)
            end do
            write(*,*)
          end do
        end if

100   continue

      end subroutine Grid3DUni

* RefinedGrid





* Output
      subroutine ShGeom1D(Nx,NPx,NIx,DelX,DeltaX)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: Nx, I
        double precision :: NPx(0:999),NIx(0:999)
        double precision :: DelX(0:999), DeltaX(0:999)

100     FORMAT(' ', A,I1,A,F6.3)
200     FORMAT(' ', A,I2,A,F7.3)
300     FORMAT(' ', A,I3,A,F8.3)
        if (debugmode) then
          do I = 0, Nx-1,1
            if (I .LT. 10) then
              write(*,100) "NPx(", I,")=", NPx(I)
              write(*,100) "NIx(", I,")=", NIx(I)
              write(*,100) "DelX(",I,")=", DelX(I)
              write(*,100) "DeltaX(",I,")=", DeltaX(I)
            end if
            if ((I .GE. 10) .AND. (I .LT. 100)) then
              write(*,200) "NPx(", I,")=", NPx(I)
              write(*,200) "NIx(", I,")=", NIx(I)
              write(*,200) "DelX(",I,")=", DelX(I)
              write(*,200) "DeltaX(",I,")=", DeltaX(I)
            end if
            if (I .GE. 100) then
              write(*,300) "NPx(", I,")=", NPx(I)
              write(*,300) "NIx(", I,")=", NIx(I)
              write(*,300) "DelX(",I,")=", DelX(I)
              write(*,300) "DeltaX(",I,")=", DeltaX(I)
            end if

            write(*,*) ""
          end do
          write(*,300) "NIx(", Nx,")=", NIx(Nx)

        end if

      end subroutine ShGeom1D

      end module Geometry
