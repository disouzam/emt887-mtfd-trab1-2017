*******************************************************
* FILE: Coefficients.f
*******************************************************
      module Coefficients
        use Properties

        implicit none

      contains

      subroutine coeff1D(Nx,DelX,Aw,Ae,Ap)
* Calculates matrix coefficients for 1D problems
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: Nx, I
        double precision :: DelX(0:Nx-2)
        double precision :: Aw(0:999,0:0), Ae(0:999,0:0), Ap(0:999,0:0)

        do I = 1, Nx - 2, 1
          Aw(I,0) = k_T() / DelX(I-1)
          Ae(I,0) = k_T() / DelX(I)
          Ap(I,0) = Aw(I,0) + Ae(I,0)
        end do

      end subroutine coeff1D

      subroutine coeff2D(N,DelX,DelY,Aw,Ae,An,As,Ap)
* Calculates matrix coefficients for 2D problems
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: N(1:3), I, J
        double precision :: DelX(0:999), DelY(0:999)
        double precision :: Aw(0:999,0:999), Ae(0:999,0:999)
        double precision :: An(0:999,0:999), As(0:999,0:999)
        double precision :: Ap(0:999,0:999)

        do I = 1, N(1) - 2, 1
          do J = 1, N(2) - 2, 1
            Aw(I,J) = k_T() / DelX(I-1)
            Ae(I,J) = k_T() / DelX(I)
            An(I,J) = k_T() / DelY(J)
            As(I,J) = k_T() / DelY(J-1)
            Ap(I,J) = Aw(I,J) + Ae(I,J) + An(I,J) + As(I,J)
          end do
        end do

      end subroutine coeff2D

      end module Coefficients
