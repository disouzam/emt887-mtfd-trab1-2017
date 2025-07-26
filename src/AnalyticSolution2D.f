*******************************************************
* FILE: AnalyticSolution2D.f
*******************************************************
      module AnalyticSolution2D
            implicit none

      contains

      function PermSol2D(X,Y,L,W)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: I, seq
        double precision :: PermSol2D
        double precision :: PI, Iflo
        double precision :: X,Y,L,W

        double precision :: Ter1,Ter2,Ter3,Ter4,Temp,S
        double precision :: Ter3Arg, Ter4Arg

        if (Y .EQ. W) then
          PermSol2D = 1.0D0
        else

          PI =acos(-1.0D0)

          S = 0.0D0
          Temp = 1
          I=1

          seq = 0

          do while ((seq .LE. 10) .AND. (I .LT. 1.0E6))
              Iflo=float(I)

              Ter1=((-1.0D0)**(I+1)+1.0D0)/Iflo
              Ter2=sin(Iflo*PI*X/L)

              Ter3Arg = Iflo*PI*Y/L
              Ter4Arg = Iflo*PI*W/L

              if ((Ter3Arg .LT. 710) .AND. (Ter4Arg .LT. 710)) then
                Ter3=sinh(Ter3Arg)
                Ter4=sinh(Ter4Arg)
                Temp = Ter1*Ter2*Ter3/Ter4
                S = S + Temp
              else
                seq = 11
              end if

              I=I+2

              if (abs(temp) .LT. 1.0D-12) then
                seq = seq + 1
              end if

          end do

          S = float(2)/PI*S

          PermSol2D = S

        end if

      end function PermSol2D

      subroutine c_exSol2D(N,Texac,NPx,NPy,Tbou)

          implicit none

          integer :: N(1:3), I, J
          double precision :: Texac(0:999,0:999)
          double precision :: NPx(0:999), NPy(0:999)
          double precision :: Tbou(1:4)
          double precision :: L, W

          L = NPx(N(1)-1)
          W = NPy(N(2)-1)

          do I = 0, N(1) - 1, 1
            do J = 0, N(2) - 1, 1

              Texac(I,J) = PermSol2D(NPx(I), NPy(J),L,W) * (Tbou(3) - Tbou(4)) + Tbou(4)

            end do
          end do

        end subroutine c_exSol2D

      end module AnalyticSolution2D
