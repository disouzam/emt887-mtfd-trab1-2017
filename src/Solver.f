*******************************************************
* FILE: Solver.f
*******************************************************

      module Solver
        use Properties
        use AnalyticSolution2D

        implicit none

      contains

*     Gauss-Seidel
      subroutine GaSe1D(Nx,Aw,Ap,Ae,TitV,TitN,Tw,Te,tol)
        implicit none

        integer :: Nx, I, iter
        double precision :: maxdiff, tol, temp
        double precision :: Tw, Te
        double precision :: TitV(0:999,0:0), TitN(0:999,0:0)
        double precision :: Aw(0:999,0:0), Ae(0:999,0:0), Ap(0:999,0:0)

        iter = 0

*       Boundary conditions (constant temperature at boundaries)

*         Left border
        TitV(0,0) = Tw
        TitN(0,0) = Tw

*         Right border
        TitV(Nx - 1,0) = Te
        TitN(Nx - 1,0) = Te

*       Setting all inner node temperature to 0 (an initial guess to iterative process, not an initial condition)
        do I = 1, Nx - 2, 1
          TitV(I,0) = 0.0D0
          TitN(I,0) = 0.0D0
        end do

*       Implementation of Jacobi Algorithm
        temp = 1
        maxdiff = 0.1D0
        do while (maxdiff .GT. tol)

*         Solving linear system
          do I = 1, Nx - 2, 1
            TitN(I,0) = (Aw(I,0)*TitV(I-1,0) + Ae(I,0)*TitV(I+1,0)) / Ap(I,0)
          end do

*         Calculating stop criteria
          maxdiff = abs(TitN(1,0) - TitV(1,0))

          do I = 1, Nx - 2, 1
            temp = abs(TitN(I,0) - TitV(I,0))
            if (temp .GT. maxdiff) then
              maxdiff = temp
            end if
*           Updating TitV values from the newly calculated ones (TitN)
            TitV(I,0) = TitN(I,0)
          end do

          iter = iter + 1

        end do

      end subroutine GaSe1D




      subroutine Jacobi(N,Aw,Ap,Ae,An,As,TitN,Tbou,tol,alpha)
        implicit none

        integer :: N(1:3), I, J
        double precision :: maxdiff, tol, temp, total, alpha

        integer :: iter
        common /control/ iter

        double precision :: NPx(0:999)
        double precision :: NPy(0:999)
        common /Nodes/ NPx, NPy

        double precision :: Tbou(1:4)
*       Tbou stores the boundary temperatures
*       T(1): temperature at left boundary (West)
*       T(2): temperature at right boundary (East)
*       T(3): temperature at top boundary (North)
*       T(4): temperature at bottom boundary (South)


        double precision :: TitV(0:999,0:999)
        double precision :: TitN(0:999,0:999)
        double precision :: Aw(0:999,0:999), Ae(0:999,0:999)
        double precision :: An(0:999,0:999), As(0:999,0:999)
        double precision :: Ap(0:999,0:999)

        iter = 0

*       Boundary conditions

*       Setting all inner node temperature to 0 (an initial guess to iterative process, not an initial condition)
        TitV = 0.0D0
        TitN = 0.0D0

        do J = 0, N(2) - 2, 1
*         Setting the temperature in the left border (except top left node where temperature is equal to Tn)
          TitV(0,J) = Tbou(1)
          TitN(0,J) = Tbou(1)

*         Setting the temperature in the right border (except top left node where temperature is equal to Tn)
          TitV(N(1)-1,J) = Tbou(2)
          TitN(N(1)-1,J) = Tbou(2)
        end do

        do I = 0, N(1) - 1, 1
*         Setting the temperature in the bottom border
          TitV(I,0) = Tbou(4)
          TitN(I,0) = Tbou(4)

*         Setting the temperature in the top border
          TitV(I,N(2)-1) = Tbou(3)
          TitN(I,N(2)-1) = Tbou(3)
        end do

*       Implementation of Jacobi Algorithm
        temp = 1
        maxdiff = 100.0D0
        do while (maxdiff .GT. tol)

*         Solving linear system
          do I = 1, N(1) - 2 , 1
            do J = 1, N(2) - 2, 1
              total = Aw(I,J) * TitV(I-1,J)
              total = total + Ae(I,J) * TitV(I+1,J)
              total = total + As(I,J) * TitV(I,J-1)
              total = total + An(I,J) * TitV(I,J+1)

              TitN(I,J) = TitV(I,J) + alpha * (total / Ap(I,J) - TitV(I,J))
            end do
          end do

*         Calculating stop criteria
          maxdiff = abs(TitN(0,0) - TitV(0,0))

          do I = 0, N(1) - 1 , 1
            do J = 0, N(2) - 1, 1

              temp = abs(TitN(I,J) - TitV(I,J))

              if (temp .GT. maxdiff) then
                maxdiff = temp
              end if

*             Updating TitV values from the newly calculated ones (TitN)
              TitV(I,J) = TitN(I,J)
            end do
          end do

          iter = iter + 1

          if (mod(iter,100) .EQ. 0) then
            print*,"CONVERGENCE Monitor: Jacobi Method - Iterations = ", iter
            print*,"Maximum Difference = ", maxdiff
            print*,"Tolerance:", tol
            print*," "
          end if

        end do

      end subroutine Jacobi




      subroutine GaSe2D(N,Aw,Ap,Ae,An,As,TitN,Tbou,tol,alpha)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: iter
        common /control/ iter

        double precision :: NPx(0:999)
        double precision :: NPy(0:999)
        common /Nodes/ NPx, NPy

        integer :: N(1:3), I, J
        double precision :: maxdiff, tol, temp, total,alpha

        double precision :: Tbou(1:4)
*       Tbou stores the boundary temperatures
*       T(1): temperature at left boundary (West)
*       T(2): temperature at right boundary (East)
*       T(3): temperature at top boundary (North)
*       T(4): temperature at bottom boundary (South)


        double precision :: TitN(0:999,0:999)
        double precision :: Aw(0:999,0:999), Ae(0:999,0:999)
        double precision :: An(0:999,0:999), As(0:999,0:999)
        double precision :: Ap(0:999,0:999)

        iter = 0

*       Boundary conditions

*       Setting all inner node temperature to 0 (an initial guess to iterative process, not an initial condition)
        TitN = 0.0D0

        do J = 0, N(2) - 2, 1
*         Setting the temperature in the left border (except top left node where temperature is equal to Tn)
          TitN(0,J) = Tbou(1)

*         Setting the temperature in the right border (except top left node where temperature is equal to Tn)
          TitN(N(1)-1,J) = Tbou(2)
        end do

        do I = 0, N(1) - 1, 1
*         Setting the temperature in the bottom border
          TitN(I,0) = Tbou(4)

*         Setting the temperature in the top border
          TitN(I,N(2)-1) = Tbou(3)
        end do

*       Implementation of Gauss-Seidel Algorithm (properly)
        temp = 1
        maxdiff = 100.0D0
        do while (maxdiff .GT. tol)

*         Calculating stop criteria
          maxdiff = 0.0D0

*         Solving linear system
          do I = 1, N(1) - 2 , 1
            do J = 1, N(2) - 2, 1
              total = Aw(I,J) * TitN(I-1,J)
              total = total + Ae(I,J) * TitN(I+1,J)
              total = total + As(I,J) * TitN(I,J-1)
              total = total + An(I,J) * TitN(I,J+1)
              temp = TitN(I,J)

              if (alpha .NE. 1.0D0) then
                TitN(I,J) = temp + alpha * (total / Ap(I,J) - temp)
              else
                TitN(I,J) = total / Ap(I,J)
              end if

              temp = abs(TitN(I,J) - temp)
              if (temp .GT. maxdiff) then
                maxdiff = temp
              end if
            end do
          end do

*         Monitoring convergence
          iter = iter + 1

          if (mod(iter,100) .EQ. 0) then
            print*,"CONVERGENCE Monitor: Gauss-Seidel - Iterations = ", iter
            print*,"Maximum Difference = ", maxdiff
            print*,"Tolerance:", tol
            print*," "
          end if

        end do

      end subroutine GaSe2D




*     TDMA
      subroutine TDMA1D(Nx,Aw,Ap,Ae,TitV,TitN,Tw,Te)
        implicit none

        integer :: Nx, I, iter
        double precision :: Tw, Te
        double precision :: TitV(0:999,0:0), TitN(0:999,0:0)
        double precision :: Aw(0:999,0:0), Ae(0:999,0:0), Ap(0:999,0:0)

        double precision :: P(0:999),Q(0:999)

        iter = 0
        TitV(0,0) = Tw
        TitN(0,0) = Tw

        print*, "Tw = ", Tw

        do I = 1, Nx - 2, 1
          TitV(I,0) = 0.0D0
          TitN(I,0) = 0.0D0
        end do

        TitV(Nx - 1,0) = Te
        TitN(Nx - 1,0) = Te

        P(1) = Ae(1,0) / Ap(1,0)
        Q(1) = Aw(1,0) / Ap(1,0) * TitV(0,0)

        do I = 2, Nx - 2, 1
          P(I) = Ae(I,0) / (Ap(I,0) - Aw(I,0) * P(I-1))
          Q(I) = (Aw(I,0) * Q(I-1)) / (Ap(I,0) - Aw(I,0) * P(I-1))

        end do

        Q(Nx-1) =  TitN(Nx-1,0)

        TitN(Nx-1,0) = Q(Nx-1)

        do I = Nx-2, 1, -1
          TitN(I,0) = P(I) * TitN(I+1,0) + Q(I)
        end do

*       Not implemented the case where coefficients change with temperature (and hence there is need of iterations)

      end subroutine TDMA1D




      subroutine TDMA2D(N,Aw,Ap,Ae,An,As,TitN,Tbou,tol,alpha)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: iter
        common /control/ iter

        integer :: N(1:3), I, J
        double precision :: mxdifM, tol, temp, total,alpha

        double precision :: Tbou(1:4)
*       Tbou stores the boundary temperatures
*       T(1): temperature at left boundary (West)
*       T(2): temperature at right boundary (East)
*       T(3): temperature at top boundary (North)
*       T(4): temperature at bottom boundary (South)

        double precision :: TitN(0:999,0:999)


        double precision :: Aw(0:999,0:999), Ae(0:999,0:999)
        double precision :: An(0:999,0:999), As(0:999,0:999)
        double precision :: Ap(0:999,0:999)

*       Solve each horizontal lines (nodes along X axis are solved at each line)
*       Temperatures above and below the line being calculated are considered to be known
        double precision :: P(0:999),Q(0:999)

        double precision :: NPx(0:999)
        double precision :: NPy(0:999)
        common /Nodes/ NPx, NPy

*       Boundary conditions
        do J = 0, N(2) - 2, 1
*         Setting the temperature in the left border (except top left node where temperature is equal to Tn)
          TitN(0,J) = Tbou(1)

*         Setting the temperature in the right border (except top left node where temperature is equal to Tn)
          TitN(N(1)-1,J) = Tbou(2)
        end do


        do I = 0, N(1) - 1, 1
*         Setting the temperature in the bottom border
          TitN(I,0) = Tbou(4)

*         Setting the temperature in the top border
          TitN(I,N(2)-1) = Tbou(3)

        end do

*       Setting all inner node temperature to 0 (an initial guess to iterative process, not an initial condition)
        do I = 1, N(1) - 2, 1
          do J = 1, N(2) - 2, 1
            TitN(I,J) = 0.0D0
          end do
        end do


        iter = 0
        temp = 1
        mxdifM = 100.0D0
        do while (mxdifM .GT. tol)

          mxdifM = 0.0D0
          do J = 1, N(2) - 2, 1

*           define P1, Q1
            P(1) = Ae(1,J) / Ap(1,J)
            total = Aw(1,J) * TitN(0,J) + An(1,J) * TitN(1,J+1) + As(1,J) * TitN(1,J-1)
            Q(1) = total / Ap(1,J)

*           define Pi, Qi for i = 2 up to N(1) - 3
            do I = 2, N(1) - 3, 1
              P(I) = Ae(I,J) / (Ap(I,J) - Aw(I,J) * P(I-1))
              total = An(I,J) * TitN(I,J+1) + As(I,J) * TitN(I,J-1) + Aw(I,J) * Q(I-1)
              Q(I) = total / (Ap(I,J) - Aw(I,J) * P(I-1))
            end do

*           define Pi, Qi for i = N(1) - 2, the last unknown
            P(N(1)-2) = 0
            total = An(N(1)-2,J) * TitN(N(1)-2,J+1)
            total = total + As(N(1)-2,J) * TitN(N(1)-2,J-1)
            total = total + Ae(N(1)-2,J) * TitN(N(1)-1,J)
            total = total + Aw(N(1)-2,J) * Q(N(1)-3)
            Q(N(1)-2) = total / (Ap(N(1)-2,J) - Aw(N(1)-2,J) * P(N(1)-3))

*           Calculate TitN for the last unknown node
            temp = TitN(N(1)-2,J)

            TitN(N(1)-2,J) = temp + alpha *(Q(N(1)-2) - temp)

            temp = abs(TitN(N(1)-2,J) - temp)

            if (temp .GT. mxdifM) then
              mxdifM = temp
            end if

            do I = N(1) - 3, 1, -1
              temp = TitN(I,J)
              TitN(I,J) = temp + alpha * (P(I) * TitN(I+1,J) + Q(I) - temp)
              temp = abs(TitN(I,J)-temp)
              if (temp .GT. mxdifM) then
                mxdifM = temp
              end if
            end do
          end do


*         Monitoring convergence
          iter = iter + 1
          if (mod(iter,100) .EQ. 0) then
            print*,"CONVERGENCE Monitor: TDMA - Iterations = ", iter
            print*,"Maximum Difference = ", mxdifM
            print*,"Tolerance:", tol
            print*," "
          end if

        end do

      end subroutine TDMA2D





      subroutine TdmaAdi2D(N,Aw,Ap,Ae,An,As,TitN,Tbou,tol,alpha)
        implicit none

        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: iter
        common /control/ iter

        integer :: N(1:3), I, J
        double precision :: mxdifM, tol, temp, total,alpha
*        double precision :: Tw, Te, Tn, Ts

        double precision :: Tbou(1:4)
*       Tbou stores the boundary temperatures
*       T(1): temperature at left boundary (West)
*       T(2): temperature at right boundary (East)
*       T(3): temperature at top boundary (North)
*       T(4): temperature at bottom boundary (South)

        double precision :: TitN(0:999,0:999)
        double precision :: Aw(0:999,0:999), Ae(0:999,0:999)
        double precision :: An(0:999,0:999), As(0:999,0:999)
        double precision :: Ap(0:999,0:999)

*       Solve each horizontal lines (nodes along X axis are solved at each line)
*       Temperatures above and below the line being calculated are considered to be known
        double precision :: P(0:999),Q(0:999)

        double precision :: NPx(0:999)
        double precision :: NPy(0:999)
        common /Nodes/ NPx, NPy

*       Boundary conditions

*       Setting all inner node temperature to 0 (an initial guess to iterative process, not an initial condition)
        TitN = 0.0D0

        do J = 0, N(2) - 2, 1
*       Setting the temperature in the left border (except top left node where temperature is equal to Tn)
          TitN(0,J) = Tbou(1)

*       Setting the temperature in the right border (except top left node where temperature is equal to Tn)
          TitN(N(1)-1,J) = Tbou(2)
        end do

        do I = 0, N(1) - 1, 1
*       Setting the temperature in the bottom border
          TitN(I,0) = Tbou(4)

*       Setting the temperature in the top border
          TitN(I,N(2)-1) = Tbou(3)
        end do

        iter = 0
        temp = 1
        mxdifM = 100.0D0
        do while (mxdifM .GT. tol)

          mxdifM = 0.0D0


          if (mod(iter,4) .EQ. 0) then
*           Scan the domain from low Y (lower J) to high Y (higher J)
            do J = 1, N(2) - 2, 1

*             define P1, Q1
              P(1) = Ae(1,J) / Ap(1,J)

                total = Aw(1,J) * TitN(0,J) + An(1,J) * TitN(1,J+1) + As(1,J) * TitN(1,J-1)

              Q(1) = total / Ap(1,J)

*             define Pi, Qi for I = 2 up to N(1) - 3
              do I = 2, N(1) - 3, 1
                P(I) = Ae(I,J) / (Ap(I,J) - Aw(I,J) * P(I-1))

                  total = An(I,J) * TitN(I,J+1) + As(I,J) * TitN(I,J-1) + Aw(I,J) * Q(I-1)

                Q(I) = total / (Ap(I,J) - Aw(I,J) * P(I-1))
              end do

*             define Pi, Qi for I = N(1) - 2, the last unknown
              P(N(1)-2) = 0

                total = An(N(1)-2,J) * TitN(N(1)-2,J+1)
                total = total + As(N(1)-2,J) * TitN(N(1)-2,J-1)
                total = total + Ae(N(1)-2,J) * TitN(N(1)-1,J)
                total = total + Aw(N(1)-2,J) * Q(N(1)-3)

              Q(N(1)-2) = total / (Ap(N(1)-2,J) - Aw(N(1)-2,J) * P(N(1)-3))

*             Calculate TitN for the last unknown node
              temp = TitN(N(1)-2,J)

              TitN(N(1)-2,J) = temp + alpha * (Q(N(1)-2) - temp)

              temp = abs(TitN(N(1)-2,J) - temp)

              if (temp .GT. mxdifM) then
                mxdifM = temp
              end if

              do I = N(1) - 3, 1, -1
                temp = TitN(I,J)
                TitN(I,J) = temp + alpha* (P(I) * TitN(I+1,J) + Q(I) - temp)
                temp = abs(TitN(I,J)-temp)
                if (temp .GT. mxdifM) then
                  mxdifM = temp
                end if
              end do
            end do

          end if

          if (mod(iter,4) .EQ. 1) then
*           Scan the domain from high Y (higher J) to low Y (lower J)
            do J = N(2) - 2, 1, -1

*             define P1, Q1
              P(1) = Ae(1,J) / Ap(1,J)

                total = Aw(1,J) * TitN(0,J) + An(1,J) * TitN(1,J+1) + As(1,J) * TitN(1,J-1)

              Q(1) = total / Ap(1,J)

*             define Pi, Qi for I = 2 up to N(1) - 3
              do I = 2, N(1) - 3, 1
                P(I) = Ae(I,J) / (Ap(I,J) - Aw(I,J) * P(I-1))

                  total = An(I,J) * TitN(I,J+1) + As(I,J) * TitN(I,J-1) + Aw(I,J) * Q(I-1)

                Q(I) = total / (Ap(I,J) - Aw(I,J) * P(I-1))
              end do

*             define Pi, Qi for I = N(1) - 2, the last unknown
              P(N(1)-2) = 0

                total = An(N(1)-2,J) * TitN(N(1)-2,J+1)
                total = total + As(N(1)-2,J) * TitN(N(1)-2,J-1)
                total = total + Ae(N(1)-2,J) * TitN(N(1)-1,J)
                total = total + Aw(N(1)-2,J) * Q(N(1)-3)

              Q(N(1)-2) = total / (Ap(N(1)-2,J) - Aw(N(1)-2,J) * P(N(1)-3))

*             Calculate TitN for the last unknown node
              temp = TitN(N(1)-2,J)

              TitN(N(1)-2,J) = temp + alpha * (Q(N(1)-2) - temp)

              temp = abs(TitN(N(1)-2,J) - temp)

              if (temp .GT. mxdifM) then
                mxdifM = temp
              end if

              do I = N(1) - 3, 1, -1
                temp = TitN(I,J)
                TitN(I,J) = temp + alpha* (P(I) * TitN(I+1,J) + Q(I) - temp)
                temp = abs(TitN(I,J)-temp)
                if (temp .GT. mxdifM) then
                  mxdifM = temp
                end if
              end do
            end do

          end if

          if (mod(iter,4) .EQ. 2) then
*           Scan the domain from low X (lower I) to high X (higher I)
            do I = 1, N(1) - 2, 1

*             define P1, Q1
              P(1) = An(I,1) / Ap(I,1)

                total = As(I,1) * TitN(I,0) + Ae(I,1) * TitN(I+1,1) + Aw(I,1) * TitN(I-1,1)

              Q(1) = total / Ap(I,1)

*             define Pi, Qi for J = 2 up to N(2) - 3
              do J = 2, N(2) - 3, 1
                P(J) = An(I,J) / (Ap(I,J) - As(I,J) * P(J-1))

                  total = Ae(I,J) * TitN(I+1,J) + Aw(I,J) * TitN(I-1,J) + As(I,J) * Q(J-1)

                Q(J) = total / (Ap(I,J) - As(I,J) * P(J-1))
              end do

*             define Pi, Qi for J = N(2) - 2, the last unknown
              P(N(2)-2) = 0

                total = An(I,N(2)-2) * TitN(I,N(2)-1)
                total = total + Aw(I,N(2)-2) * TitN(I-1,N(2)-2)
                total = total + Ae(I,N(2)-2) * TitN(I+1,N(2)-2)
                total = total + As(I,N(2)-2) * Q(N(2)-3)

              Q(N(2)-2) = total / (Ap(I,N(2)-2) - As(I,N(2)-2) * P(N(2)-3))

*             Calculate TitN for the last unknown node
              temp = TitN(I,N(2)-2)

              TitN(I,N(2)-2) = temp + alpha * (Q(N(2)-2) - temp)

              temp = abs(TitN(I,N(2)-2) - temp)

              if (temp .GT. mxdifM) then
                mxdifM = temp
              end if

              do J = N(2) - 3, 1, -1
                temp = TitN(I,J)
                TitN(I,J) = temp + alpha*(P(J) * TitN(I,J+1) + Q(J) - temp)
                temp = abs(TitN(I,J)-temp)
                if (temp .GT. mxdifM) then
                  mxdifM = temp
                end if
              end do

            end do

          end if

          if (mod(iter,4) .EQ. 3) then
*          Scan the domain from high X (higher I) to low X (lower I)
            do I = N(1) - 2, 1, -1

*             define P1, Q1
              P(1) = An(I,1) / Ap(I,1)

                total = As(I,1) * TitN(I,0) + Ae(I,1) * TitN(I+1,1) + Aw(I,1) * TitN(I-1,1)

              Q(1) = total / Ap(I,1)

*             define Pi, Qi for J = 2 up to N(2) - 3
              do J = 2, N(2) - 3, 1
                P(J) = An(I,J) / (Ap(I,J) - As(I,J) * P(J-1))

                  total = Ae(I,J) * TitN(I+1,J) + Aw(I,J) * TitN(I-1,J) + As(I,J) * Q(J-1)

                Q(J) = total / (Ap(I,J) - As(I,J) * P(J-1))
              end do

*             define Pi, Qi for J = N(2) - 2, the last unknown
              P(N(2)-2) = 0

                total = An(I,N(2)-2) * TitN(I,N(2)-1)
                total = total + Aw(I,N(2)-2) * TitN(I-1,N(2)-2)
                total = total + Ae(I,N(2)-2) * TitN(I+1,N(2)-2)
                total = total + As(I,N(2)-2) * Q(N(2)-3)

              Q(N(2)-2) = total / (Ap(I,N(2)-2) - As(I,N(2)-2) * P(N(2)-3))

*             Calculate TitN for the last unknown node
              temp = TitN(I,N(2)-2)

              TitN(I,N(2)-2) = temp + alpha * (Q(N(2)-2) - temp)

              temp = abs(TitN(I,N(2)-2) - temp)

              if (temp .GT. mxdifM) then
                mxdifM = temp
              end if

              do J = N(2) - 3, 1, -1
                temp = TitN(I,J)
                TitN(I,J) = temp + alpha*(P(J) * TitN(I,J+1) + Q(J) - temp)
                temp = abs(TitN(I,J)-temp)
                if (temp .GT. mxdifM) then
                  mxdifM = temp
                end if
              end do

            end do

          end if

*         Monitoring convergence
          iter = iter + 1
          if (mod(iter,100) .EQ. 0) then
            print*,"CONVERGENCE Monitor: TDMA with ADI - Iterations = ", iter
            print*,"Maximum Difference = ", mxdifM
            print*,"Tolerance:", tol
            print*," "
          end if

        end do

      end subroutine TdmaAdi2D

      end module Solver
