*******************************************************
* FILE: Results.f
*******************************************************
      module Results
        implicit none

      contains

      subroutine print_res1D(Nx,Temp)
        implicit none

        integer :: Nx, I
        double precision :: Temp(0:999,0:0)

        print*, "========================================================"
        print*, "Results for 1D calculation will be printed out on screen"
        print*, "========================================================"
        print*,""
        print*,""

        do I = 0, Nx - 1, 1
          print*,"T(",I,",0) = ", Temp(I,0)
        end do

      end subroutine print_res1D




      subroutine print_res2D(name,N,Temp,NPx,NPy)
        implicit none

        integer :: N(1:3), I, J
        double precision :: Temp(0:999,0:999)
        double precision :: NPx(0:999), NPy(0:999)
        character(LEN=20) :: name

        print*, "========================================================"
        print*, "Results for 2D calculation will be printed out on screen"
        print*, "              ", name, "                                "
        print*, "========================================================"
        print*,""
        print*,""

100     FORMAT(' ', A,I3,A,F7.2)
200     FORMAT('',F7.2,A)

        do I = 0, N(1) - 1, 1
          write(*,200,advance='no')  NPx(I),"  "
        end do
        print*, ""

        do J = N(2) - 1, 0, -1
          write(*,100) "J = ", J, "  ============ Y = ", NPy(J)
          print*, ""

          do I = 0, N(1) - 1, 1
            write(*,200,advance='no')  Temp(I,J),"  "
          end do
          print*, ""
          print*, ""
        end do
      end subroutine print_res2D

      subroutine compSol2D(N,TGaSe,TTDMA)
        implicit none

        integer :: N(1:3), I, J
        integer :: Imax,Jmax
        double precision :: maxdiff, temp
        double precision :: TGaSe(0:999,0:999)
        double precision :: TTDMA(0:999,0:999)

        maxdiff = 0
        Imax = 0
        Jmax = 0
        do I = 1, N(1) - 1, 1
          do J = 1, N(2) - 1, 1
            temp = abs(TGaSe(I,J)-TTDMA(I,J))
            if (temp .GT. maxdiff) then
              maxdiff = temp
              Imax = I
              Jmax = J
            end if
          end do
        end do

        print*,"Max difference = ", maxdiff
        print*,"Position =", Imax,",",Jmax
        print*,"Temp Gauss Seidel b: ", TGaSe(Imax,Jmax)
        print*,"Temp TDMA: ", TTDMA(Imax,Jmax)

      end subroutine compSol2D



      subroutine res2Dmod(name,N,Temp,NPx,NPy)
        implicit none

        integer :: N(1:3), I, J, interN, totalN
        double precision :: Temp(0:999,0:999)
        double precision :: NPx(0:999), NPy(0:999)
        character(LEN=20) :: name

        interN = 0
        totalN = 11
        do while (totalN .NE. N(1))
          interN = interN + 1
          totalN = 11 + interN * 10
        end do

        print*, "=================================================================="
        print*, "  Partial Results for 2D calculation will be printed out on screen"
        print*, "                        ", name, "                                "
        print*, "=================================================================="
        print*,""
        print*,""

100     FORMAT(' ', A,I3,A,F7.2)
200     FORMAT('',F7.2,A)

        do I = 0, N(1) - 1, interN + 1
          write(*,200,advance='no')  NPx(I),"  "
        end do
        print*, ""

        do J = N(2) - 1, 0, - interN - 1
          write(*,100) "J = ", J, "  ============ Y = ", NPy(J)
          print*, ""

          do I = 0, N(1) - 1, interN + 1
            write(*,200,advance='no')  Temp(I,J),"  "
          end do
          print*, ""
          print*, ""
        end do
      end subroutine res2Dmod


      end module Results
