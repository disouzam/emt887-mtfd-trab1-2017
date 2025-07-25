*******************************************************
* FILE: main.f
*******************************************************
      program Trab1

*Modules
        use AnalyticSolution2D
        use Coefficients
        use Geometry
        use Solver

*Settings
        implicit none

*Variable declaration
        logical :: debugmode
        common /dbgMode/ debugmode

        integer :: iter,iComp(1:4)
        common /control/ iter

        integer :: N(1:3), I, J, meth,relax_m, ADI
*       N: number of nodes in each direction
*          1,2,3 indices refer to x,y,z respectively

        double precision :: Lth(1:3)
*       Lth: Lenght along each direction
*          1,2,3 indices refer to x,y,z respectively

        double precision :: NPx(0:999), NIx(0:999), x_adim
* NPx: array containing node positions in X axis (to be calculated)
* NIx: array containing interface positions in X axis (to be calculated)

        double precision :: NPy(0:999), NIy(0:999), y_adim
* NPy: array containing node positions in Y axis (to be calculated)
* NIy: array containing interface positions in Y axis (to be calculated)

        double precision :: DelX(0:999), DeltaX(0:999)
* DelX: array containing distance between neighboring nodes in X axis(to be calculated)
* DeltaX: array containing distance between parallel interfaces of a volume element in X axis (to be calculated)

        double precision :: DelY(0:999), DeltaY(0:999)
* DelY: array containing distance between neighboring nodes in Y axis(to be calculated)
* DeltaY: array containing distance between parallel interfaces of a volume element in Y axis (to be calculated)

* Nodes: index ranges from 0 to N-1. There are Nx nodes
* Interfaces: index ranges from 0 to N. There is one interface more than the number of nodes
* Distances between nodes (DelX,DelY): index ranges from 0 to N-2. One less than the number of nodes
* Distance between interfaces (DeltaY,DeltaY): index ranges from 0 to N-1. The same number of nodes (but one less than the number of interfaces)

        common /Nodes/ NPx, NPy

* TitN: temperature at next iteration
* TitV: temperature at current iteration
        double precision :: tol
*       Tolerance value for iterative methods - it defines the stopping criteria
*       It is an absolute value, i.e. it is intented to be the maximum value
*         for the difference between temperatures at successive iterations

        double precision :: T1, T2
*       T1 and T2 are temperature values for a 2D permanent problem

        double precision :: Tbou(1:4)
*       Tbou stores the boundary temperatures
*       T(1): temperature at left boundary (West)
*       T(2): temperature at right boundary (East)
*       T(3): temperature at top boundary (North)
*       T(4): temperature at bottom boundary (South)

        double precision :: Aw(0:999,0:999), Ae(0:999,0:999)
        double precision :: An(0:999,0:999), As(0:999,0:999)
        double precision :: Ap(0:999,0:999)
*       Aw: coefficient for multiply the temperature of node at left (West)
*       Ae: coefficient for multiply the temperature of node at right (East)
*       An: coefficient for multiply the temperature of node at north (North)
*       As: coefficient for multiply the temperature of node at south (South)
*       Ap: coefficient for multiply the temperature of current node


        double precision :: Ta(0:999,0:999), Tb(0:999,0:999)
        double precision :: Tc(0:999,0:999), Td(0:999,0:999)
        double precision :: Tex(0:999,0:999)

        double precision :: alpha
*       A constant to be used as a factor for relaxation in Gauss-Seidel and TDMA methods


*       In debug mode (=.TRUE.), some info is displayed in screen
        debugmode = .FALSE.

*       Header of program
        print*,"==============================================================="
        print*,"        Program for calculation of a 2D thermal profile        "
        print*,"      in a rectangular plate in steady state conditions.       "
        print*,"                                                               "
        print*," Boundary conditions are constant temperature in borders       "
        print*,"      T1: temperatue at left, bottom and right borders         "
        print*,"      T2: temperatue at top border                             "
        print*,"==============================================================="
        print*," Author: Dickson Alves de Souza                                "
        print*," Based on lectures by professor Roberto Parreiras Tavares      "
        print*,"      and book Numerical Heat Transfer and Fluid Flow          "
        print*,"      by Suhas V. Patankar (1980)                              "
        print*,"                                                               "
        print*," Federal University of Minas Gerais                            "
        print*," September 21st, 2017                                          "
        print*,"==============================================================="
        print*,"                                                               "
        print*,"                                                               "
        print*,"                                                               "

* Instructions
        print*,"What are the dimensions of rectangular plate?"
        print*,"Lenght in X direction:"
        if (debugmode .EQV. .TRUE.) then
          Lth(1) = 2.5D0
          print*,"Lth(1) = ", Lth(1)
        else
          read*, Lth(1)
        end if

        print*
        print*,"Lenght in Y direction:"
        if (debugmode .EQV. .TRUE.) then
          Lth(2) = 2.5D0
          print*,"Lth(2) = ", Lth(2)
        else
          read*, Lth(2)
        end if
        print*
        print*
        print*,"What are the number of nodes in each direction?"
        print*,"Nodes in X direction:"
        if (debugmode .EQV. .TRUE.) then
          N(1) = 201
          print*,"N(1) = ", N(1)
        else
          read*, N(1)
        end if

        print*
        print*,"Nodes in Y direction:"
        if (debugmode .EQV. .TRUE.) then
          N(2) = 201
          print*,"N(2) = ", N(2)
        else
          read*, N(2)
        end if
        print*
        print*
        print*,"Specify the boundary conditions:"
        print*,"Temperature T1 (left, right and bottom borders):"
        if (debugmode .EQV. .TRUE.) then
          T1 = 1200.0D0
          print*,"T1 = ", T1
        else
          read*, T1
        end if

        print*
        print*,"Temperature T2 (top border):"
        if (debugmode .EQV. .TRUE.) then
          T2 = 700.0D0
          print*,"T2 = ", T2
        else
          read*, T2
        end if

        print*
        print*

        Tbou(1) = T1
        Tbou(2) = T1
        Tbou(3) = T2
        Tbou(4) = T1

        print*,"Select the solution method:"
        print*,"1 for Jacobi method"
        print*,"2 for Gauss-Seidel method"
        print*,"3 for TDMA (tri-diagonal matrix algorithm)"
        print*,"4 for comparison mode"
        print*
        if (debugmode .EQV. .TRUE.) then
          meth = 2
          print*,"meth = ", meth
        else
          read*, meth
        end if


        if (meth .EQ. 3) then
          print*,"Do you want to use ADI (alternating direction implicit) method?"
          print*,"1 for Yes"
          print*,"2 for No"
          print*
          if (debugmode .EQV. .TRUE.) then
            ADI = 2
            print*,"ADI = ", ADI
          else
            read*, ADI
          end if

        end if

        print*
        print*
        print*,"What tolerance should be considered?"
        print*,"Range: 1.0D-1 - 1.0D-8"
        if (debugmode .EQV. .TRUE.) then
          tol = 1.0D-5
          print*,"tol = ", tol
        else
          read*, tol

          if (tol .GT. 1.0D-1) then
            tol = 1.0D-1
          elseif (tol .LT. 1.0D-8) then
            tol = 1.0D-8
          end if
        end if

        print*
        print*
        print*,"Do you want to use a relaxation factor?"
        print*,"1 for Yes"
        print*,"2 for No"
        read*, relax_m
        print*

        if (relax_m .EQ. 1) then
          print*,"What relaxation factor do you want to use?"
          print*,"Range: 0.01 - 2.00"
          print*,"Be cautious: the chosen value could cause divergence of solution."
          print*
          read*,alpha

          if ((alpha .LT. 0.01) .OR. (alpha .GT. 2.0)) then
            print*,"You typed a wrong value for alpha. Default value (alpha = 1.0D0) will be considered, instead."
            alpha = 1.0D0
          end if

        else
          alpha = 1.0D0
        end if

*       Generate a 2D grid for solution
        call Grid2DUni(N,Lth,NPx,NIx,NPy,NIy,DelX,DeltaX,DelY,DeltaY)

*       Calculate coefficients for a 2D problem
        call coeff2D(N,DelX,DelY,Aw,Ae,An,As,Ap)

*       Select the solution method
        if (meth .EQ. 1) then
          call Jacobi(N,Aw,Ap,Ae,An,As,Ta,Tbou,tol,alpha)
          print*,"Jacobi method has terminated - Iterations: ", iter
          print*,"=========================================================================="
        elseif (meth .EQ. 2) then
          call GaSe2D(N,Aw,Ap,Ae,An,As,Ta,Tbou,tol,alpha)
          print*,"Gauss-Seidel method has terminated - Iterations: ", iter
          print*,"=========================================================================="
        elseif (meth .EQ. 3) then
          if (ADI .EQ. 1) then
            call TdmaAdi2D(N,Aw,Ap,Ae,An,As,Ta,Tbou,tol,alpha)
            print*,"TDMA method with ADI has terminated - Iterations: ", iter
            print*,"=========================================================================="
          elseif (ADI .EQ. 2) then
            call TDMA2D(N,Aw,Ap,Ae,An,As,Ta,Tbou,tol,alpha)
            print*,"TDMA method without ADI has terminated - Iterations: ", iter
            print*,"=========================================================================="
          end if
        elseif (meth .EQ. 4) then
          call Jacobi(N,Aw,Ap,Ae,An,As,Ta,Tbou,tol,alpha)
          iComp(1) = iter

          print*,"Jacobi method has terminated - Iterations: ", iComp(1)
          print*,"=========================================================================="
          do I = 0, 10000000, 1
            if(mod(I,2000000) .EQ. 0) then
              print*,""
            end if
          end do


          call GaSe2D(N,Aw,Ap,Ae,An,As,Tb,Tbou,tol,alpha)
          iComp(2) = iter

          print*,"Gauss-Seidel method has terminated - Iterations: ", iComp(2)
          print*,"=========================================================================="
          do I = 0, 10000000, 1
            if(mod(I,2000000) .EQ. 0) then
              print*,""
            end if
          end do


          call TDMA2D(N,Aw,Ap,Ae,An,As,Tc,Tbou,tol,alpha)
          iComp(3) = iter

          print*,"TDMA method without ADI has terminated - Iterations: ", iComp(3)
          print*,"=========================================================================="
          do I = 0, 10000000, 1
            if(mod(I,2000000) .EQ. 0) then
              print*,""
            end if
          end do


          call TdmaAdi2D(N,Aw,Ap,Ae,An,As,Td,Tbou,tol,alpha)
          iComp(4) = iter

          print*,"TDMA method with ADI has terminated - Iterations: ", iComp(3)
          print*,"=========================================================================="
          do I = 0, 10000000, 1
            if(mod(I,2000000) .EQ. 0) then
              print*,""
            end if
          end do
        end if

        call c_exSol2D(N,Tex,NPx,NPy,Tbou)

        print*,"========================================================================================"
        print*,"                          CALCULATION WAS DONE SUCCESSFULLY."
        print*,"========================================================================================"
        print*,"Absolute temperatures will be written to RESULT.DAT in the same folder of this program. "
        print*,"Adimensional values will be written to RES_ADIM.DAT in the same folder of this program. "
        print*,"========================================================================================"
        print*,""
        print*,""
        print*,"Press any key to end this program."
        read*

*       Writing results in RESULT.DAT file, stored in the same folder as the program
100     format(' ', 2I12,4F12.4)
200     format(' ', 2A12,4A12)
300     format(' ', A, F12.4)
400     format(' ', A, I12)
500     format(' ', 2I12,7F12.4)
600     format(' ', 9A12)
700     format(' ', A, E12.4)

        open(10,file='RESULT.DAT',status='UNKNOWN')
        write(10,*)"==============================================================="
        write(10,*)"        Steady state solution for a 2D rectangular plate       "
        write(10,*)"            with borders at constant temperature               "
        write(10,*)"                    Absolute temperatures                      "
        write(10,*)"==============================================================="
        write(10,*)"                Author: Dickson Alves de Souza                 "
        write(10,*)"   Based on lectures by professor Roberto Parreiras Tavares    "
        write(10,*)"      and book Numerical Heat Transfer and Fluid Flow          "
        write(10,*)"      by Suhas V. Patankar (1980)                              "
        write(10,*)"                                                               "
        write(10,*)"              Federal University of Minas Gerais               "
        write(10,*)"                     September 21st, 2017                      "
        write(10,*)"==============================================================="
        write(10,*)"                      Input parameters:                        "
        write(10,*)""
        write(10,300)"Lenght in X direction: ",Lth(1)
        write(10,300)"Lenght in Y direction: ",Lth(2)
        write(10,400)"Nodes in X direction:  ",N(1)
        write(10,400)"Nodes in Y direction:  ",N(2)
        write(10,300)"Temperature T1 (left, right and bottom borders): ", T1
        write(10,300)"Temperature T2 (top border):                     ", T2
        write(10,*)""

        if (meth .EQ. 1) then
          write(10,*) "SOLUTION of LINEAR SYSTEM: Jacobi method"
        elseif (meth .EQ. 2) then
          write(10,*) "SOLUTION of LINEAR SYSTEM: Gauss-Seidel method"
        elseif (meth .EQ. 3) then

          write(10,*) "SOLUTION of LINEAR SYSTEM: TDMA (tri-diagonal matrix algorithm)"

          if (ADI .EQ. 1) then
            write(10,*) "     ADI approach employed in solution."
          elseif (ADI .EQ. 2) then
            write(10,*) "     ADI approach NOT employed in solution."
          end if

        elseif (meth .EQ. 4) then
           write(10,*) "SOLUTION of LINEAR SYSTEM: Comparison of four methods"
           write(10,*) "1 - Jacobi Method"
           write(10,*) "2 - Gauss-Seidel Method"
           write(10,*) "3 - TDMA method without ADI"
           write(10,*) "4 - TDMA method with ADI"
           write(10,*)""
        end if

        if (relax_m .EQ. 1) then
          write(10,*)"    Relaxation factor: ", alpha
        elseif (relax_m .EQ. 2) then
          write(10,*)"    No relaxation applied to solution."
        end if

        write(10,700)"Tolerance: ", tol
        write(10,*)""

        if (meth .NE. 4) then
          write(10,*)"Iterations: ", iter
        else
          write(10,400)"Iterations -  Jacobi Method:           ", iComp(1)
          write(10,400)"Iterations -  Gauss-Seidel Method:     ", iComp(2)
          write(10,400)"Iterations -  TDMA method without ADI: ", iComp(3)
          write(10,400)"Iterations -  TDMA method with ADI:    ", iComp(4)
          write(10,*)""
        end if

        if (meth .NE. 4) then
          write(10,*)"========================================================================"
          write(10,*)"                   Calculation Results                                  "
          write(10,*)"========================================================================"
        else
          write(10,*)"============================================================================================================"
          write(10,*)"                                           Calculation Results                                              "
          write(10,*)"============================================================================================================"
        end if

        if (meth .NE. 4) then
          write(10,200)"I","J","X","Y","T", "Exact Sol"
        else
          write(10,600)"I","J","X","Y","T_Jac","T_Gaus","T_TDMA","T_TDMA_ADI", "Exact Sol"
        end if

        if (meth .NE. 4) then
          write(10,*)"========================================================================"
        else
          write(10,*)"============================================================================================================"
        end if

        if (meth .NE. 4) then

          do I = 0, N(1) - 1, 1
            do J = 0, N(2) - 1, 1
              write(10,100)I,J,NPx(I),NPy(J),Ta(I,J),Tex(I,J)
            end do
          end do

        else

          do I = 0, N(1) - 1, 1
            do J = 0, N(2) - 1, 1
              write(10,500)I,J,NPx(I),NPy(J),Ta(I,J),Tb(I,J),Tc(I,J),Td(I,J),Tex(I,J)
            end do
          end do

        end if

        endfile 10
        close(10,status='KEEP')


*       Writing results in RES_ADIM.DAT file, stored in the same folder as the program
120     format(' ', 2I12,4F12.6)
220     format(' ', 2A12,4A12)
320     format(' ', A, F12.6)
420     format(' ', A, I12)
520     format(' ', 2I12,7F12.6)
620     format(' ', 9A12)
720     format(' ', A, E12.6)

        open(20,file='RES_ADIM.DAT',status='UNKNOWN')
        write(20,*)"==============================================================="
        write(20,*)"        Steady state solution for a 2D rectangular plate       "
        write(20,*)"            with borders at constant temperature               "
        write(20,*)"                Adimensional temperatures                      "
        write(20,*)"==============================================================="
        write(20,*)"                Author: Dickson Alves de Souza                 "
        write(20,*)"   Based on lectures by professor Roberto Parreiras Tavares    "
        write(20,*)"      and book Numerical Heat Transfer and Fluid Flow          "
        write(20,*)"      by Suhas V. Patankar (1980)                              "
        write(20,*)"                                                               "
        write(20,*)"              Federal University of Minas Gerais               "
        write(20,*)"                     September 21st, 2017                      "
        write(20,*)"==============================================================="
        write(20,*)"                      Input parameters:                        "
        write(20,*)""
        write(20,320)"Lenght in X direction: ",Lth(1)
        write(20,320)"Lenght in Y direction: ",Lth(2)
        write(20,420)"Nodes in X direction:  ",N(1)
        write(20,420)"Nodes in Y direction:  ",N(2)
        write(20,320)"Temperature T1 (left, right and bottom borders): ", T1
        write(20,320)"Temperature T2 (top border):                     ", T2
        write(20,*)""

        if (meth .EQ. 1) then
          write(20,*) "SOLUTION of LINEAR SYSTEM: Jacobi method"
        elseif (meth .EQ. 2) then
          write(20,*) "SOLUTION of LINEAR SYSTEM: Gauss-Seidel method"
        elseif (meth .EQ. 3) then

          write(20,*) "SOLUTION of LINEAR SYSTEM: TDMA (tri-diagonal matrix algorithm)"

          if (ADI .EQ. 1) then
            write(20,*) "     ADI approach employed in solution."
          elseif (ADI .EQ. 2) then
            write(20,*) "     ADI approach NOT employed in solution."
          end if

        elseif (meth .EQ. 4) then
           write(20,*) "SOLUTION of LINEAR SYSTEM: Comparison of four methods"
           write(20,*) "1 - Jacobi Method"
           write(20,*) "2 - Gauss-Seidel Method"
           write(20,*) "3 - TDMA method without ADI"
           write(20,*) "4 - TDMA method with ADI"
           write(20,*)""
        end if

        if (relax_m .EQ. 1) then
          write(20,*)"    Relaxation factor: ", alpha
        elseif (relax_m .EQ. 2) then
          write(20,*)"    No relaxation applied to solution."
        end if

        write(20,720)"Tolerance: ", tol
        write(20,*)""

        if (meth .NE. 4) then
          write(20,*)"Iterations: ", iter
        else
          write(20,420)"Iterations -  Jacobi Method:           ", iComp(1)
          write(20,420)"Iterations -  Gauss-Seidel Method:     ", iComp(2)
          write(20,420)"Iterations -  TDMA method without ADI: ", iComp(3)
          write(20,420)"Iterations -  TDMA method with ADI:    ", iComp(4)
          write(20,*)""
        end if

        if (meth .NE. 4) then
          write(20,*)"========================================================================"
          write(20,*)"                   Calculation Results                                  "
          write(20,*)"                   Adimensional Values                                  "
          write(20,*)"========================================================================"
        else
          write(20,*)"============================================================================================================"
          write(20,*)"                                           Calculation Results                                              "
          write(20,*)"                                           Adimensional Values                                              "
          write(20,*)"============================================================================================================"
        end if

        if (meth .NE. 4) then
          write(20,220)"I","J","X/L","Y/W","T", "Exact Sol"
        else
          write(20,620)"I","J","X/L","Y/W","T_Jac","T_Gaus_Sei","T_TDMA","T_TDMA_ADI", "Exact Sol"
        end if

        if (meth .NE. 4) then
          write(20,*)"========================================================================"
        else
          write(20,*)"============================================================================================================"
        end if

        if (meth .NE. 4) then
          Ta = (Ta - T1) /(T2 - T1)
          Tex =(Tex - T1) /(T2 - T1)

          do I = 0, N(1) - 1, 1
            do J = 0, N(2) - 1, 1
              x_adim = NPx(I)/NPx(N(1)-1)
              y_adim = NPy(J)/NPy(N(2)-1)

              if (Ta(I,J) .LT. 1.0D-4) then
                Ta(I,J) = 0.0D0
              end if

              if (Tex(I,J) .LT. 1.0D-4) then
                Tex(I,J) = 0.0D0
              end if

              write(20,120)I,J,x_adim,y_adim,Ta(I,J),Tex(I,J)
            end do
          end do

        else
          Ta = (Ta - T1) /(T2 - T1)
          Tb = (Tb - T1) /(T2 - T1)
          Tc = (Tc - T1) /(T2 - T1)
          Td = (Td - T1) /(T2 - T1)
          Tex =(Tex - T1) /(T2 - T1)

          do I = 0, N(1) - 1, 1
            do J = 0, N(2) - 1, 1


              if (Ta(I,J) .LT. 1.0D-6) then
                Ta(I,J) = 0.0D0
              end if
              if (Tb(I,J) .LT. 1.0D-6) then
                Tb(I,J) = 0.0D0
              end if
              if (Tc(I,J) .LT. 1.0D-6) then
                Tc(I,J) = 0.0D0
              end if
              if (Td(I,J) .LT. 1.0D-6) then
                Td(I,J) = 0.0D0
              end if
              if (Tex(I,J) .LT. 1.0D-6) then
                Tex(I,J) = 0.0D0
              end if


              x_adim = NPx(I)/NPx(N(1)-1)
              y_adim = NPy(J)/NPy(N(2)-1)
              write(20,520)I,J,x_adim,y_adim,Ta(I,J),Tb(I,J),Tc(I,J),Td(I,J),Tex(I,J)
            end do
          end do

        end if

        endfile 20
        close(20,status='KEEP')

      end program Trab1

