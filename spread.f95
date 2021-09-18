! Author: Alex Schmith, aschmith2019@my.fit.edu
! Course: CSE 4250, Fall 2021
! Project: Proj1, Spread of Epidemics

! This program compiles with gfortan version :
! GNU Fortran (Ubuntu 5.4.0-6ubuntu1~16.04.10) 5.4.0 20160609


program epidemics
    implicit none
    integer :: population, population_infected, n = 1
    real :: y = 1, alpha

    read *, population, alpha
    do while(population > 0 .and. alpha > 0)
        y = infected(alpha) * 100
        population_infected = 1 + y * population / 100 
        
        print *,"Case", n, ":", population_infected, int(y),"%"
        n = n + 1
        read *, population, alpha
    end do

    contains 
        !The inverse function of x = ye^y approximated using the Halley Method
        real function lambert_w(x) result(wn)
            implicit none 
            real :: x, w = 1, difference
            
            !Using the Halley Method until difference between W(x) and W(x+1) is less than .000001
            wn = w - (w * EXP(w) - x)/( EXP(w) * (w + 1) - ((w+2) * (w * EXP(w) - x))/(2 * w + 2)  )
            difference = wn - w


            do while ((ABS(difference))> 0.00000001)
                w = wn 
                wn = w - (w * EXP(w) - x)/( EXP(w) * (w + 1) - ((w+2) * (w * EXP(w) - x))/(2 * w + 2)  )
                difference = wn - w
            end do
        end function lambert_w

        !function that returns the percentage of people that are infected in a population 
        !given the population and the expected amount of people a person is in contact with (alpha)
        real function infected(alpha) result(y)
            implicit none
            real :: alpha, num
            num = -1 * alpha * EXP(-1 * alpha)
            y = 1 + lambert_w(num)/alpha
            

        end function infected

end program epidemics