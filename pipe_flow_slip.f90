program main
    USE iso_fortran_env
    IMPLICIT NONE
    integer(INT32) , parameter :: n =500 ! number of nodes
    INTEGER(INT32) :: info , i
    integer(INT32) , dimension ( n ) :: ipiv
    REAL(REAL64) , dimension (n , n ) :: A ! An array of real numbers of size nxn
    REAL(REAL64) , dimension ( n ) :: b,f,f_ana,r
    REAL(REAL64) , dimension (n) :: absErr, relErr
    REAL(REAL64) :: avg_rel_error, rel_err_flux
    ! Note solution is now
    !contained in f and x is spatial coordinates
    REAL(REAL64) :: h, G, nu, W, flux_num, flux_ana ,l, pi! Mesh step
    ! Distance between nodes , note we make ’n ’ the integer into a real number
    !using ’ real ’
    pi = 3.14159265359_REAL64
    l=0.0000014_REAL64
    G=3.0_REAL64
    !nu=0.00089_REAL64
    nu=0.000001_REAL64
    W=0.000000001_REAL64
    h = ((W) /( REAL(n, kind=REAL64) -1.0_REAL64 ))
    
!---------------! Identify nodal positions r ( i ) and create analytic solution f_ana

    DO i =1 , n
        r(i)=0+(i-1)*h
    END DO
    DO i=1,n
        f_ana(i) = ((G/(4.0_REAL64*nu))*(W**2-(r(i)**2))) +((G*l*W)/(2.0_REAL64*nu))
    END DO

    ! Initialise all values of A and b to be zero
    !--------------Calc velocity numerically--------------------------------------
    A =0.0_REAL64 ; b =0.0_REAL64
    ! Fill the matrices A and b with known values by looping through the
    !bulk / non - boundary nodes
    do i =2 ,n-1 
        
        A (i ,i -1) = (-h+(2*r(i)))
        A (i , i ) = -4.0_REAL64 *r(i) 
        A (i , i +1) = (h+(2*r(i)))
        b(i) = ((-2*G*r(i)*(h**2))/nu)
    end do
    ! Insert boundary values
   !A (1 ,2) =A(2,1)                    
    A(1,1) = -1.0_REAL64
    A(1,2) = 1.0_REAl64
    A (n , n ) = ((-h/l)-1.0_REAL64)
    A(n,n-1) =1.0_REAL64

    ! ‘ Call ’ the linear solver provided by Lapack
    ! Form of inputs and outputs can be found on Lapack website
    call dgesv (n ,1 ,A ,n , ipiv ,b ,n , info )
    ! If dgesv subroutine fails , then the integer ’ info ’ will be non - zero .
    ! Then the if statement will be true and the write statement will occur
    if ( info /= 0) then
        PRINT*,'Inversion Failed'
        PRINT*,info
    else
        PRINT*,'Successful Inversion' ! info =0 so dgesv ran as planned
    end if
    ! dgesv outputs the solution to b , which we assign to f for neatness
    f=b
    ! Write the solution to a text file
    open (10 , file = "slp.txt" , form = 'formatted')

    23 FORMAT (3 ( ES23 .12 E3 ) )
    ! Choose format for text file ( real numbers )
    do i =1 , n
        write (10 ,23) r ( i ) , f(i), f_ana(i)
        
    end do

    
    close(10)

    do i = 1, n
        absErr(i) = abs(f(i) - f_ana(i))
        if (abs(f_ana(i)) < 1.0e-9) then
          relErr(i) = 0.0_real64
        else
          relErr(i) = absErr(i)/f_ana(i)
        end if
    end do
    avg_rel_error = SUM(relErr)/n

!-------------------------------------------------------------------------------------------
    flux_num = 0.0_REAL64
    flux_ana = ((pi*G*r(n)**3)/(8.0_REAL64*nu))*(((r(n)))+ (4.0_REAL64*l)) !calc flux analytically
    DO i=1, (N-1)
        flux_num = flux_num + ((r(i)*f(i))+(r(i+1)*f(i+1))) !calculate flux numerically
    END DO 
    flux_num = 2.0_REAL64*pi*(0.5_REAL64)*h*flux_num
    
    print*,flux_ana
    print*,flux_num
    !PRINT*,relErr
    !PRINT*,absErr
    !PRINT*,avg_rel_error
    rel_err_flux=(flux_num-flux_ana)/flux_ana
    !print*,rel_err_flux
    !PRINT*, shape(absErr)
    !PRINT*, r
    !PRINT*,absErr
end program main