program main
    USE iso_fortran_env
    IMPLICIT NONE
    integer(INT32) , parameter :: n =100 ! number of nodes
    INTEGER(INT32) :: info , i
    integer(INT32) , dimension ( n ) :: ipiv
    REAL(REAL64) , dimension (n , n ) :: A ! An array of real numbers of size nxn
    REAL(REAL64) , dimension ( n ) :: x ,b,y,f,f_ana
    REAL(REAL64) :: h, G, nu, W, flux_num, flux_ana ! Mesh step
    REAL(REAL64) , dimension (n) :: absErr, relErr
    REAL(REAL64) :: avg_rel_error, rel_err_flux
    ! Distance between nodes , note we make ’n ’ the integer into a real number
    !using ’ real ’
    G=1.0_REAL64
    nu=0.000001_REAL64
    !nu=1.0_REAL64
    W=10.0_REAL64
    h = ((2.0_REAL64*W) /( REAL(n, kind=REAL64) -1.0_REAL64 ))
    
    ! Identify nodal positions y ( i ) and create analytic solution f_ana
    do i =1 , n
        y(i) = -W + ((i-1)*h)
        x ( i ) =( REAL( i, kind=REAL64 ) -1.0_REAL64 ) * h
!-----------------------!Analytic solution for velocity----------------
        f_ana ( i ) = (G/(2*nu))*((W**2_REAL64)-(y(i)**2.0_REAL64)) 
    end do

!-----------------Numerical Solution workflow----------------------------
    ! Initialise all values of A and b to be zero
    A =0.0_REAL64 ; b =0.0_REAL64

    do i =2 ,n -1
        
        A (i ,i -1) =1.0_REAL64 
        A (i , i ) = -2.0_REAL64 
        A (i , i +1) =1.0_REAL64 
        b(i) = ((-G*(h**2))/nu)
    end do
    ! Insert boundary values
    A (1 ,1) =1.0_REAL64
    b(1) = 0.0_REAL64
    A (n , n ) = 1.0_REAL64
    b(n) = 0.0_REAL64


    ! ‘ Call ’ the linear solver provided by Lapack
    ! Form of inputs and outputs can be found on Lapack website
    call dgesv (n ,1 ,A ,n , ipiv ,b ,n , info )
    ! If dgesv subroutine fails , then the integer ’ info ’ will be non - zero .
    ! Then the if statement will be true and the write statement will occur
    if ( info /= 0) then
        PRINT*,'Inversion Failed'
    else
        PRINT*,'Successful Inversion' ! info =0 so dgesv ran as planned
    end if
    ! dgesv outputs the solution to b , which we assign to f for neatness
    f=b
    ! Write the solution to a text file
    open (9 , file = "ana.txt" , form = 'formatted')

    23 FORMAT (3 ( ES23 .12 E3 ) )
    ! Choose format for text file ( real numbers )


    do i = 1, n
        absErr(i) = abs(f(i) - f_ana(i))
        if (abs(f_ana(i)) < 1.0e-9) then
          relErr(i) = 0.0_real64
        else
          relErr(i) = absErr(i)/f_ana(i)
        end if
    end do
    avg_rel_error = SUM(relErr)/n
    
    
    do i =1 , n
        write (9 ,23) y ( i ) , f_ana ( i ), f(i)
    end do
    close (9)
    flux_num = 0.0_REAL64 !initialise value of flux before summation
    flux_ana = ((2*(W**3)*G)/(3*nu))  !Calculte flux analytically
    DO i=2, N          
        flux_num = flux_num + (h * 0.5_REAL64 * (f(i)+f(i-1))) !calculate flux numerically
    END DO 
    !PRINT*,flux_ana
    !print*,flux1
    PRINT*,avg_rel_error
    rel_err_flux=(flux_num-flux_ana)/flux_ana
    print*,rel_err_flux
    !print*,(flux_num-flux_ana)
    !print*,flux_ana
    !print*,flux_num
end program main