program main
    USE iso_fortran_env
    IMPLICIT NONE
    integer(INT32) , parameter :: n =5000 ! number of nodes
    INTEGER(INT32) :: info , i, j
    integer(INT32) , dimension ( n ) :: ipiv
    REAL(REAL64) , dimension (n , n ) :: A ! An array of real numbers of size nxn
    REAL(REAL64) , dimension ( n ) :: x ,b,y,f,f_ana,f_num_temp, f_num
    REAL(REAL64) , dimension(:,:), ALLOCATABLE :: time_vel,time_vel_ana
    REAL(REAL64) , dimension(:), ALLOCATABLE :: time, time_vel_ana_final,time_vel_final
    INTEGER :: leng,r ,d, uno 
    ! Note solution is now
    !contained in f and x is spatial coordinates
    REAL(REAL64) :: h, G, nu, W, flux, flux1,T, dT, uana1, uana2, uana3, pi ! Mesh step
    ! Distance between nodes , note we make ’n ’ the integer into a real number
    !using ’ real ’
    !uno = 1_INT32
    pi=3.1415926535_REAL64
    G=1.0_REAL64
    nu=1.0_REAL64
    W=1.0_REAL64
    h = 2.0_REAL64*W /(REAL(n, kind=REAL64) -1.0_REAL64 )
    dt = 0.1_REAL64*(h**2)/nu
    T=(5.0_REAL64*(w**2))/((pi**2)*nu)    
    !T=50.0_REAL64*dt
    leng = FLOOR(T/dt, kind=INT32)
    !leng=50.0_REAL64

    !ALLOCATE(time(leng+1))
    !time(1) = 0.0_REAL64

    !do i = 2,(leng+1)
    !    time(i)= 0.0_REAL64 +((i-1)*dt)
    !end do



    
   !  Identify nodal positions x ( i ) and create analytic solution f_ana
    do i =1 , n
        y(i) = -W + ((i-1)*h)
        x ( i ) =( (REAL( i, kind=REAL64 ) -1.0_REAL64 ) * h)
        
        !x ( i ) =( REAL( i, kind=REAL64 ) -1.0_REAL64 ) * h
        f_ana ( i ) = (w**2 * G/(2*nu))*(((1.0_REAL64-y(i)**2)/W**2)-((32.0_REAL64/pi**3)*&
        COS(pi*y(i)/(2.0_REAL64*w))*EXP(-pi**2*nu*t/(4.0_REAL64*w**2))))
    end do
    

    f_num_temp = 0.0_REAL64
    f_num = 0.0_REAL64
    DO j=2, leng
        DO i = 2, (n-1)
            f_num(i) = f_num_temp(i) + dt*(g+nu*((f_num_temp(i-1) + f_num_temp(i+1) - 2*f_num_temp(i))/h**2))
        END DO
        f_num_temp=f_num
    END DO



    ! Write the solution to a text file
    open (9 , file = "ana.txt" , form = 'formatted')
!
    23 FORMAT (3 ( ES23 .12 E3 ) )
!    ! Choose format for text file ( real numbers )
    do i =1 , n
        write (9 ,23) y ( i ) , f_num(i), f_ana(i)
    end do

    !print*,time
end program main




    
    !do i =1 , n
    !    y(i) = -W + ((i-1)*h)
    !    x ( i ) =( REAL( i, kind=REAL64 ) -1.0_REAL64 ) * h
    !    f_ana ( i ) = (w**2 * g/(2*nu))*(((1.0_REAL64 -x(i)**2)/W**2)-((32.0_REAL64/pi**3)* &
    !    COS(pi*x(i)/(2.0_REAL64*W))*EXP(-pi**2*nu*T/(4.0_REAL64*W**2))))
    !end do


!    ALLOCATE(time_vel(n,(leng+1)))
!    DO i =1,(leng+1)
!        time_vel(1,i) = 0.0_REAL64
!        time_vel(n,i) =0.0_REAL64
!    END DO!

!    DO i=1,n
!        time_vel(i,1) = 0.0_REAL64
!    END DO
    

    ! Initialise all values of A and b to be zero
    !A =0.0_REAL64 ; b =0.0_REAL64
    ! Fill the matrices A and b with known values by looping through the
    !bulk / non - boundary nodes
!    do i =2 ,n -1
        
!        A (i ,i -1) =1.0_REAL64 
!        A (i , i ) = -2.0_REAL64 
!        A (i , i +1) =1.0_REAL64 
!        b(i) = ((-G*(h**2))/nu)
!    end do
    ! Insert boundary values
!    A (1 ,1) =1.0_REAL64
!    A (n , n ) = 1.0_REAL64


    ! ‘ Call ’ the linear solver provided by Lapack
    ! Form of inputs and outputs can be found on Lapack website
!    call dgesv (n ,1 ,A ,n , ipiv ,b ,n , info )
    ! If dgesv subroutine fails , then the integer ’ info ’ will be non - zero .
    ! Then the if statement will be true and the write statement will occur
 !   if ( info /= 0) then
 !       PRINT*,'Inversion Failed'
 !   else
 !       PRINT*,'Successful Inversion' ! info =0 so dgesv ran as planned
 !   end if
    ! dgesv outputs the solution to b , which we assign to f for neatness
 !   f=b



!    DO r = 1,leng
!        DO d= 2,n-1
!!            time_vel(d,r+1) = time_vel(d,r) + dt*(G + ((nu*((time_vel(d-1,r)) +time_vel(d+1,r) &
 !           -(2.0_REAL64*time_vel(d,r))))/(h**2)))
!!        END DO
 !   END DO


  !  ALLOCATE(time_vel_final(n))
  !  DO i=1,n
  !      time_vel_final(i) = time_vel(i, leng+1)
  !  END DO





!    close (9)
!    flux = 0.0_REAL64
!    DO i=2, N
!        flux = flux + (h * 0.5_REAL64 * (f_ana(i)+f_ana(i-1)))
!        flux1 = flux1 + (h * 0.5_REAL64 * (f(i)+f(i-1)))
!    END DO 
!    print*,flux
!    print*,flux1
