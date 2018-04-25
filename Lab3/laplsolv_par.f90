program laplsolv
    use omp_lib
    !-----------------------------------------------------------------------
    ! Serial program for solving the heat conduction problem 
    ! on a square using the Jacobi method. 
    ! Written by Fredrik Berntsson (frber@math.liu.se) March 2003
    ! Modified by Berkant Savas (besav@math.liu.se) April 2006
    !-----------------------------------------------------------------------
    integer, parameter                  :: n=100, maxiter=1000
    double precision,parameter          :: tol=1.0E-3
    double precision,dimension(0:n+1,0:n+1) :: T
    double precision,dimension(n)       :: left, curr, right
    double precision                    :: error
    real                                :: t1,t0
    integer                             :: i,j,k,l_mainsize_block,from,to
    character(len=20)                   :: str

    ! Set boundary conditions and initial values for the unknowns
    T=0.0D0
    T(0:n+1 , 0)     = 1.0D0
    T(0:n+1 , n+1)   = 1.0D0
    T(n+1   , 0:n+1) = 2.0D0


    ! Solve the linear system of equations using the Jacobi method
    t0 = omp_get_wtime()

    do k=1,maxiter

        error = 0.0D0

        !$omp parallel private(j, size_block, from, to, left, curr, right) shared (T,k) reduction(max : error)

        size_block = CEILING(real(n) / omp_get_num_threads()) ! should us ceiling, but didn't work

        from = omp_get_thread_num() * size_block ! first column excluded
        to = from + size_block ! last column included
        if (to > n) then
            to = n
        end if
        
        left  = T(1:n, from) ! first column outside the section on the left
        right = T(1:n, to+1) ! first column outside the section on the right
        
        !$omp barrier
        
        do j=from+1,to
            curr=T(1:n,j) ! column currently working on
            if (j < to) then
                T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+T(1:n,j+1)+left)/4.0D0
            else ! if last column of section need to use data from next section saved at the start
                T(1:n,j)=(T(0:n-1,j)+T(2:n+1,j)+right+left)/4.0D0
            end if
            error=max(error,maxval(abs(curr-T(1:n,j))))
            left=curr
        end do

        !$omp end parallel 
        
        if (error<tol) then
            exit
        end if
        
    end do
    
    t1 = omp_get_wtime()

    write(unit=*,fmt=*) 'Time:',t1-t0,'Number of Iterations:',k
    write(unit=*,fmt=*) 'Temperature of element T(1,1)  =',T(1,1)

    ! Uncomment the next part if you want to write the whole solution
    ! to a file. Useful for plotting. 

    !open(unit=7,action='write',file='result.dat',status='unknown')
    !write(unit=str,fmt='(a,i6,a)') '(',N,'F10.6)'
    !do i=0,n+1
    !   write (unit=7,fmt=str) T(i,0:n+1)  
    !end do
    !close(unit=7)

end program laplsolv
