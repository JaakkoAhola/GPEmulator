program gp_in
    use m_util
    use m_gp
    use m_gp_sparse
    use m_gp_dense
    use m_gp_optim
    use m_cov_all
    use m_noise_all
    use readIO, only : readArray, readFileDimensions

    implicit none

    class(BaseGP), allocatable :: gp

    ! number of training points, dimension of the input, dimension of design
    integer :: input_sample_size, input_dimension, input_design_dimension

    ! loop counter
    integer :: ind
    ! number of hyperparameters and noise parameters required
    integer :: ntheta, nnu

    integer :: ioStatusCode

    real(dp), dimension(:), allocatable :: tmpArray
    real(dp), dimension(:), allocatable :: theta, nu, response, lbounds, ubounds
    real(dp), dimension(:,:), allocatable :: design
    real(dp), dimension(:,:), allocatable :: array

    integer, dimension(:), allocatable :: obs_type

    real(dp), dimension(12) :: tmpTheta = (/0.9010, 0.9650, 0.6729, 3.5576, 4.7418, 1.2722, 4.0612, 1.5, 0.5, 2.4, 4.3, 3.2 /)

    logical :: optimize

    integer :: optimize_max_iter
    real(dp) :: optimize_ftol

    character(len=max_name_len) :: covariance_function = 'LINSQEXP'
    character(len=max_name_len) :: noise_model_name    = 'VAL'

    character(len=1000) :: inputfile
    character(len=1000) :: outputfile
    character(len=1) :: separator
    logical :: debugFlag

    NAMELIST /inputoutput/  &
                          inputfile, & ! training data
                          outputfile, & ! output of executable (trained emulator)
                          separator,  & ! separator of data file
                          debugFlag

    class(cov_fn), allocatable :: cf
    class(noise_model), allocatable :: nm


    open  (1,status='old',file='NAMELIST.nml', iostat = ioStatusCode)
    if ( ioStatusCode /= 0 ) stop "Error opening NAMELIST.nml file"
    read  (1, nml=inputoutput)
    close (1)

    write(*,*) "inputfile: ", trim(inputfile)
    write(*,*) "outputfile: ", trim(outputfile)

    call readFileDimensions( inputfile, input_sample_size, input_dimension, separator, debugFlag)

    write(*,*) "rows: ", input_sample_size
    write(*,*) "columns: ", input_dimension
    input_design_dimension = input_dimension-2

    allocate( real(dp) :: array( input_sample_size, input_dimension ) )

    call string_to_cov_fn(covariance_function, cf)
    call string_to_noise_model(noise_model_name, nm)

    nnu = nm%nparams_required(input_dimension)
    ntheta = cf%ntheta_required(input_dimension)

    allocate(real(dp) :: nu(nnu))
    allocate(real(dp) :: theta(ntheta))
    allocate(real(dp) :: lbounds(nnu + ntheta))
    allocate(real(dp) :: ubounds(nnu + ntheta))
    allocate(real(dp) :: response(input_sample_size))
    allocate(real(dp) :: design(input_sample_size,input_design_dimension))
    allocate(integer :: obs_type(input_sample_size))
    allocate(real(dp) :: tmpArray(input_sample_size))

    nu = 0.001
    theta = tmpTheta(1:input_design_dimension+1)
    if (debugFlag) print*, "theta", theta
    lbounds(1) = 0.001
    lbounds(2:) = 0.01
    ubounds(:) = 100.0
    ubounds(1) = 0.001

    optimize = .true.
    optimize_max_iter = 10000
    optimize_ftol = 1.0d-7

    call readArray(inputfile, array, separator, debugFlag)

    if (debugFlag) then
        print*, "array(1,:)", array(1,:)
        print*, "array(2,:)", array(2,:)
        print*, "shape design", shape(design)
        print*, "shape array", shape(array)
    end if

    design(:,:) = array(:, 1:input_design_dimension)
    obs_type(:) = int(array(:, input_dimension-1 ))
    response(:) = array(:, input_dimension)

    if ( debugFlag ) then
        print*, "100's design second variable", design(100, 2)
        print*, "100's obs type", obs_type(100)
        print*, "100's response", response(100)
        print*, "100's array", array(100,:)
    end if

    ! Transform design and response here
    response = standardize(response,input_sample_size)
    do ind=1,input_design_dimension
        tmpArray = standardize( design(:,ind), input_sample_size )
        design(:,ind) = tmpArray
    end do

    allocate(gp, source=DenseGP(nu, theta, design, obs_type, response, cf, nm))

    
    if (optimize) then
        call log_lik_optim(nnu + ntheta, gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
    else
        gp%theta = (/ 0.9010,0.9650,0.6729,3.5576,4.7418,1.2722,4.0612, 1.5 /)
        gp%nu = 0.001
        call gp%update_matrices
    end if

    print *, gp%nu,' and ', gp%theta

    call gp%write_out(outputfile)

contains
    function mean(x,dmn) result(res)
        integer dmn
        real(dp) x(dmn)
        real(dp) :: res

        res = sum(x)/dmn
    end function mean

    function std(x,meanx,dmn) result(res)
        integer dmn
        real(dp) x(dmn)
        real(dp) :: meanx
        real(dp) :: res

        res = sqrt(sum((x - meanx)**2)/(dmn-1))
    end function std

    function standardize(x,dmn) result(res)
        integer dmn
        real(dp) x(dmn)

        real(dp), dimension(dmn) :: res

        real(dp) :: meanx
        real(dp) :: stdx


        meanx= mean(x,dmn)

        stdx = std(x,meanx,dmn)

        res = (x - meanx)/stdx
    end function standardize

    function logistic(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res = 0.9 * x + 0.05;
        res = log(res) - log(1-res)
    end function logistic

    function logistic_vector(x,dmn) result(res)
        integer dmn
        real(dp)  x(dmn)
        real(dp),dimension(dmn) :: res
        res = 0.9 * x + 0.05
        res = log(res) - log(1-res)
    end function logistic_vector

    function inv_logistic(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res =  ( (1.0 / (1.0 + exp(-x))) - 0.05) / 0.9
    end function inv_logistic

    function inv_logistic_vector(x,dmn) result(res)
        integer dmn
        real(dp), intent(in),dimension(dmn) :: x
        real(dp),dimension(dmn) :: res
        res = ( (1.0 / (1.0 + exp(-x))) - 0.05) / 0.9
    end function inv_logistic_vector

end program gp_in
