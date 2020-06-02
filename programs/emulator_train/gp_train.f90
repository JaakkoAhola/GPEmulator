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
    integer :: trainingSampleSize, trainingDimensionSize, designDimensionSize

    ! loop counter
    integer :: ind
    ! number of hyperparameters and noise parameters required
    integer :: ntheta, nnu

    integer :: ioStatusCode

    real(dp), dimension(:), allocatable :: tmpArray
    real(dp), dimension(:), allocatable :: theta, nu, responseVector, lbounds, ubounds
    real(dp), dimension(:,:), allocatable :: designMatrix
    real(dp), dimension(:,:), allocatable :: trainingDataMatrix

    integer, dimension(:), allocatable :: trainingDataObservationTypeVector

    real(dp), dimension(16) :: tmpTheta = (/0.9010, 0.9650, 0.6729, 3.5576, 4.7418, 1.2722, &
     4.0612, 0.5, 2.4, 4.3, 3.2, 1.5, 0.5, 2.4, 4.3, 3.2 /)

    logical :: optimize

    integer :: optimize_max_iter
    real(dp) :: optimize_ftol

    character(len=max_name_len) :: covariance_function = 'LINSQEXP'
    character(len=max_name_len) :: noise_model_name    = 'VAL'

    character(len=1000) :: trainingDataInputFile
    character(len=1000) :: trainedEmulatorGPFile
    character(len=1) :: separator
    logical :: debugFlag

    NAMELIST /inputoutput/  &
                          trainingDataInputFile, & ! training data
                          trainedEmulatorGPFile, & ! output of executable (trained emulator)
                          separator,  & ! separator of data file
                          debugFlag

    class(cov_fn), allocatable :: cf
    class(noise_model), allocatable :: nm

    !!!!!!!!!!!!!!!!! END OF DECLARATIONS

    ! READ NAMELIST
    open  (1,status='old',file='train.nml', iostat = ioStatusCode)
        if ( ioStatusCode /= 0 ) stop "Error opening train.nml file"
        read  (1, nml=inputoutput)
    close (1)

    ! set dimension integers
    call readFileDimensions( trainingDataInputFile, trainingSampleSize, trainingDimensionSize, separator, debugFlag)

    designDimensionSize = trainingDimensionSize-2

    ! allocate matrices & vectors according to dimensions
    allocate( real(dp) :: trainingDataMatrix( trainingSampleSize, trainingDimensionSize ) )
    allocate(real(dp) :: designMatrix(trainingSampleSize,designDimensionSize))
    allocate(integer :: trainingDataObservationTypeVector(trainingSampleSize))
    allocate(real(dp) :: responseVector(trainingSampleSize))
    allocate(real(dp) :: tmpArray(trainingSampleSize))

    ! assign array values
    call readArray(trainingDataInputFile, trainingDataMatrix, separator, debugFlag)

    designMatrix(:,:) = trainingDataMatrix(:, 1:designDimensionSize)
    trainingDataObservationTypeVector(:) = int(trainingDataMatrix(:, trainingDimensionSize-1 ))
    responseVector(:) = trainingDataMatrix(:, trainingDimensionSize)

    ! Transform response and design here
    responseVector = standardize(responseVector,trainingSampleSize)
    do ind=1,designDimensionSize
        tmpArray = standardize( designMatrix(:,ind), trainingSampleSize )
        designMatrix(:,ind) = tmpArray
    end do
    !!!!!!!!!!!!!!!!!!!!!!!! END OF READING INPUT DATA

    call string_to_cov_fn(covariance_function, cf)
    call string_to_noise_model(noise_model_name, nm)

    nnu = nm%nparams_required(designDimensionSize)
    ntheta = cf%ntheta_required(designDimensionSize)

    allocate(real(dp) :: nu(nnu))
    allocate(real(dp) :: theta(ntheta))
    allocate(real(dp) :: lbounds(nnu + ntheta))
    allocate(real(dp) :: ubounds(nnu + ntheta))

    nu = 0.001
    theta = tmpTheta(1:designDimensionSize+1)
    lbounds(1) = 0.001
    lbounds(2:) = 0.01
    ubounds(:) = 100.0
    ubounds(1) = 0.001

    optimize = .true.
    optimize_max_iter = 10000
    optimize_ftol = 1.0d-7


    if ( debugFlag ) then
        write(*,*) "trainingDataInputFile: ", trim(trainingDataInputFile)
        write(*,*) "trainedEmulatorGPFile: ", trim(trainedEmulatorGPFile)
        write(*,*) "rows: ", trainingSampleSize
        write(*,*) "columns: ", trainingDimensionSize
        print*, "trainingDataMatrix(1,:)", trainingDataMatrix(1,:)
        print*, "trainingDataMatrix(2,:)", trainingDataMatrix(2,:)
        print*, "shape designMatrix", shape(designMatrix)
        print*, "shape trainingDataMatrix", shape(trainingDataMatrix)
        print*, "100's design second variable", designMatrix(100, 2)
        print*, "100's obs type", trainingDataObservationTypeVector(100)
        print*, "100's response", responseVector(100)
        print*, "100's trainingDataMatrix", trainingDataMatrix(100,:)
        print*, "nnu", nnu
        print*, "ntheta", ntheta
        print*, "lbounds", lbounds
        print*, "ubounds", ubounds
        print*, "optimize_max_iter", optimize_max_iter
        print*, "optimize_ftol", optimize_ftol
    end if

    allocate(gp, source=DenseGP(nu, theta, designMatrix, trainingDataObservationTypeVector, responseVector, cf, nm))


    if (optimize) then
        call log_lik_optim(nnu + ntheta, gp, lbounds, ubounds, optimize_max_iter, optimize_ftol)
    else
        gp%theta = (/ 0.9010,0.9650,0.6729,3.5576,4.7418,1.2722,4.0612, 1.5 /)
        gp%nu = 0.001
        call gp%update_matrices
    end if

    print *, gp%nu,' and ', gp%theta

    call gp%write_out(trainedEmulatorGPFile)

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
