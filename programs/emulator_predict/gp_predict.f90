program gp_predict
    use m_util
    use m_gp
    use m_gp_dense
    use readIO, only : readArray, readFileDimensions

    implicit none

    class(BaseGP), allocatable :: gp

    character(len=1000) :: trainingDataInputFile
    character(len=1000) :: trainedEmulatorGPFile
    character(len=1000) :: predictionDataInputFile
    character(len=1) :: separator
    logical :: debugFlag

    NAMELIST /inputoutput/  &
                      trainingDataInputFile, & ! training data
                      trainedEmulatorGPFile, & ! gp data of emulator training
                      predictionDataInputFile, & ! prediction data, i.e. points where you want your emulator to be used
                      separator,  & ! separator of data file
                      debugFlag

    integer :: ioStatusCode
    integer u,i

    real(dp), dimension(:,:), allocatable :: predictionDataMatrix,designMatrix
    integer, dimension(:), allocatable :: trainingDataObservationTypeVector,predictionDataObservationTypeVector
    real(dp), dimension(:), allocatable :: responseVector,predictionLastVector,lt
    integer :: predictionSampleSize, predictionDimensionSize
    integer :: trainingSampleSize, designDimensionSize
    real(dp) :: meanx
    real(dp) :: stdx
    real(dp) :: meany
    real(dp) :: stdy
    real(dp) :: rmse
    real(dp) :: meanlt
    real(dp) :: stdlt

    ! predictionSampleSize = 2
    ! predictionDimensionSize = 7
    ! trainingSampleSize = 497
    ! designDimensionSize = 7

    allocate(gp, source = DenseGP(trainedEmulatorGPFile))
    allocate(real(dp) :: predictionDataMatrix( predictionSampleSize,predictionDimensionSize ) )
    allocate(real(dp) :: predictionLastVector(predictionSampleSize))
    allocate(real(dp) :: responseVector(trainingSampleSize))
    allocate(real(dp) :: lt(trainingSampleSize))
    allocate(real(dp) :: designMatrix(trainingSampleSize,designDimensionSize))
    allocate(integer :: trainingDataObservationTypeVector(trainingSampleSize))
    allocate(integer :: predictionDataObservationTypeVector(predictionSampleSize))

    open  (1,status='old',file='NAMELIST.nml', iostat = ioStatusCode)
    if ( ioStatusCode /= 0 ) stop "Error opening NAMELIST.nml file"
    read  (1, nml=inputoutput)
    close (1)

    call readFileDimensions( predictionDataInputFile, predictionSampleSize, predictionDimensionSize, separator, debugFlag)

    call readFileDimensions( trainingDataInputFile, trainingSampleSize, designDimensionSize, separator, debugFlag)

    !!! REFACTORING HERE
    open(newunit=u, file=predictionDataInputFile)

    read (u,*) (predictionDataMatrix(i,:), predictionDataObservationTypeVector(i),&
      predictionLastVector(i), i=1,predictionSampleSize)

    close(u)

    open(newunit=u, file="../train/data/DATA_TRAIN")

    read (u,*) (designMatrix(i,:), trainingDataObservationTypeVector(i), responseVector(i), i=1,trainingSampleSize)

    close(u)

    meanx = mean(designMatrix(:,1),trainingSampleSize)
    meany = mean(responseVector,trainingSampleSize)
    stdx  = std(designMatrix(:,1),meanx,trainingSampleSize)
    stdy  = std(responseVector,meany,trainingSampleSize)

    lt = logistic_vector(responseVector,trainingSampleSize)
    meanlt = mean(lt,trainingSampleSize)
    stdlt  = std(lt,meanlt,trainingSampleSize)

    predictionLastVector = logistic_vector(predictionLastVector,predictionSampleSize)
    predictionLastVector = standardize(predictionLastVector,meanlt,stdlt,predictionSampleSize)

    predictionDataMatrix(:,1) = standardize(predictionDataMatrix(:,1),meanx,stdx,predictionSampleSize)

    meanx = mean(designMatrix(:,2),trainingSampleSize)
    stdx  = std(designMatrix(:,2),meanx,trainingSampleSize)
    predictionDataMatrix(:,2) = standardize(predictionDataMatrix(:,2),meanx,stdx,predictionSampleSize)

    meanx = mean(designMatrix(:,3),trainingSampleSize)
    stdx  = std(designMatrix(:,3),meanx,trainingSampleSize)
    predictionDataMatrix(:,3) = standardize(predictionDataMatrix(:,3),meanx,stdx,predictionSampleSize)

    meanx = mean(designMatrix(:,4),trainingSampleSize)
    stdx  = std(designMatrix(:,4),meanx,trainingSampleSize)
    predictionDataMatrix(:,4) = standardize(predictionDataMatrix(:,4),meanx,stdx,predictionSampleSize)

    meanx = mean(designMatrix(:,5),trainingSampleSize)
    stdx  = std(designMatrix(:,5),meanx,trainingSampleSize)
    predictionDataMatrix(:,5) = standardize(predictionDataMatrix(:,5),meanx,stdx,predictionSampleSize)

    meanx = mean(designMatrix(:,6),trainingSampleSize)
    stdx  = std(designMatrix(:,6),meanx,trainingSampleSize)
    predictionDataMatrix(:,6) = standardize(predictionDataMatrix(:,6),meanx,stdx,predictionSampleSize)

    rmse = predicttestset()
    print *, 'rmse: ',rmse

    print *, 'thetas: ',gp%theta

    print *, 'nu: ',gp%nu

    print *, 'design nrows: ',size(gp%x(:,1))


contains

    function predicttestset() result(rmse)
        real(dp) :: rmse
        real(dp) prdct
        rmse = 0
        do i = 1,predictionSampleSize
            prdct = gp%predict(predictionDataMatrix(i,:), 0)
            rmse = rmse + ( inv_logistic(unstandardize_s(prdct,meanlt,stdlt))&
            -inv_logistic(unstandardize_s(predictionLastVector(i),meanlt,stdlt) ) )**2
            print *, predictionDataMatrix(i,:)
            print *, inv_logistic(unstandardize_s(prdct,meanlt,stdlt))
            print *, unstandardize_s(prdct,meanlt,stdlt)
            print *, inv_logistic(unstandardize_s(predictionLastVector(i),meanlt,stdlt ) )
        end do
        rmse = rmse / predictionSampleSize
        rmse = sqrt(rmse)
    end function

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

    function standardize(x,meanx,stdx,dmn) result(res)
        integer dmn
        real(dp) x(dmn)

        real(dp), dimension(dmn) :: res

        real(dp) :: meanx
        real(dp) :: stdx


        res = (x - meanx)/stdx
    end function standardize

    function unstandardize(x,meanx,stdx,dmn) result(res)
        integer dmn
        real(dp) x(dmn)

        real(dp), dimension(dmn) :: res

        real(dp) :: meanx
        real(dp) :: stdx

        res = (x * stdx) + meanx
    end function unstandardize

    function unstandardize_s(x,meanx,stdx) result(res)
        real(dp) x

        real(dp) :: res

        real(dp) :: meanx
        real(dp) :: stdx


        res = (x * stdx) + meanx
    end function unstandardize_s

    function logistic(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res = 0.9_dp * x + 0.05_dp;
        res = log(res) - log(1-res)
    end function logistic

    function logistic_vector(x,dmn) result(res)
        integer dmn
        real(dp)  x(dmn)
        real(dp),dimension(dmn) :: res
        res = 0.9_dp * x + 0.05_dp
        res = log(res) - log(1-res)
    end function logistic_vector

    function inv_logistic(x) result(res)
        real(dp), intent(in) :: x
        real(dp) :: res
        res =  ( (1.0_dp / (1.0_dp + exp(-x))) - 0.05_dp) / 0.9_dp
    end function inv_logistic

    function inv_logistic_vector(x,dmn) result(res)
        integer dmn
        real(dp), intent(in),dimension(dmn) :: x
        real(dp),dimension(dmn) :: res
        res =  ( (1.0_dp / (1.0_dp + exp(-x))) - 0.05_dp) / 0.9_dp
    end function inv_logistic_vector


end program gp_predict
