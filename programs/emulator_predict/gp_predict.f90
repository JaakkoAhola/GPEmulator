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
    character(len=1000) :: predictionOutputFile
    character(len=1) :: separator
    logical :: debugFlag

    NAMELIST /inputoutput/  &
                      trainingDataInputFile, & ! training data
                      trainedEmulatorGPFile, & ! gp data of emulator training
                      predictionDataInputFile, & ! prediction input data, i.e. points where you want your emulator to be used
                      predictionOutputFile, & ! output of prediction
                      separator,  & ! separator of data file
                      debugFlag

    integer :: ioStatusCode
    integer u,i,ind

    real(dp), dimension(:,:), allocatable :: trainingDataMatrix, designMatrix, predictionDataMatrix,predictionSubMatrix
    integer, dimension(:), allocatable :: trainingDataObservationTypeVector,predictionDataObservationTypeVector
    real(dp), dimension(:), allocatable :: responseVector,predictionLastVector,logisticVector

    integer :: predictionSampleSize, predictionDimensionSize, predictionSubMatrixDimensionsSize
    integer :: trainingSampleSize, trainingDimensionSize, designDimensionSize
    real(dp) :: meanx
    real(dp) :: stdx
    real(dp) :: meany
    real(dp) :: stdy
    real(dp) :: rmse
    real(dp) :: meanlt
    real(dp) :: stdlt

    ! END OF DECLARATIONS

    ! predictionSampleSize = 26
    ! predictionDimensionSize = 6
    ! trainingSampleSize = 500
    ! designDimensionSize = 6

    ! READ NAMELIST
    open  (1,status='old',file='predict.nml', iostat = ioStatusCode)
        if ( ioStatusCode /= 0 ) stop "Error opening predict.nml file"
        read  (1, nml=inputoutput)
    close (1)

    allocate(gp, source = DenseGP(trainedEmulatorGPFile))

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! READ TRAINING DATA !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set dimension integers for training data
    call readFileDimensions( trainingDataInputFile, &
                            trainingSampleSize, trainingDimensionSize, &
                            separator, debugFlag)

    designDimensionSize = trainingDimensionSize-2

    ! allocate matrices & vectors according to dimensions for training data
    allocate( real(dp) :: trainingDataMatrix( trainingSampleSize, trainingDimensionSize ) )
    allocate(real(dp) :: designMatrix(trainingSampleSize,designDimensionSize))
    allocate(integer :: trainingDataObservationTypeVector(trainingSampleSize))
    allocate(real(dp) :: responseVector(trainingSampleSize))

    allocate(real(dp) :: logisticVector(trainingSampleSize))


    ! extract values from training data to suitable arrays
    call readArray(trainingDataInputFile, trainingDataMatrix, separator, debugFlag)

    designMatrix(:,:) = trainingDataMatrix(:, 1:designDimensionSize)
    trainingDataObservationTypeVector(:) = int(trainingDataMatrix(:, trainingDimensionSize-1 ))
    responseVector(:) = trainingDataMatrix(:, trainingDimensionSize)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!! READ PREDICTION DATA !!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! set dimension integers for prediction data
    call readFileDimensions( predictionDataInputFile, &
                            predictionSampleSize, predictionDimensionSize, &
                            separator, debugFlag)

    predictionSubMatrixDimensionsSize = predictionDimensionSize - 2

    ! allocate matrices & vectors according to dimensions for prediction data
    allocate(real(dp) :: predictionDataMatrix( predictionSampleSize,predictionDimensionSize ) )
    allocate(real(dp) :: predictionSubMatrix( predictionSampleSize,predictionSubMatrixDimensionsSize ) )
    allocate(integer :: predictionDataObservationTypeVector(predictionSampleSize))
    allocate(real(dp) :: predictionLastVector(predictionSampleSize))


    ! extract values from prediction data to suitable arrays
    call readArray(predictionDataInputFile, predictionDataMatrix, separator, debugFlag)

    predictionSubMatrix(:,:) = predictionDataMatrix(:, 1:predictionSubMatrixDimensionsSize)
    predictionDataObservationTypeVector(:) = int(predictionDataMatrix(:, predictionDimensionSize-1 ))
    predictionLastVector(:) = predictionDataMatrix(:, predictionDimensionSize)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! END OF READING DATA



    meany = mean(responseVector,trainingSampleSize)
    stdy  = std(responseVector,meany,trainingSampleSize)

    logisticVector = responseVector
    meanlt = mean(logisticVector,trainingSampleSize)
    stdlt  = std(logisticVector,meanlt,trainingSampleSize)

    predictionLastVector = logistic_vector(predictionLastVector,predictionSampleSize)
    predictionLastVector = standardize(predictionLastVector,meanlt,stdlt,predictionSampleSize)

    do ind=1,predictionSubMatrixDimensionsSize
        meanx = mean(designMatrix(:,ind),trainingSampleSize)
        stdx  = std(designMatrix(:,ind),meanx,trainingSampleSize)
        predictionSubMatrix(:,ind) = standardize(predictionSubMatrix(:,ind),meanx,stdx,predictionSampleSize)
    end do

    rmse = predicttestset(predictionOutputFile)
    print *, 'rmse: ',rmse

    print *, 'thetas: ',gp%theta

    print *, 'nu: ',gp%nu

    print *, 'design nrows: ',size(gp%x(:,1))


contains

    function predicttestset(predictionOutputFile) result(rmse)
        real(dp) :: rmse
        real(dp) prdct
        integer :: ioStatusCode
        character(len=1000) :: predictionOutputFile
        rmse = 0
        open  (14,status='new',file=predictionOutputFile, iostat = ioStatusCode)
            do i = 1,predictionSampleSize
                prdct = gp%predict(predictionSubMatrix(i,:), 0)
                rmse = rmse + ( inv_logistic(unstandardize_s(prdct,meanlt,stdlt))&
                -inv_logistic(unstandardize_s(predictionLastVector(i),meanlt,stdlt) ) )**2
                print *, predictionSubMatrix(i,:)
                print *, inv_logistic(unstandardize_s(prdct,meanlt,stdlt))
                print *, unstandardize_s(prdct,meanlt,stdlt)
                print *, inv_logistic(unstandardize_s(predictionLastVector(i),meanlt,stdlt ) )

                write(14,*) unstandardize_s(prdct,meanlt,stdlt)
            end do
        close(14)
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
