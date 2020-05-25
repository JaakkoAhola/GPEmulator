program gp_predict
  use m_util
  use m_gp
  use m_gp_dense
  use readIO, only : readArray, readFileDimensions

  implicit none

  class(BaseGP), allocatable :: gp

  NAMELIST /inputoutput/  &
                      inputFileTrainingData, & ! training data
                      gpTrainingDataFile, & ! gp data of emulator training
                      separator,  & ! separator of data file
                      debugFlag

  integer :: ioStatusCode
  integer u,i

  character(len=max_name_len) label
  character(len=*), parameter :: filename = '../train/out.gp'
  real(dp), dimension(:,:), allocatable :: x_p,x
  integer, dimension(:), allocatable :: obs_type,obs_type_p
  real(dp), dimension(:), allocatable :: t,t_p,lt
  integer :: prediction_sample_size, prediction_dimension
  integer :: n, input_dimension
  real(dp) :: meanx
  real(dp) :: stdx
  real(dp) :: meany
  real(dp) :: stdy
  real(dp) :: rmse
  real(dp) :: meanlt
  real(dp) :: stdlt

  ! prediction_sample_size = 2
  ! prediction_dimension = 7
  ! n = 497
  ! input_dimension = 7

  allocate(gp, source = DenseGP(filename))
  allocate(real(dp) :: x_p(prediction_sample_size,prediction_dimension))
  allocate(real(dp) :: t_p(prediction_sample_size))
  allocate(real(dp) :: t(n))
  allocate(real(dp) :: lt(n))
  allocate(real(dp) :: x(n,input_dimension))
  allocate(integer :: obs_type(n))
  allocate(integer :: obs_type_p(prediction_sample_size))

  open  (1,status='old',file='NAMELIST.nml', iostat = ioStatusCode)
  if ( ioStatusCode /= 0 ) stop "Error opening NAMELIST.nml file"
  read  (1, nml=inputoutput)
  close (1)

  call readFileDimensions( inputFilePredictionData, prediction_sample_size, prediction_dimension, separator, debugFlag)

  call readFileDimensions( inputFileTrainingData, input_sample_size, input_dimension, separator, debugFlag)

  !!! REFACTORING HERE
  open(newunit=u, file="./data/DATA_predict")

  read (u,*) (x_p(i,:), obs_type_p(i), t_p(i), i=1,prediction_sample_size)

  close(u)

  open(newunit=u, file="../train/data/DATA_TRAIN")

  read (u,*) (x(i,:), obs_type(i), t(i), i=1,n)

  close(u)

  meanx = mean(x(:,1),n)
  meany = mean(t,n)
  stdx  = std(x(:,1),meanx,n)
  stdy  = std(t,meany,n)

  lt = logistic_vector(t,n)  ! lt = identity_vector(t,n),  which is not used in matlab, currently matlab code is using identify fun, change for all logistic and inv_logistic_vector in fortrune code.
  meanlt = mean(lt,n)
  stdlt  = std(lt,meanlt,n)

  t_p = logistic_vector(t_p,prediction_sample_size)
  t_p = standardize(t_p,meanlt,stdlt,prediction_sample_size)

  x_p(:,1) = standardize(x_p(:,1),meanx,stdx,prediction_sample_size)

  meanx = mean(x(:,2),n)
  stdx  = std(x(:,2),meanx,n)
  x_p(:,2) = standardize(x_p(:,2),meanx,stdx,prediction_sample_size)

  meanx = mean(x(:,3),n)
  stdx  = std(x(:,3),meanx,n)
  x_p(:,3) = standardize(x_p(:,3),meanx,stdx,prediction_sample_size)

  meanx = mean(x(:,4),n)
  stdx  = std(x(:,4),meanx,n)
  x_p(:,4) = standardize(x_p(:,4),meanx,stdx,prediction_sample_size)

  meanx = mean(x(:,5),n)
  stdx  = std(x(:,5),meanx,n)
  x_p(:,5) = standardize(x_p(:,5),meanx,stdx,prediction_sample_size)

  meanx = mean(x(:,6),n)
  stdx  = std(x(:,6),meanx,n)
  x_p(:,6) = standardize(x_p(:,6),meanx,stdx,prediction_sample_size)

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
    do i = 1,prediction_sample_size
       prdct = gp%predict(x_p(i,:), 0)
       rmse = rmse + ( inv_logistic(unstandardize_s(prdct,meanlt,stdlt))-inv_logistic(unstandardize_s(t_p(i),meanlt,stdlt) ) )**2
       print *, x_p(i,:)
      print *, inv_logistic(unstandardize_s(prdct,meanlt,stdlt))
       print *, unstandardize_s(prdct,meanlt,stdlt)
       print *, inv_logistic(unstandardize_s(t_p(i),meanlt,stdlt ) )
    end do
    rmse = rmse / prediction_sample_size
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
