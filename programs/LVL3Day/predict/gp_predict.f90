program gp_predict
  use m_util
  use m_gp
  use m_gp_dense

  implicit none

  class(BaseGP), allocatable :: gp

  integer u,i

  character(len=max_name_len) label
  character(len=*), parameter :: filename = '../train/out.gp'
  real(dp), dimension(:,:), allocatable :: x_p,x
  integer, dimension(:), allocatable :: obs_type,obs_type_p
  real(dp), dimension(:), allocatable :: t,t_p,lt
  integer :: n_p, input_dimension_p
  integer :: n, input_dimension
  real(dp) :: meanx
  real(dp) :: stdx
  real(dp) :: meany
  real(dp) :: stdy
  real(dp) :: rmse
  real(dp) :: meanlt
  real(dp) :: stdlt

  n_p = 2
  input_dimension_p = 7
  n = 497
  input_dimension = 7

  allocate(gp, source = DenseGP(filename))
  allocate(real(dp) :: x_p(n_p,input_dimension_p))
  allocate(real(dp) :: t_p(n_p))
  allocate(real(dp) :: t(n))
  allocate(real(dp) :: lt(n))
  allocate(real(dp) :: x(n,input_dimension))
  allocate(integer :: obs_type(n))
  allocate(integer :: obs_type_p(n_p))



  open(newunit=u, file="./data/DATA_predict")

  read (u,*) (x_p(i,:), obs_type_p(i), t_p(i), i=1,n_p)

  close(u)

  open(newunit=u, file="../train/data/DATA_TRAIN")

  read (u,*) (x(i,:), obs_type(i), t(i), i=1,n)

  close(u)

  meanx = mean(x(:,1),n)
  meany = mean(t,n)
  stdx  = std(x(:,1),meanx,n)
  stdy  = std(t,meany,n)

  !  lt = logistic_vector(t,n)
  lt = t
  meanlt = mean(lt,n)
  stdlt  = std(lt,meanlt,n)

  !  t_p = logistic_vector(t_p,n_p)
  t_p = t_p
  t_p = standardize(t_p,meanlt,stdlt,n_p)

  x_p(:,1) = standardize(x_p(:,1),meanx,stdx,n_p)

  meanx = mean(x(:,2),n)
  stdx  = std(x(:,2),meanx,n)
  x_p(:,2) = standardize(x_p(:,2),meanx,stdx,n_p)

  meanx = mean(x(:,3),n)
  stdx  = std(x(:,3),meanx,n)
  x_p(:,3) = standardize(x_p(:,3),meanx,stdx,n_p)

  meanx = mean(x(:,4),n)
  stdx  = std(x(:,4),meanx,n)
  x_p(:,4) = standardize(x_p(:,4),meanx,stdx,n_p)

  meanx = mean(x(:,5),n)
  stdx  = std(x(:,5),meanx,n)
  x_p(:,5) = standardize(x_p(:,5),meanx,stdx,n_p)

  meanx = mean(x(:,6),n)
  stdx  = std(x(:,6),meanx,n)
  x_p(:,6) = standardize(x_p(:,6),meanx,stdx,n_p)

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
    do i = 1,n_p
       prdct = gp%predict(x_p(i,:), 0)
       rmse = rmse + ( inv_logistic(unstandardize_s(prdct,meanlt,stdlt))-inv_logistic(unstandardize_s(t_p(i),meanlt,stdlt) ) )**2
       print *, x_p(i,:)
      print *, inv_logistic(unstandardize_s(prdct,meanlt,stdlt))
       print *, unstandardize_s(prdct,meanlt,stdlt)
       print *, inv_logistic(unstandardize_s(t_p(i),meanlt,stdlt ) )
    end do
    rmse = rmse / n_p
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
