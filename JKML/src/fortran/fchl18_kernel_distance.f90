module fchl18_kernel_distance
   use ffchl_module, only: scalar, get_angular_norm2, get_pmax, get_ksi, init_cosp_sinp, get_selfscalar

   use ffchl_kernels, only: kernel
   implicit none

   private
   ! use target to declare as pointer-able (to avoid copying x1 and x2 a million times)
   !f2py intent(out) :: x1, x2
   double precision, dimension(:, :, :, :), allocatable :: X_tr, X_te
   ! Number of molecules
   integer :: num_molecules1, num_molecules2
   ! List of numbers of atoms in each molecule/cluster
   integer, dimension(:), allocatable :: num_atoms1, num_atoms2
   integer :: maxatoms1, maxatoms2
   ! Number of sigmas
   integer :: nsigmas ! save
   ! Number of neighbors
   integer, dimension(:, :), allocatable :: nneigh1, nneigh2
   ! Max neighbors
   integer :: maxneigh1, maxneigh2
   double precision :: t_width, d_width, ang_norm2
   ! Max index in the periodic table
   integer :: pmax1 ! can precalculate
   integer :: pmax2 ! can precalculate

   ! Fraction of cut_distance at which cut-off starts
   double precision :: cut_start, cut_distance
   ! Truncation order for Fourier terms
   integer :: order
   ! Periodic table distance matrix
   double precision, dimension(:, :), allocatable :: pd
   ! Scaling for angular and distance terms
   double precision :: distance_scale, angular_scale

   ! Switch alchemy on or off
   logical :: alchemy
   ! Decaying power laws for two- and three-body terms
   double precision :: two_body_power ! save
   double precision :: three_body_power ! save

   ! Kernel ID and corresponding parameters
   integer:: kernel_idx ! save
   double precision, dimension(:, :), allocatable :: parameters ! save

   logical :: verbose
   public :: init_train, init_test
   public :: train_kernel_element, test_kernel_element, test_train_kernel_element, test_distance, train_distance
   private :: X_tr, X_te
   private :: kernel_element

contains

   subroutine clean_train()
      ! clean up training data if it exists to avoid double alloc
      if (allocated(X_tr)) deallocate(X_tr)
      if (allocated(num_atoms1)) deallocate(num_atoms1)
      if (allocated(nneigh1)) deallocate(nneigh1)
      if (allocated(pd)) deallocate(pd)
      if (allocated(parameters)) deallocate(parameters)
   end subroutine clean_train

   subroutine init_train(x_in1, verbose_in, n1, nneigh1_in, nm1, nsigmas_in, &
   & t_width_in, d_width_in, cut_start_in, cut_distance_in, order_in, pd_in, &
   & distance_scale_in, angular_scale_in, alchemy_in, two_body_power_in, three_body_power_in, &
   & kernel_idx_in, parameters_in)
      ! fchl descriptors for the training set, format (i,maxatoms,5,maxneighbors)
      double precision, dimension(:, :, :, :), intent(in) :: x_in1 ! can be saved

      ! Whether to be verbose with output
      logical, intent(in) :: verbose_in

      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n1 ! save

      ! Number of molecules
      integer, intent(in) :: nm1 ! save

      ! Number of sigmas
      integer, intent(in) :: nsigmas_in ! save

      ! Number of neighbors for each atom in each compound
      integer, dimension(:, :), intent(in) :: nneigh1_in ! save

      ! Angular Gaussian width
      double precision, intent(in) :: t_width_in ! save

      ! Distance Gaussian width
      double precision, intent(in) :: d_width_in ! save

      ! Fraction of cut_distance at which cut-off starts
      double precision, intent(in) :: cut_start_in ! save
      double precision, intent(in) :: cut_distance_in ! save

      ! Truncation order for Fourier terms
      integer, intent(in) :: order_in ! save

      ! Periodic table distance matrix
      double precision, dimension(:, :), intent(in) :: pd_in ! save

      ! Scaling for angular and distance terms
      double precision, intent(in) :: distance_scale_in ! save
      double precision, intent(in) :: angular_scale_in ! save

      ! Switch alchemy on or off
      logical, intent(in) :: alchemy_in ! save

      ! Decaying power laws for two- and three-body terms
      double precision, intent(in) :: two_body_power_in ! save
      double precision, intent(in) :: three_body_power_in ! save

      ! Kernel ID and corresponding parameters
      integer, intent(in) :: kernel_idx_in ! save
      double precision, dimension(:, :), intent(in) :: parameters_in ! save

      ! clean up old data
      call clean_train()
      ! allocate raw datas
      allocate(X_tr(size(x_in1,1), size(x_in1,2), size(x_in1,3), size(x_in1,4)))
      X_tr = x_in1

      ! save sigmas
      nsigmas = nsigmas_in
      ! allocate atom number lists
      allocate(num_atoms1(size(n1,1)))
      num_atoms1 = n1

      maxatoms1 = maxval(n1)

      num_molecules1 = nm1

      ! allocate neighbour number lists
      allocate(nneigh1(size(nneigh1_in,1),size(nneigh1_in,2)))
      nneigh1 = nneigh1_in


      ! Get max number of neighbors
      maxneigh1 = maxval(nneigh1)

      t_width = t_width_in
      d_width = d_width_in
      ! Calculate angular normalization constant
      ang_norm2 = get_angular_norm2(t_width)

      ! pmax = max nuclear charge
      pmax1 = get_pmax(X_tr, n1)

      cut_start = cut_start_in
      cut_distance = cut_distance_in

      order = order_in

      allocate(pd(size(pd_in,1), size(pd_in,2)))
      pd = pd_in

      distance_scale = distance_scale_in
      angular_scale = angular_scale_in

      alchemy = alchemy_in

      two_body_power = two_body_power_in
      three_body_power = three_body_power_in

      kernel_idx = kernel_idx_in

      allocate(parameters(size(parameters_in,1),size(parameters_in,2)))
      parameters = parameters_in

      verbose = verbose_in
   end subroutine init_train

   subroutine clean_test()
      ! deallocate everything related to X_test
      if (allocated(X_te)) then
         deallocate(X_te)
      end if
      if (allocated(num_atoms2)) then
         deallocate(num_atoms2)
      end if
      if (allocated(nneigh2)) then
         deallocate(nneigh2)
      end if


   end subroutine clean_test

   subroutine init_test(x_in2, n2, nneigh2_in, nm2)
      ! fchl descriptors for the training set, format (i,maxatoms,5,maxneighbors)
      double precision, dimension(:, :, :, :), intent(in) :: x_in2 ! can be saved

      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n2 ! save

      ! Number of molecules
      integer, intent(in) :: nm2 ! save

      ! Number of neighbors for each atom in each compound
      integer, dimension(:, :), intent(in) :: nneigh2_in ! save

      ! clean up previously allocated values
      call clean_test()
      ! allocate raw datas
      allocate(X_te(size(x_in2,1), size(x_in2,2), size(x_in2,3), size(x_in2,4)))
      X_te = x_in2

      ! allocate atom number lists
      allocate(num_atoms2(size(n2,1)))
      num_atoms2 = n2

      maxatoms2 = maxval(n2)

      num_molecules2 = nm2

      ! allocate neighbour number lists
      allocate(nneigh2(size(nneigh2_in,1),size(nneigh2_in,2)))
      nneigh2 = nneigh2_in

      pmax2 = get_pmax(X_te, n2)

      ! Get max number of neighbors
      maxneigh2 = maxval(nneigh2)

   end subroutine init_test

   double precision function kernel_element(i, j, mode)

      ! Internal counters
      character(len=*), intent(in) :: mode
      integer, intent(in) :: i, j
      integer :: ii, jj

      ! Temporary variables necessary for parallelization
      double precision :: s12

      ! Pre-computed terms in the full distance matrix
      double precision, allocatable, dimension(:, :) :: self_scalar1 ! calculate on-the-fly
      double precision, allocatable, dimension(:, :) :: self_scalar2 ! calculate on-the-fly

      ! Pre-computed two-body weights
      double precision, allocatable, dimension(:, :, :) :: ksi1 ! calculate on-the-fly
      double precision, allocatable, dimension(:, :, :) :: ksi2 ! calculate on-the-fly

      ! Pre-computed terms for the Fourier expansion of the three-body term
      double precision, allocatable, dimension(:, :, :, :, :) :: sinp1 ! otf
      double precision, allocatable, dimension(:, :, :, :, :) :: sinp2 ! otf
      double precision, allocatable, dimension(:, :, :, :, :) :: cosp1 ! otf
      double precision, allocatable, dimension(:, :, :, :, :) :: cosp2 ! otf

      ! single elements from x1, x2
      double precision, dimension(:, :, :, :), allocatable :: v1
      double precision, dimension(:, :, :, :), allocatable :: v2

      integer :: n1(1)
      integer :: n2(1)
      integer, dimension(:, :), allocatable :: nneigh1_i
      integer, dimension(:, :), allocatable :: nneigh2_j
      integer :: mneigh1, mneigh2, pm1, pm2

      ! Work kernel space
      integer :: n
      double precision, allocatable, dimension(:) :: ktmp
      double precision, allocatable, dimension(:) :: kernels

      ! copy relevant elements to the pointers
      select case (trim(mode))
       case ("train")
         ! train x train kernel
         allocate(v1(1,size(X_tr,2),5,size(X_tr,4)))
         allocate(v2(1,size(X_tr,2),5,size(X_tr,4)))
         v1(1,:,:,:) = X_tr(i, :, :, :)
         v2(1,:,:,:) = X_tr(j, :, :, :)
         n1 = [num_atoms1(i)]
         n2 = [num_atoms1(j)]
         nneigh1_i = reshape(nneigh1(i, :), [1, size(nneigh1, 2)])
         nneigh2_j = reshape(nneigh1(j, :), [1, size(nneigh1, 2)])
         mneigh1 = maxval(nneigh1_i)
         mneigh2 = maxval(nneigh2_j)
         pm1 = pmax1
         pm2 = pmax1
         allocate (ksi1(1, num_atoms1(i), mneigh1))
         allocate (ksi2(1, num_atoms1(j), mneigh2))
         allocate (cosp1(1, num_atoms1(i), pmax1, order, mneigh1))
         allocate (sinp1(1, num_atoms1(i), pmax1, order, mneigh1))
         allocate (cosp2(1, num_atoms1(j), pmax1, order, mneigh2))
         allocate (sinp2(1, num_atoms1(j), pmax1, order, mneigh2))
         allocate (self_scalar1(1, num_atoms1(i)))
         allocate (self_scalar2(1, num_atoms1(j)))
       case ("test")
         ! test x test kernel
         allocate(v1(1,size(X_te,2),5,size(X_te,4)))
         allocate(v2(1,size(X_te,2),5,size(X_te,4)))
         v1(1,:,:,:) = X_te(i, :, :, :)
         v2(1,:,:,:) = X_te(j, :, :, :)
         n1 = [num_atoms2(i)]
         n2 = [num_atoms2(j)]
         nneigh1_i = reshape(nneigh2(i, :), [1, size(nneigh2, 2)])
         nneigh2_j = reshape(nneigh2(j, :), [1, size(nneigh2, 2)])
         mneigh1 = maxval(nneigh1_i)
         mneigh2 = maxval(nneigh2_j)
         pm1 = pmax2
         pm2 = pmax2
         allocate (ksi1(1, num_atoms2(i), mneigh1))
         allocate (ksi2(1, num_atoms2(j), mneigh2))
         allocate (cosp1(1, num_atoms2(i), pmax2, order, mneigh1))
         allocate (sinp1(1, num_atoms2(i), pmax2, order, mneigh1))
         allocate (cosp2(1, num_atoms2(j), pmax2, order, mneigh2))
         allocate (sinp2(1, num_atoms2(j), pmax2, order, mneigh2))
         allocate (self_scalar1(1, num_atoms2(i)))
         allocate (self_scalar2(1, num_atoms2(j)))
       case ("test_train")
         ! test x train kernel
         allocate(v1(1,size(X_te,2),5,size(X_te,4)))
         allocate(v2(1,size(X_tr,2),5,size(X_tr,4)))
         v1(1,:,:,:) = X_te(i, :, :, :)
         v2(1,:,:,:) = X_tr(j, :, :, :)
         ! admittedly confusing, but now 1 = 2
         n1 = [num_atoms2(i)]
         n2 = [num_atoms1(j)]
         nneigh1_i = reshape(nneigh2(i, :), [1, size(nneigh2, 2)])
         nneigh2_j = reshape(nneigh1(j, :), [1, size(nneigh1, 2)])
         mneigh1 = maxval(nneigh1_i)
         mneigh2 = maxval(nneigh2_j)
         pm1 = pmax2
         pm2 = pmax1
         allocate (ksi1(1, num_atoms2(i), mneigh1))
         allocate (ksi2(1, num_atoms1(j), mneigh2))
         allocate (cosp1(1, num_atoms2(i), pmax2, order, mneigh1))
         allocate (sinp1(1, num_atoms2(i), pmax2, order, mneigh1))
         allocate (cosp2(1, num_atoms1(j), pmax1, order, mneigh2))
         allocate (sinp2(1, num_atoms1(j), pmax1, order, mneigh2))
         allocate (self_scalar1(1, num_atoms2(i)))
         allocate (self_scalar2(1, num_atoms1(j)))
       case default
         print *, "Unknown kernel mode: ", mode
         stop
      end select

      ! compute ksis
      call get_ksi(v1, n1, nneigh1_i, two_body_power, cut_start, cut_distance, verbose, ksi1)
      call get_ksi(v2, n2, nneigh2_j, two_body_power, cut_start, cut_distance, verbose, ksi2)

      n = size(parameters, dim=1)
      allocate (ktmp(n))
      allocate(kernels(n))
      kernels = 0.0d0

      ! Initialize and pre-calculate three-body Fourier terms
      call init_cosp_sinp(v1, n1, nneigh1_i, three_body_power, order, cut_start, cut_distance, pm1, mneigh1, &
      & cosp1, sinp1, verbose)

      ! Initialize and pre-calculate three-body Fourier terms
      call init_cosp_sinp(v2, n2, nneigh2_j, three_body_power, order, cut_start, cut_distance, pm2, mneigh2, &
      & cosp2, sinp2, verbose)

      ! Pre-calculate self-scalar terms
      call get_selfscalar(v1, 1, n1, nneigh1_i, ksi1, sinp1, cosp1, t_width, d_width, &
      & cut_distance, order, pd, ang_norm2, distance_scale, angular_scale, alchemy, verbose,&
      &self_scalar1)

      ! Pre-calculate self-scalar terms
      call get_selfscalar(v2, 1, n2, nneigh2_j, ksi2, sinp2, cosp2, t_width, d_width, &
      & cut_distance, order, pd, ang_norm2, distance_scale, angular_scale, alchemy, verbose,&
      &self_scalar2)

      do ii = 1, n1(1)
         do jj = 1, n2(1)

            s12 = scalar(v1(1, ii, :, :), v2(1, jj, :, :), &
            & nneigh1_i(1, ii), nneigh2_j(1, jj), ksi1(1, ii, :), ksi2(1, jj, :), &
            & sinp1(1, ii, :, :, :), sinp2(1, jj, :, :, :), &
            & cosp1(1, ii, :, :, :), cosp2(1, jj, :, :, :), &
            & t_width, d_width, cut_distance, order, &
            & pd, ang_norm2, distance_scale, angular_scale, alchemy)

            !kernels(:, a, b) = kernels(:, a, b) &
            !    & + kernel(self_scalar1(a, i), self_scalar2(b, j), s12, &
            !    & kernel_idx, parameters)

            ktmp = 0.0d0
            call kernel(self_scalar1(1, ii), self_scalar2(1, jj), s12, &
            & kernel_idx, parameters, ktmp)
            kernels = kernels + ktmp

         end do
      end do

      kernel_element = sum(kernels)

      deallocate (v1)
      deallocate (v2)
      deallocate (kernels)
      deallocate (ktmp)
      deallocate (self_scalar1)
      deallocate (self_scalar2)
      deallocate (ksi1)
      deallocate (ksi2)
      deallocate (cosp1)
      deallocate (cosp2)
      deallocate (sinp1)
      deallocate (sinp2)

   end function kernel_element

   function train_kernel_element(i, j) result(k)
      integer, intent(in) :: i, j
      double precision :: k
      k = kernel_element(i, j, 'train')
   end function train_kernel_element

   function test_kernel_element(i, j) result(k)
      integer, intent(in) :: i, j
      double precision :: k
      k = kernel_element(i, j, 'test')
   end function test_kernel_element

   function test_train_kernel_element(i, j) result(k)
      integer, intent(in) :: i, j
      double precision :: k
      k = kernel_element(i, j, 'test_train')
   end function test_train_kernel_element

   function train_distance(i, j) result(d)
      integer, intent(in) :: i, j
      double precision :: d, d_squared

      d_squared = train_kernel_element(i, i) + train_kernel_element(j, j) - 2 * train_kernel_element(i, j)
      d_squared = max(d_squared, 0.0d0)
      d = sqrt(d_squared)

   end function train_distance

   function test_distance(i, j) result(d)
      integer, intent(in) :: i, j
      double precision :: d, d_squared

      d_squared = test_kernel_element(i, i) + train_kernel_element(j, j) - 2 * test_train_kernel_element(i, j)
      d_squared = max(d_squared, 0.0d0)
      d = sqrt(d_squared)

   end function test_distance

end module fchl18_kernel_distance
