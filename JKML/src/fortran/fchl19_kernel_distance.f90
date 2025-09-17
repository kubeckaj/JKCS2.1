module fchl19_kernel_distance
   implicit none

   private
   ! use target to declare as pointer-able (to avoid copying x1 and x2 a million times)
   !f2py intent(out) :: x1, x2
   double precision, dimension(:, :, :), allocatable :: X_tr, X_te
   integer, dimension(:, :), allocatable :: Q_tr, Q_te
   ! Number of molecules
   integer :: num_molecules1, num_molecules2
   ! List of numbers of atoms in each molecule/cluster
   integer, dimension(:), allocatable :: num_atoms1, num_atoms2
   integer :: maxatoms1, maxatoms2
   ! sigma of the gaussian kernel
   double precision :: sigma
   double precision :: inv_sigma2

   logical :: verbose
   public :: init_train, init_test
   public :: train_kernel_element, test_kernel_element, test_train_kernel_element, test_distance, train_distance
   private :: X_tr, X_te
   private :: kernel_element

contains

   subroutine clean_train()
      ! clean up training data if it exists to avoid double alloc
      if (allocated(X_tr)) deallocate(X_tr)
      if (allocated(Q_tr)) deallocate(Q_tr)
      if (allocated(num_atoms1)) deallocate(num_atoms1)
   end subroutine clean_train

   subroutine init_train(x_in1, q_tr_in, verbose_in, n1, nm1, sigma_in)
      ! fchl19 descriptors for the training set
      double precision, dimension(:, :, :), intent(in) :: x_in1 ! can be saved

      integer, dimension(:, :), intent(in) :: q_tr_in

      ! Whether to be verbose with output
      logical, intent(in) :: verbose_in

      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n1 ! save

      ! Number of molecules
      integer, intent(in) :: nm1 ! save

      ! sigma in
      double precision, intent(in) :: sigma_in

      ! clean up old data
      call clean_train()
      ! allocate raw datas
      allocate(X_tr(size(x_in1,1), size(x_in1,2), size(x_in1,3)))
      X_tr = x_in1

      ! allocate charges
      allocate(Q_tr(size(q_tr_in,1), size(q_tr_in,2)))
      Q_tr = q_tr_in

      ! save sigmas
      allocate(num_atoms1(size(n1,1)))
      num_atoms1 = n1

      maxatoms1 = maxval(n1)

      num_molecules1 = nm1

      sigma = sigma_in
      inv_sigma2 = -1.0d0 / (2 * sigma ** 2)

      verbose = verbose_in
   end subroutine init_train

   subroutine clean_test()
      ! deallocate everything related to X_test
      if (allocated(X_te)) then
         deallocate(X_te)
      end if
      if (allocated(Q_te)) deallocate(Q_te)
      if (allocated(num_atoms2)) then
         deallocate(num_atoms2)
      end if
   end subroutine clean_test

   subroutine init_test(x_in2, q_te_in, n2, nm2)
      ! fchl descriptors for the training set, format (i,maxatoms,5,maxneighbors)
      double precision, dimension(:, :, :), intent(in) :: x_in2 ! can be saved

      integer, dimension(:, :), intent(in) :: q_te_in

      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n2 ! save

      ! Number of molecules
      integer, intent(in) :: nm2 ! save

      ! clean up previously allocated values
      call clean_test()
      ! allocate raw datas
      allocate(X_te(size(x_in2,1), size(x_in2,2), size(x_in2,3)))
      X_te = x_in2

      ! allocate charges
      allocate(Q_te(size(q_te_in,1), size(q_te_in,2)))
      Q_te = q_te_in

      ! allocate atom number lists
      allocate(num_atoms2(size(n2,1)))
      num_atoms2 = n2

      maxatoms2 = maxval(n2)

      num_molecules2 = nm2
   end subroutine init_test

   double precision function kernel_element(i, j, mode)

      ! Internal counters
      character(len=*), intent(in) :: mode
      integer, intent(in) :: i, j
      integer :: ii, jj

      ! single elements from x1, x2
      double precision, dimension(:, :, :), allocatable :: v1
      double precision, dimension(:, :, :), allocatable :: v2
      ! charges for molecules i and j
      integer, dimension(:), allocatable :: q1
      integer, dimension(:), allocatable :: q2

      integer :: n1(1)
      integer :: n2(1)

      ! Work kernel space
      double precision :: l2

      ! copy relevant elements to the pointers
      select case (trim(mode))
       case ("train")
         ! train x train kernel
         allocate(v1(1,size(X_tr,2),size(X_tr,3)))
         allocate(v2(1,size(X_tr,2),size(X_tr,3)))
         v1(1,:,:) = X_tr(i, :, :)
         v2(1,:,:) = X_tr(j, :, :)
         n1 = [num_atoms1(i)]
         n2 = [num_atoms1(j)]
         allocate(q1(num_atoms1(i)))
         allocate(q2(num_atoms1(j)))
         q1 = Q_tr(:, i)
         q2 = Q_tr(:, j)
       case ("test")
         ! test x test kernel
         allocate(v1(1,size(X_te,2),size(X_te,3)))
         allocate(v2(1,size(X_te,2),size(X_te,3)))
         v1(1,:,:) = X_te(i, :, :)
         v2(1,:,:) = X_te(j, :, :)
         n1 = [num_atoms2(i)]
         n2 = [num_atoms2(j)]
         allocate(q1(num_atoms2(i)))
         allocate(q2(num_atoms2(j)))
         q1 = Q_te(:, i)
         q2 = Q_te(:, j)
       case ("test_train")
         ! test x train kernel
         allocate(v1(1,size(X_te,2),size(X_te,3)))
         allocate(v2(1,size(X_tr,2),size(X_tr,3)))
         v1(1,:,:) = X_te(i, :, :)
         v2(1,:,:) = X_tr(j, :, :)
         ! admittedly confusing, but now 1 = 2
         n1 = [num_atoms2(i)]
         n2 = [num_atoms1(j)]
         allocate(q1(num_atoms2(i)))
         allocate(q2(num_atoms1(j)))
         q1 = Q_te(:, i)
         q2 = Q_tr(:, j)
       case default
         print *, "Unknown kernel mode: ", mode
         stop
      end select

      kernel_element = 0.0d0

      do ii = 1, n1(1)
         do jj = 1, n2(1)
            if (q1(ii) == q2(jj)) then
               l2 = sum((v1(1, ii, :) - v2(1, jj, :)) ** 2)
               kernel_element = kernel_element + exp(l2 * inv_sigma2)
            endif
         end do
      end do


      deallocate (v1)
      deallocate (v2)
      deallocate (q1)
      deallocate (q2)

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

end module fchl19_kernel_distance
