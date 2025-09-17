module fchl19_vp_tree
   use fchl19_kernel_distance
   use omp_lib
   implicit none

   ! flag for having test data
   logical :: has_test = .false.
   ! training indices, left and right children arrays
   integer, allocatable :: vp_index(:)
   integer, allocatable :: vp_left(:)
   integer, allocatable :: vp_right(:)
   ! threshold array
   double precision, allocatable :: vp_threshold(:)
   integer :: vp_next = 1

   ! train, test data
   double precision, allocatable :: X_train(:,:,:), Y_train(:), X_test(:,:,:)
contains

   subroutine init_tree(max_nodes)
      implicit none
      integer, intent(in) :: max_nodes

      ! clean up old
      if (allocated(vp_index)) deallocate(vp_index)
      if (allocated(vp_left)) deallocate(vp_left)
      if (allocated(vp_right)) deallocate(vp_right)
      if (allocated(vp_threshold)) deallocate(vp_threshold)

      allocate(vp_index(max_nodes))
      allocate(vp_left(max_nodes))
      allocate(vp_right(max_nodes))
      allocate(vp_threshold(max_nodes))

      vp_next = 1
   end subroutine init_tree

   subroutine clean_training_data()
      if (allocated(X_train)) deallocate(X_train)
      if (allocated(Y_train)) deallocate(Y_train)
      if (allocated(vp_index)) deallocate(vp_index)
      if (allocated(vp_left)) deallocate(vp_left)
      if (allocated(vp_right)) deallocate(vp_right)
      if (allocated(vp_threshold)) deallocate(vp_threshold)

   end subroutine clean_training_data

   subroutine load(X_in, Q_in, Y_in, vp_index_in, vp_left_in, vp_right_in, vp_threshold_in, verbose_in, n1, nm1, sigma_in)
      ! raw input data
      double precision, intent(in) :: X_in(:, :, :), Y_in(:)
      integer, intent(in) :: Q_in(:, :)

      ! previously learned VP params
      integer, intent(in) :: vp_index_in(:), vp_left_in(:), vp_right_in(:)
      double precision, intent(in) :: vp_threshold_in(:)

      ! Whether to be verbose with output
      logical, intent(in) :: verbose_in

      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n1 ! save

      ! Number of molecules
      integer, intent(in) :: nm1 ! save

      ! Number of sigmas
      double precision, intent(in) :: sigma_in ! save

      integer :: n_train

      ! clean up training data (if already allocated)
      call clean_training_data()
      allocate(X_train(size(X_in,1), size(X_in,2), size(X_in,3)))
      X_train = X_in
      allocate(Y_train(size(Y_in,1)))
      Y_train = Y_in
      ! get size of X_train
      n_train = size(X_train, 1)

      ! initialize vp-tree
      call init_tree(n_train)
      ! copy learned vp tree parameters to collections
      vp_index = vp_index_in
      vp_left = vp_left_in
      vp_right = vp_right_in
      vp_threshold = vp_threshold_in

      ! initialise kernel
      call init_train(X_train, Q_in, verbose_in, n1, nm1, sigma_in)

   end subroutine load

   subroutine train(X_in, Q_in, Y_in, verbose_in, n1, nm1, sigma_in, root_id)

      double precision, intent(in) :: X_in(:, :, :), Y_in(:)
      integer, intent(in) :: Q_in(:, :)

      ! Whether to be verbose with output
      logical, intent(in) :: verbose_in

      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n1 ! save

      ! Number of molecules
      integer, intent(in) :: nm1 ! save

      ! Number of sigmas
      double precision, intent(in) :: sigma_in ! save

      integer, intent(out) :: root_id

      integer :: n_train, i
      integer, allocatable :: indices(:)

      ! openmp parameters
      integer :: nthreads, max_depth

      ! clean up training data (if already allocated)
      call clean_training_data()
      allocate(X_train(size(X_in,1), size(X_in,2), size(X_in,3)))
      X_train = X_in
      allocate(Y_train(size(Y_in,1)))
      Y_train = Y_in
      ! get size of X_train
      n_train = size(X_train, 1)

      ! initialize vp-tree
      call init_tree(n_train)

      ! initialise kernel
      call init_train(X_train, Q_in, verbose_in, n1, nm1, sigma_in)

      allocate(indices(n_train))
      do i = 1, n_train
         indices(i) = i
      end do

      ! set up multithreading parameters
      nthreads = omp_get_num_threads()
      ! max_depth = log2(nthreads) + 1 as the tree always splits into at most 2 branches
      max_depth = int(log(real(nthreads))/log(2.0d0)) + 1
      call build_vptree(indices, n_train, 0, max_depth, root_id)

      deallocate(indices)
   end subroutine train

   recursive subroutine build_vptree(indices, n, depth, max_parallel_depth, node_id)
      implicit none
      integer, intent(in) :: indices(n)
      integer, intent(in) :: n
      integer, intent(in) :: depth, max_parallel_depth
      integer, intent(out) :: node_id

      integer :: i, left_id, right_id
      double precision, allocatable :: dists(:)
      integer, allocatable :: left_set(:), right_set(:)
      integer :: count_left, count_right
      double precision :: median

      if (n == 0) then
         node_id = -1
         return
      end if

      ! Ensure only one thread updates vp_next at a time
      !$omp atomic capture
      node_id = vp_next
      vp_next = vp_next + 1
      !$omp end atomic

      vp_index(node_id) = indices(n)  ! choose last as vantage point

      if (n == 1) then
         vp_threshold(node_id) = 0.0
         vp_left(node_id) = -1
         vp_right(node_id) = -1
         return
      end if

      allocate(dists(n-1))
      do i = 1, n-1
         dists(i) = train_distance(indices(i), indices(n))  ! dist to vantage point
      end do

      median = quick_median(dists)

      allocate(left_set(n-1), right_set(n-1))
      count_left = 0
      count_right = 0
      do i = 1, n-1
         if (dists(i) <= median) then
            count_left = count_left + 1
            left_set(count_left) = indices(i)
         else
            count_right = count_right + 1
            right_set(count_right) = indices(i)
         end if
      end do

      vp_threshold(node_id) = median

      if (depth < max_parallel_depth) then
         !$omp parallel sections default(shared)
         !$omp section
         call build_vptree(left_set(1:count_left), count_left, depth+1, max_parallel_depth, left_id)
         !$omp section
         call build_vptree(right_set(1:count_right), count_right, depth+1, max_parallel_depth, right_id)
         !$omp end parallel sections
      else
         call build_vptree(left_set(1:count_left), count_left, depth+1, max_parallel_depth, left_id)
         call build_vptree(right_set(1:count_right), count_right, depth+1, max_parallel_depth, right_id)
      end if

      vp_left(node_id) = left_id
      vp_right(node_id) = right_id

      deallocate(dists, left_set, right_set)
   end subroutine build_vptree


   subroutine search_vp_tree(i, k, neighbors, distances, return_distances)
      integer, intent(in) :: i, k
      integer, intent(out) :: neighbors(k)
      double precision, optional, intent(out) :: distances(k)
      logical, optional, intent(in) :: return_distances

      logical :: ret_dists
      integer, allocatable :: heap_idx(:)
      double precision, allocatable :: heap_dist(:)
      integer :: heap_size

      ret_dists = .false.
      if (present(return_distances)) ret_dists = return_distances
      allocate(heap_idx(k))
      allocate(heap_dist(k))
      heap_size = 0

      call search_node(1, i, k, heap_idx, heap_dist, heap_size)

      neighbors = heap_idx

      if (ret_dists) distances = heap_dist

      deallocate(heap_idx, heap_dist)
   end subroutine search_vp_tree

   subroutine set_up_test(X_test_in, Q_test_in, n2, nm2)
      ! List of numbers of atoms in each molecule
      integer, dimension(:), intent(in) :: n2 ! save

      ! Number of molecules
      integer, intent(in) :: nm2 ! save


      double precision, intent(in) :: X_test_in(:,:,:)
      integer, intent(in) :: Q_test_in(:, :)

      ! free previous X_test
      if (allocated(X_test)) then
         deallocate(X_test)
      end if
      allocate(X_test(size(X_test_in,1),size(X_test_in,2),size(X_test_in,3)))
      X_test = X_test_in
      has_test = .true.
      call init_test(X_test, Q_test_in, n2, nm2)

   end subroutine set_up_test

   subroutine predict(k, n_test, Y_pred, weight_by_distance)
      integer, intent(in) :: k, n_test
      double precision, dimension(n_test), intent(out) :: Y_pred
      logical, optional, intent(in) :: weight_by_distance

      integer :: i, j
      integer, allocatable :: neighbors(:,:)
      double precision, allocatable :: distances(:,:)
      double precision :: sum, weight_sum, w
      logical :: use_weights

      if (.not. has_test) then
         print *, "Test data has not been provided! Call set_up_test first!"
         stop
      end if
      use_weights = .false.
      if (present(weight_by_distance)) use_weights = weight_by_distance

      allocate(neighbors(k,n_test))

      if (use_weights) then
         allocate(distances(k,n_test))
         call kneighbors(k, n_test, neighbors,.true.,distances)
      else
         call kneighbors(k, n_test, neighbors)
      end if
      do i = 1, n_test
         sum = 0.0d0
         if (use_weights) then
            weight_sum = 0.0d0
            do j = 1, k
               w = 1.0d0 / max(distances(j,i), 1.0d-12) ! avoid div-by-0
               sum = sum + w * Y_train(neighbors(j, i))
               weight_sum = weight_sum + w
            end do
            Y_pred(i) = sum / weight_sum
         else
            do j = 1, k
               sum = sum + Y_train(neighbors(j, i))
            end do
            Y_pred(i) = sum / k
         end if
      end do

      deallocate(neighbors)
      if (use_weights) deallocate(distances)

   end subroutine predict

   subroutine kneighbors(k, n_test, neighbor_indices, return_distances, neighbor_distances)
      integer, intent(in) :: k, n_test
      integer, dimension(k,n_test), intent(out) :: neighbor_indices
      logical, intent(in), optional :: return_distances
      double precision, intent(out), optional, dimension(k,n_test) :: neighbor_distances

      integer :: i

      logical :: do_return_distances

      if (.not. has_test) then
         print *, "Test data has not been provided! Call set_up_test first!"
         stop
      end if

      do_return_distances = .false.
      if (present(return_distances)) do_return_distances = return_distances

      !$omp parallel do default(shared) private(i)
      do i = 1, n_test
         if (do_return_distances) then
            call search_vp_tree(i, k, neighbor_indices(:, i), neighbor_distances(:, i), do_return_distances)
         else
            call search_vp_tree(i, k, neighbor_indices(:, i))
         end if
      end do
      !$omp end parallel do
   end subroutine kneighbors

   recursive subroutine search_node(node, test_idx, k, heap_idx, heap_dist, heap_size)
      integer, intent(in) :: node, test_idx, k
      integer, intent(inout) :: heap_idx(:), heap_size
      double precision, intent(inout) :: heap_dist(:)

      integer :: tr_idx, left, right
      double precision :: d, diff, threshold
      integer :: i, worst_idx
      double precision :: worst_dist

      if (node == -1) return

      tr_idx = vp_index(node)
      d = test_distance(test_idx, tr_idx)

      ! Insert into heap (brute force k-NN buffer)
      if (heap_size < k) then
         heap_size = heap_size + 1
         heap_idx(heap_size) = tr_idx
         heap_dist(heap_size) = d
      else
         worst_idx = 1
         worst_dist = heap_dist(1)
         do i = 2, k
            if (heap_dist(i) > worst_dist) then
               worst_idx = i
               worst_dist = heap_dist(i)
            end if
         end do
         if (d < worst_dist) then
            heap_dist(worst_idx) = d
            heap_idx(worst_idx) = tr_idx
         end if
      end if

      threshold = vp_threshold(node)
      diff = d - threshold
      left = vp_left(node)
      right = vp_right(node)

      if (diff < 0) then
         call search_node(left, test_idx, k, heap_idx, heap_dist, heap_size)
         if (abs(diff) < maxval(heap_dist(1:heap_size))) then
            call search_node(right, test_idx, k, heap_idx, heap_dist, heap_size)
         end if
      else
         call search_node(right, test_idx, k, heap_idx, heap_dist, heap_size)
         if (abs(diff) < maxval(heap_dist(1:heap_size))) then
            call search_node(left, test_idx, k, heap_idx, heap_dist, heap_size)
         end if
      end if
   end subroutine search_node

   recursive function quickselect(arr, left, right, k) result(pivot)
      double precision, intent(inout) :: arr(:)
      integer, intent(in) :: left, right, k
      double precision :: pivot
      integer :: pivot_index, i, j
      double precision :: temp, pivot_value

      if (left == right) then
         pivot = arr(left)
         return
      end if

      ! Deterministic pivot: middle element
      pivot_index = (left + right) / 2
      pivot_value = arr(pivot_index)

      ! Swap pivot to end
      temp = arr(pivot_index)
      arr(pivot_index) = arr(right)
      arr(right) = temp

      ! Partition
      j = left
      do i = left, right - 1
         if (arr(i) < pivot_value) then
            temp = arr(i)
            arr(i) = arr(j)
            arr(j) = temp
            j = j + 1
         end if
      end do

      ! Move pivot to final position
      temp = arr(j)
      arr(j) = arr(right)
      arr(right) = temp

      ! Recursive selection
      if (k == j) then
         pivot = arr(j)
      else if (k < j) then
         pivot = quickselect(arr, left, j - 1, k)
      else
         pivot = quickselect(arr, j + 1, right, k)
      end if
   end function quickselect

   function quick_median(arr) result(median)
      double precision, intent(in) :: arr(:)
      double precision :: median
      double precision, allocatable :: temp(:)
      integer :: n, mid

      n = size(arr)
      allocate(temp(n))
      temp = arr

      mid = (n + 1) / 2

      if (mod(n, 2) == 0) then
         median = 0.5 * (quickselect(temp, 1, n, mid) + quickselect(temp, 1, n, mid + 1))
      else
         median = quickselect(temp, 1, n, mid)
      end if

      deallocate(temp)
   end function quick_median

end module fchl19_vp_tree
