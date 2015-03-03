module node_class
    use IR_Precision
    implicit none
    private

    type, public :: node
    private
        complex(R_P), allocatable, dimension(:) :: elements
        integer, allocatable, dimension(:) :: node_shape
    contains
        procedure, private :: tuple_to_index
        procedure, private :: index_to_tuple
        generic :: tuple => index_to_tuple, tuple_to_index
        procedure, private :: tensor_product
        generic :: operator(*) => tensor_product
        procedure :: contract
        procedure :: get_elements
        procedure :: get_shape
        final :: delete
    end type node

    private :: module_tuple_to_index
    private :: module_index_to_tuple

    public :: merge
    public :: transpose

    interface module_tuple
        module procedure module_tuple_to_index
        module procedure module_index_to_tuple
    end interface module_tuple

    interface node
        module procedure new_node
        module procedure new_node_rank1
        module procedure new_node_rank2
        module procedure new_node_rank3
        module procedure new_node_rank4
    end interface node

    interface merge
        module procedure node_merge
    end interface merge

    interface transpose
        module procedure module_transpose
    end interface transpose

contains

! Type Constructors

    function new_node(elements,node_shape)

        complex(R_P), dimension(:), intent(in) :: elements
        integer, dimension(:), intent(in) :: node_shape
        type(node) :: new_node

        if( size(elements) /= product(node_shape)) then
            print *, "Non-compatible number of array elements: ", size(elements), " and node shape: ", node_shape
            call abort()
        end if

        allocate(new_node%elements(size(elements)), source=elements)
        allocate(new_node%node_shape(size(node_shape)), source=node_shape)

    end function new_node


    function new_node_rank1(rank1_array)

        complex(R_P), dimension(:), intent(in) :: rank1_array
        type(node) :: new_node_rank1

        new_node_rank1 = new_node(reshape(rank1_array,[size(rank1_array)]),shape(rank1_array))

    end function new_node_rank1


    function new_node_rank2(rank2_array)

        complex(R_P), dimension(:,:), intent(in) :: rank2_array
        type(node) :: new_node_rank2

        new_node_rank2 = new_node(reshape(rank2_array,[size(rank2_array)]),shape(rank2_array))

    end function new_node_rank2


    function new_node_rank3(rank3_array)

        complex(R_P), dimension(:,:,:), intent(in) :: rank3_array
        type(node) :: new_node_rank3

        new_node_rank3 = new_node(reshape(rank3_array,[size(rank3_array)]),shape(rank3_array))

    end function new_node_rank3


    function new_node_rank4(rank4_array)

        complex(R_P), dimension(:,:,:,:), intent(in) :: rank4_array
        type(node) :: new_node_rank4

        new_node_rank4 = new_node(reshape(rank4_array,[size(rank4_array)]),shape(rank4_array))

    end function new_node_rank4


    function index_to_tuple(self,ind) result(tup)

        class(node), intent(in) :: self
        integer, intent(in) :: ind
        integer, dimension(size(self%node_shape)) :: tup

        tup = module_index_to_tuple(self%node_shape,ind)

    end function index_to_tuple


    function tuple_to_index(self,tup) result(ind)

        class(node), intent(in) :: self
        integer, dimension(:) :: tup
        integer :: ind

        ind = module_tuple_to_index(self%node_shape,tup)

    end function tuple_to_index


    function tensor_product(node1,node2) result(res_node)

        class(node), intent(in) :: node1
        class(node), intent(in) :: node2
        type(node) :: res_node

        integer, dimension(size(node1%node_shape) + size(node2%node_shape)) :: res_node_shape
        complex(R_P), dimension(product(node1%node_shape)*product(node2%node_shape)) :: res_node_elements
        integer :: i, size_node1, size_node2

        size_node1 = size(node1%elements)
        size_node2 = size(node2%elements)

        res_node_shape(:size(node1%node_shape)) = node1%node_shape
        res_node_shape((size(node1%node_shape) + 1):) = node2%node_shape

        do i = 1,size_node2
            res_node_elements((((i-1) * size_node1) + 1):(i * size_node1) ) = node1%elements * node2%elements(i)
        end do

        res_node = new_node(res_node_elements,res_node_shape)

    end function


    subroutine contract(self,contract_pair)

        class(node), intent(inout) :: self
        integer, dimension(2), intent(in) :: contract_pair

        integer, dimension(size(self%node_shape) - 2) :: new_node_shape
        complex(R_P), dimension(product(self%node_shape)/self%node_shape(contract_pair(1))**2) :: new_elements
        logical, dimension(size(self%node_shape)) :: dim_mask
        integer, dimension(size(self%node_shape)) :: old_index
        integer, dimension(size(self%node_shape) - 2) :: new_index
        integer :: i

        do i=1,2
            if((contract_pair(i) < 1) .or. (contract_pair(i) > size(self%node_shape))) then
                print *, "contract dimension not compatible with rank of node"
                print *, "contract dimension: ", contract_pair(i)
                print *, "rank of node: ", size(self%node_shape)
                call exit()
            end if
        end do

        if( self%node_shape(contract_pair(1)) /= self%node_shape(contract_pair(2))) then
            print *, "unable to contract given dimensions because of differing range"
            print *, "dimension: ", contract_pair(1), " range: 1 ..", self%node_shape(contract_pair(1))
            print *, "dimension: ", contract_pair(2), " range: 1 ..", self%node_shape(contract_pair(2))
            call exit()
        end if

        dim_mask = .true.
        dim_mask(contract_pair(1)) = .false.
        dim_mask(contract_pair(2)) = .false.
        new_elements = (0.0, 0.0_R_P)
        new_node_shape = pack(self%node_shape,dim_mask)

        do i = 1, size(self%elements)
            old_index = self%tuple(i)
            if(old_index(contract_pair(1))==old_index(contract_pair(2))) then
                new_index = pack(self%tuple(i),dim_mask)
                new_elements(module_tuple(new_node_shape,new_index)) = new_elements(module_tuple(new_node_shape,new_index)) &
                         + self%elements(i)
            end if
        end do

        deallocate(self%elements)
        deallocate(self%node_shape)
        allocate(self%elements(size(new_elements)),source=new_elements)
        allocate(self%node_shape(size(new_node_shape)),source=new_node_shape)

    end subroutine


    function node_merge(node1,node2,contract_pairs) result(res_node)

        type(node), intent(in) :: node1
        type(node), intent(in) :: node2
        integer, dimension(:,:), intent(in) :: contract_pairs
        type(node) :: res_node

        integer :: i
        integer, dimension(2) :: contract_pair
        integer, dimension(size(contract_pairs,1),size(contract_pairs,2)) :: local_contract_pairs

        local_contract_pairs(1,:) = contract_pairs(1,:)
        local_contract_pairs(2,:) = contract_pairs(2,:) + size(node1%get_shape())

        res_node = node1 * node2

        do i = 1,size(contract_pairs,2)
            contract_pair = local_contract_pairs(:,i)
            call contract(res_node,contract_pair)
            where(local_contract_pairs > contract_pair(1)) local_contract_pairs = local_contract_pairs - 1
            where(local_contract_pairs > contract_pair(2)) local_contract_pairs = local_contract_pairs - 1
        end do

    end function node_merge


    function module_transpose(node1,transpose_pairs) result(res_node)

        type(node), intent(in) :: node1
        integer, dimension(:,:) :: transpose_pairs
        type(node) :: res_node

        complex(R_P), dimension(size(node1%elements)) :: res_elements
        integer(R_P), dimension(size(node1%node_shape)) :: res_node_shape
        integer :: i, new_i

        res_elements = (0.0,0.0_R_P)
        res_node_shape = exchange(node1%node_shape,transpose_pairs)
        res_node = node(res_elements,res_node_shape)

        do i = 1,size(node1%elements)
            new_i = res_node%tuple(exchange(node1%tuple(i),transpose_pairs))
            res_node%elements(new_i) = node1%elements(i)
        end do

    contains

        function exchange(index_vec,transpose_pairs) result(trans_index_vec)

            integer, dimension(:), intent(in) :: index_vec
            integer, dimension(:,:), intent(in) :: transpose_pairs
            integer, dimension(size(index_vec)) :: trans_index_vec

            integer :: i

            trans_index_vec = index_vec

            do i = 1,size(transpose_pairs,2)
                trans_index_vec(transpose_pairs(2,i)) = index_vec(transpose_pairs(1,i))
                trans_index_vec(transpose_pairs(1,i)) = index_vec(transpose_pairs(2,i))
            end do

        end function exchange

    end function module_transpose


! Type Accessors

    function get_elements(self)

        class(node) :: self
        complex(R_P), allocatable, dimension(:) :: get_elements

        allocate(get_elements(size(self%elements)),source=self%elements)

    end function get_elements


    function get_shape(self)

        class(node) ::  self
        integer, allocatable, dimension(:) :: get_shape

        allocate(get_shape(size(self%node_shape)),source=self%node_shape)

    end function

! Destructor

    elemental subroutine delete(self)

        type(node), intent(inout) :: self

        deallocate(self%elements)
        deallocate(self%node_shape)

    end subroutine delete

! Helper Functions

    function module_tuple_to_index(node_shape,tup) result(ind)

        integer, dimension(:) :: node_shape
        integer, dimension(:) :: tup
        integer :: ind

        integer, dimension(size(node_shape)) :: place_value
        integer :: i

        if(size(tup) /= size(node_shape)) then
            print *, "Tuple rank does not match rank of node"
            print *, "tuple rank: ", size(tup)
            print *, "node rank: ", size(node_shape)
            call exit()
        end if

        do i = 1, size(tup)
            if((tup(i) < 1) .or. (tup(i) > node_shape(i))) then
                print *, "Tuple index at dimension ", i," is out of valid range"
                print *, "tuple index: ", tup(i)
                print *, "valid range: 1 to ", node_shape(i)
                call exit()
            end if
        end do

        place_value = 1
        do i = 2, size(place_value)
            place_value(i) = product(node_shape(1:(i - 1)))
        end do

        ind = 1
        do i = 1, size(tup)
            ind = ind + (tup(i) - 1) * place_value(i)
        end do

    end function module_tuple_to_index


    function module_index_to_tuple(node_shape,ind) result(tup)

        integer, dimension(:), intent(in) :: node_shape
        integer, intent(in) :: ind
        integer, dimension(size(node_shape)) :: tup

        integer :: i, local_ind
        integer, dimension(size(node_shape)) :: place_value

        if(ind > product(node_shape)) then
            print *, "index is larger than number of elements"
            print *, "index: ", ind
            print *, "number of elements: ", product(node_shape)
        end if

        place_value = 1
        do i = 2, size(place_value)
            place_value(i) = product(node_shape(1:(i - 1)))
        end do

        local_ind = ind - 1
        do i = size(node_shape),2,-1
            tup(i) = 1 + local_ind / place_value(i)
            local_ind = modulo(local_ind,place_value(i))
        end do
        tup(1) = local_ind + 1

    end function module_index_to_tuple




end module node_class
