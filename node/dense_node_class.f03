module dense_node_class
    use IR_Precision
    use abstract_node_class
    use random
    implicit none
    private

    type, public, extends(node) :: dense_node
    private
        complex(R_P), allocatable, dimension(:) :: node_elements
    contains
        procedure :: randInit
        procedure :: elements
        procedure :: safe_setElements
        procedure :: safe_trace
        procedure, private :: tuple_to_index
        procedure, private :: index_to_tuple
    end type dense_node

    interface dense_node
        module procedure new_dense_node
    end interface dense_node

contains

! TYPE CONSTRUCTOR

    function new_dense_node(outlets,elements) result(self)

        integer, dimension(:) :: outlets
        complex(R_P), dimension(:) :: elements
        type(dense_node) :: self

       call self%setNode(outlets,elements)

    end function new_dense_node

! RANDOM INITIALIZER

    subroutine randInit(self,nde_outlets,sparsity)

        class(dense_node), intent(inout) :: self
        integer, dimension(:), intent(in) :: nde_outlets
        real(R_P), intent(in) :: sparsity

        complex(R_P), dimension(product(nde_outlets)) :: nde_elements

        call normal(nde_elements,sparsity)

        call self%setNode(nde_outlets,nde_elements)

    end subroutine randInit


! ACCESSORS


    pure function elements(self) result(nde_elements)

        class(dense_node), intent(in) :: self
        complex(R_P), allocatable, dimension(:) :: nde_elements

        allocate(nde_elements(size(self%node_elements)), source = self%node_elements)

    end function elements


    subroutine safe_setElements(self,nde_elements)

        class(dense_node), intent(inout) :: self
        complex(R_P), dimension(:), intent(in) :: nde_elements

        self%node_elements = nde_elements

    end subroutine safe_setElements


! TRACE METHOD

    subroutine safe_trace(nodeRes,node1,tracePairs,node2)

        class(dense_node), intent(out) :: nodeRes
        class(node), intent(in) :: node1
        integer, dimension(:,:), intent(in) :: tracePairs
        class(node), intent(in) :: node2



        integer :: i,j, ind
        integer, dimension(size(node1%outlets())) :: tuple1
        integer, dimension(size(node2%outlets())) :: tuple2
        integer, allocatable, dimension(:) :: noderes_outlets
        complex(R_P), allocatable, dimension(:) :: noderes_elements
        logical, dimension(size(node1%outlets()) + size(node2%outlets())) :: merge_mask

        merge_mask = .true.
        merge_mask(tracePairs(1,:)) = .false.
        merge_mask(tracePairs(2,:) + size(node1%outlets())) = .false.
        tuple1 = 1
        tuple2 = 1
        noderes_outlets = pack([node1%outlets(),node2%outlets()], merge_mask)
        allocate(noderes_elements(product(noderes_outlets)))
        noderes_elements = (0.0,0.0_R_P)
        !allocate(nodeRes, source = dense_node(noderes_outlets,noderes_elements))
        call nodeRes%setNode(noderes_outlets,noderes_elements)

        select type(node1)
            type is (dense_node)
        select type(node2)
            type is (dense_node)


        do i = 1,size(node1%node_elements)

            do j = 1,size(node2%node_elements)

                if(all(tuple1(tracePairs(1,:)) == tuple2(tracePairs(2,:)))) then
                    ind = nodeRes%tuple_to_index(pack([tuple1,tuple2],merge_mask))
                    nodeRes%node_elements(ind) = nodeRes%node_elements(ind) + &
                        node1%node_elements(i) * node2%node_elements(j)
                end if
                call increment_node2(tuple2)

            end do
            call increment_node1(tuple1)

        end do

        end select
        end select

    contains

        subroutine increment_node1(tuple)

            integer, dimension(:), intent(inout) :: tuple

            integer :: i

            do i = 1, size(node1%outlets())
                if (tuple(i) < node1%outlets(i)) then
                    tuple(i) = tuple(i) + 1
                    exit
                else
                    tuple(i) = 1
                end if
            end do

        end subroutine increment_node1


        subroutine increment_node2(tuple)

            integer, dimension(:), intent(inout) :: tuple

            integer :: i

            do i = 1, size(node2%outlets())
                if (tuple(i) < node2%outlets(i)) then
                    tuple(i) = tuple(i) + 1
                    exit
                else
                    tuple(i) = 1
                end if
            end do

        end subroutine increment_node2

    end subroutine safe_trace


! SPECIFIC DENSE_NODE  METHODS

    function tuple_to_index(self,tup) result(ind)

        class(dense_node), intent(in) :: self
        integer, dimension(:) :: tup
        integer :: ind

        integer, dimension(size(self%outlets())) :: outlets
        integer, dimension(size(self%outlets())) :: place_value
        integer :: i

        outlets =  self%outlets()

        if(size(tup) /= size(outlets)) then
            print *, "Tuple rank does not match rank of node"
            print *, "tuple rank: ", size(tup)
            print *, "node rank: ", size(outlets)
            call exit()
        end if

        do i = 1, size(tup)
            if((tup(i) < 1) .or. (tup(i) > outlets(i))) then
                print *, "Tuple index at dimension ", i," is out of valid range"
                print *, "tuple index: ", tup(i)
                print *, "valid range: 1 to ", outlets(i)
                call exit()
            end if
        end do

        place_value = 1
        do i = 2, size(place_value)
            place_value(i) = product(outlets(1:(i - 1)))
        end do

        ind = 1
        do i = 1, size(tup)
            ind = ind + (tup(i) - 1) * place_value(i)
        end do

    end function tuple_to_index


    function index_to_tuple(self,ind) result(tup)

        class(dense_node), intent(in) :: self
        integer, intent(in) :: ind
        integer, dimension(size(self%outlets())) :: tup

        integer :: i, local_ind
        integer, dimension(size(self%outlets())) :: place_value
        integer, dimension(size(self%outlets())) :: outlets

        if(ind > product(outlets)) then
            print *, "index is larger than number of elements"
            print *, "index: ", ind
            print *, "number of elements: ", product(outlets)
        end if

        place_value = 1
        do i = 2, size(place_value)
            place_value(i) = product(outlets(1:(i - 1)))
        end do

        local_ind = ind - 1
        do i = size(outlets),2,-1
            tup(i) = 1 + local_ind / place_value(i)
            local_ind = modulo(local_ind,place_value(i))
        end do
        tup(1) = local_ind + 1

    end function index_to_tuple

end module dense_node_class
