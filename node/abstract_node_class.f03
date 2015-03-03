module abstract_node_class
    use IR_Precision
    implicit none
    private


    type, public, abstract :: node
    private
        integer, allocatable, dimension(:) :: node_outlets
    contains
        ! RANDOM INITIALIZER
        procedure(abs_randInit), deferred :: randInit
        ! DENSE GETTERS
        procedure, private :: one_outlet
        procedure, private :: all_outlets
        generic :: outlets => one_outlet, all_outlets
        procedure(abs_elements), deferred  :: elements
        ! DENSE SETTERS
        procedure, non_overridable :: setElements
        procedure(abs_safe_setElements), deferred :: safe_setElements
        procedure, non_overridable :: setNode
        ! METHODS
        procedure(abs_safe_trace), deferred :: safe_trace
        procedure, non_overridable :: trace
        ! procedure(abs_svd), deferred :: svd
    end type node

    abstract interface
        subroutine abs_randInit(self,nde_outlets,sparsity)
            import :: node, R_P
            class(node), intent(inout) :: self
            integer, dimension(:), intent(in) :: nde_outlets
            real(R_P), intent(in) :: sparsity
        end subroutine abs_randInit
    end interface

    abstract interface
        pure function abs_elements(self) result(nde_elements)
            import :: node, R_P
            class(node), intent(in) :: self
            complex(R_P), allocatable, dimension(:) :: nde_elements
        end function abs_elements
    end interface

    abstract interface
        subroutine abs_safe_setElements(self,nde_elements)
            import :: node, R_P
            class(node), intent(inout) :: self
            complex(R_P), dimension(:), intent(in) :: nde_elements
        end subroutine abs_safe_setElements
    end interface

    abstract interface
        subroutine abs_safe_trace(nodeRes,node1,tracePairs,node2)
            import :: node
            class(node), intent(out) :: nodeRes
            class(node), intent(in) :: node1
            integer, dimension(:,:), intent(in) :: tracePairs
            class(node), intent(in) :: node2
        end subroutine abs_safe_trace
    end interface


contains

    function one_outlet(self,ind) result(nde_outlet)

        class(node), intent(in) :: self
        integer, intent(in) :: ind
        integer :: nde_outlet

        if((ind < 1) .or.  (ind > size(self%node_outlets))) then
            print *, "Error in Outlet Getter:"
            print *, ""
            print *, "Outlet requested: ", ind, " is out of range: 1.. ", size(self%node_outlets)
            call exit()
        end if

        nde_outlet = self%node_outlets(ind)

    end function one_outlet


    pure function all_outlets(self) result(nde_outlets)

        class(node), intent(in) :: self
        integer, allocatable, dimension(:) :: nde_outlets

        allocate(nde_outlets(size(self%node_outlets)), source = self%node_outlets)

    end function all_outlets


    subroutine setElements(self,nde_elements)

        class(node), intent(inout) :: self
        complex(R_P), dimension(:), intent(in) :: nde_elements

        if( size(nde_elements) /= product(self%node_outlets)) then
            print *, "Error in Element Setter:"
            print *, ""
            print *, "Number of elements ", size(nde_elements), " not compatible with outlets: ", self%node_outlets
            call exit()
        end if

        call self%safe_setElements(nde_elements)

    end subroutine setElements


    subroutine setNode(self,nde_outlets,nde_elements)

        class(node), intent(inout) :: self
        integer, dimension(:), intent(in) :: nde_outlets
        complex(R_P), dimension(:), intent(in) :: nde_elements

        if(size(nde_elements) /= product(nde_outlets)) then
            print *, "Error in Node Setter:"
            print *, ""
            print *, "Number of elements ", size(nde_elements), " not compatible with outlets: ", nde_outlets
            call exit()
        end if

        self%node_outlets = nde_outlets
        call self%safe_setElements(nde_elements)

    end subroutine setNode


    subroutine trace(nodeRes,node1,tracePairs,node2)

        class(node), intent(out) :: nodeRes
        class(node), intent(in) :: node1
        integer, dimension(:,:), intent(in) :: tracePairs
        class(node), intent(in) :: node2


        integer, dimension(2) :: tracePair
        integer :: j

        do j=1,size(tracePairs,2)
            tracePair = tracePairs(:,j)

            if((tracePair(1) < 1) .or. (tracePair(1) > size(node1%node_outlets))) then
                print *, "Error in Trace Operation:"
                print *, ""
                print *, "contract dimension not compatible with rank of node"
                print *, "contract dimension: ", tracePair(1)
                print *, "rank of node: ", size(node1%node_outlets)
                call exit()
            end if

            if((tracePair(2) < 1) .or. (tracePair(2) > size(node2%node_outlets))) then
                print *, "Error in Trace Operation:"
                print *, ""
                print *, "contract dimension not compatible with rank of node"
                print *, "contract dimension: ", tracePair(2)
                print *, "rank of node: ", size(node2%node_outlets)
                call exit()
            end if

            if( node1%node_outlets(tracePair(1)) /= node2%node_outlets(tracePair(2))) then
                print *, "Error in Trace Operation:"
                print *, ""
                print *, "unable to contract given dimensions because of differing range"
                print *, "dimension: ", tracePair(1), " range: 1 ..", node1%node_outlets(tracePair(1))
                print *, "dimension: ", tracePair(2), " range: 1 ..", node2%node_outlets(tracePair(2))
                call exit()
            end if

        end do

        call nodeRes%safe_trace(node1,tracePairs,node2)

    end subroutine trace

end module abstract_node_class
