module dense_node_test_class
    use IR_Precision
    use abstract_test_class
    use dense_node_class
    implicit none

    type, public, extends(test_case) :: dense_node_test
    private
        type(dense_node) :: node1, node2, nodeRes
        complex(R_P), allocatable, dimension(:,:) :: mat1, mat2, matRes
    contains
        procedure :: init
        procedure :: run_test
        procedure :: evaluate
    end type dense_node_test

contains

    subroutine init(self,test_size)

        class(dense_node_test), intent(inout) :: self
        integer, intent(in) :: test_size

        call self%node1%randInit([test_size , test_size], 1.0_R_P)
        call self%node2%randInit([test_size , test_size], 1.0_R_P)
        self%mat1 = reshape(self%node1%elements(),[test_size, test_size])
        self%mat2 = reshape(self%node2%elements(),[test_size, test_size])
        self%matres = matmul(self%mat1,self%mat2)

    end subroutine init


    subroutine run_test(self)

        class(dense_node_test), intent(inout) :: self

        call self%nodeRes%trace(self%node1,reshape([2,1],[2,1]),self%node2)

    end subroutine run_test


    function evaluate(self) result(test_res)

        class(dense_node_test), intent(in) :: self
        real(R_P) :: test_res

        integer, dimension(2) :: outlets

        outlets = self%nodeRes%outlets()
        test_res = sqrt(sum(abs(self%matRes - reshape(self%nodeRes%elements(),outlets))**2))

    end function evaluate


end module dense_node_test_class
