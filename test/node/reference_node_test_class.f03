module reference_node_test_class
    use IR_Precision
    use abstract_test_class
    use random
    use lapack_wrap
    implicit none

    type, public, extends(test_case) :: reference_node_test
    private
        complex(R_P), allocatable, dimension(:,:) :: mat1, mat2, matRes
    contains
        procedure :: init
        procedure :: run_test
        procedure :: evaluate
    end type reference_node_test

contains

    subroutine init(self,test_size)

        class(reference_node_test), intent(inout) :: self
        integer, intent(in) :: test_size

        if(allocated(self%mat1)) then
            deallocate(self%mat1)
        end if
        if(allocated(self%mat2)) then
            deallocate(self%mat2)
        end if
        if(allocated(self%matRes)) then
            deallocate(self%matRes)
        end if
        allocate(self%mat1(test_size,test_size),self%mat2(test_size,test_size),self%matRes(test_size,test_size))
        call normal(self%mat1,1.0_R_P)
        call normal(self%mat2,1.0_R_P)

    end subroutine init


    subroutine run_test(self)

        class(reference_node_test), intent(inout) :: self

        self%matRes = multiply(self%mat1,self%mat2)

    end subroutine run_test


    function evaluate(self) result(test_res)

        class(reference_node_test), intent(in) :: self
        real(R_P) :: test_res

        test_res = 0.0_R_P

    end function evaluate

end module reference_node_test_class
