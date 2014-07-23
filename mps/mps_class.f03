module mps_class_module
    implicit none
    private

    type, public :: mps
        private
        type(mp_node),dimension(L) :: node_array
    contains
        function :: mps
        final :: delete
    end type mps

    type :: mp_node
        complex, allocatable, dimension(:,:) :: matrix
    end type mp_node

    interface mps
        module procedure zero_mps
    end interface mps

contains

    function zero_mps(d,D,L) result(state)

        integer, intent(in) :: d   !physical on-site dimension
        integer, intent(in) :: D   !maximal on-site virtual dimension
        integer, intent(in) :: L    ! Length of one-dimensional system
        class(mps), intent(out) :: state

        do i = 1, L
            allocate(state(i),)
        end do
    end function
    
    subroutine delete(self)
        type(mps), intent(in) :: self
    end subroutine

end module mps_class_module
