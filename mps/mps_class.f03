module mps_class_module
    implicit none
    private

    type, public :: mps_class
        private
        integer :: field_name
    contains
        procedure :: method_name
        final :: destructor
    end type mps_class

contains

    subroutine method_name(self)
        class(mps_class), intent(in) :: self
    end subroutine
    
    subroutine destructor(self)
        type(mps_class), intent(in) :: self
    end subroutine

end module mps_class_module
