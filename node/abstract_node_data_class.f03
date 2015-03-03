!module abstract_node_data_class
!    implicit none
!    private
!
!    type, public, abstract :: node_data
!    private
!    contains
!        procedure(abs_ToDense), deferred :: ToDense
!        procedure(abs_FromDense), deferred :: FromDense
!        procedure(abs_consistent), deferred :: int_consistent
!        generic, public :: operator(.consistent.) => int_consistent
!    end type node_data
!
!    abstract interface
!        subroutine toDense(self,nde_outlets,nde_elements)
!            import node_data, R_P
!            class(node_data), intent(in) :: self
!            integer, dimension(:), intent(in) :: nde_outlets
!            complex(R_P), allocatable, dimension(:), intent(out) :: nde_elements
!        end subroutine toDense
!    end interface
!
!    abstract interface
!        subroutine FromDense(self,nde_outlets,nde_elements)
!            import node_data, R_P
!            class(node_data), intent(out) :: self
!            integer, dimension(:), intent(in) :: nde_outlets
!            complex(R_P), dimension(:), intent(in) :: nde_elements
!        end subroutine FromDense
!    end interface
!
!    abstract interface
!        pure function abs_consistent(self,outlets) result(res)
!            import :: node_data
!            class(node_data), intent(in) :: self
!            integer, dimension(:), intent(in) :: outlets
!            logical :: res
!        end function abs_consistent
!    end interface
!
!end module abstract_node_data_class
