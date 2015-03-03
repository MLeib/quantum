module mpo_node_module
    use IR_Precision
    use mps_parameters
    use lapack_wrap
    implicit none
    private

    type, public :: mpo_node
        private
        complex(R_P), allocatable, dimension(:,:,:,:) :: state
        integer :: vDim1, vDim2
    contains
        procedure, private :: all_node
        procedure, private :: single_node
        generic :: node => all_node, single_node
        final :: delete
    end type mpo_node

    public :: size

    interface mpo_node
        module procedure new_mpo_node
    end interface mpo_node

    interface size
        module procedure mpo_node_size
        module procedure mpo_node_size_ind
    end interface size


contains

    function new_mpo_node(state) result(self)

        complex(R_P), dimension(:,:,:,:) :: state
        type(mpo_node) :: self

        self%vDim1 = size(state,1)
        self%vDim2 = size(state,2)

        if(allocated(self%state)) deallocate(self%state)
        allocate(self%state(self%vDim1,self%vDim2,pDim,pDim), source=state)

    end function


    pure function all_node(self) result(state)

        class(mpo_node), intent(in) :: self
        complex(R_P), allocatable, dimension(:,:,:,:) :: state

        allocate(state(self%vDim1,self%vDim2,pDim,pDim),source=self%state)

    end function all_node


    pure function single_node(self, pInd1, pInd2) result(state)

        class(mpo_node), intent(in) :: self
        integer, intent(in) :: pInd1
        integer, intent(in) :: pInd2
        complex(R_P), allocatable, dimension(:,:) :: state

        allocate(state(self%vDim1,self%vDim2),source=self%state(:,:,pInd1,pInd2))

    end function single_node


    pure function mpo_node_size(self) result(sze)

        type(mpo_node), intent(in) :: self
        integer :: sze

        sze = size(self%state)

    end function mpo_node_size


    pure function mpo_node_size_ind(self,dim_ind) result(sze)

        type(mpo_node), intent(in) :: self
        integer, intent(in) :: dim_ind
        integer :: sze

        sze = size(self%state,dim_ind)

    end function mpo_node_size_ind


    elemental subroutine delete(self)

        type(mpo_node), intent(inout) :: self

        deallocate(self%state)

    end subroutine delete



end module mpo_node_module

