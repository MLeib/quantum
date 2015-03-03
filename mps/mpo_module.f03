module mpo_module
    use IR_Precision
    use mps_parameters
    use mpo_node_module
    implicit none
    private

    type, public :: mpo
        private
        type(mpo_node), allocatable, dimension(:) :: nodes
    contains
        procedure, private :: all_pDim_node
        procedure, private :: single_pDim_node
        generic :: node => all_pDim_node, single_pDim_node
    end type mpo

    public :: size

    interface mpo
        module procedure dedicated_new_mpo
    end interface mpo

    interface size
        module procedure mpo_size
    end interface size

contains

    function dedicated_new_mpo(start_matrix, middle_matrix, end_matrix) result(self)

        complex(R_P), dimension(:,:,:,:), intent(in) :: start_matrix
        complex(R_P), dimension(:,:,:,:), intent(in) :: middle_matrix
        complex(R_P), dimension(:,:,:,:), intent(in) :: end_matrix
        type(mpo) :: self
        integer :: i

        ! Consistancy Checks physical dimension
        if(size(start_matrix,3) /= pDim) then
            print *, "physical row dimension pDim1: ", size(start_matrix,3) ,&
                     " of start_matrix incorrect (it should be: ", pDim, " )"
            call abort()
        end if
        if(size(start_matrix,4) /= pDim) then
            print *, "physical column dimension pDim2: ", size(start_matrix,4) ,&
                     " of start_matrix incorrect (it should be: ", pDim, " )"
            call abort()
        end if
        if(size(middle_matrix,3) /= pDim) then
            print *, "physical row dimension pDim1: ", size(middle_matrix,3) ,&
                     " of middle_matrix incorrect (it should be: ", pDim, " )"
            call abort()
        end if
        if(size(middle_matrix,4) /= pDim) then
            print *, "physical column dimension pDim2: ", size(middle_matrix,4) ,&
                     " of middle_matrix incorrect (it should be: ", pDim, " )"
            call abort()
        end if
        if(size(end_matrix,3) /= pDim) then
            print *, "physical row dimension pDim1: ", size(end_matrix,3) ,&
                     " of end_matrix incorrect (it should be: ", pDim, " )"
            call abort()
        end if
        if(size(end_matrix,4) /= pDim) then
            print *, "physical column dimension pDim2: ", size(end_matrix,4) ,&
                     " of end_matrix incorrect (it should be: ", pDim, " )"
            call abort()
        end if
        if(size(start_matrix,1) /= 1) then
            print *, "start_matrix is not virtual row-vector (vDim1: ", size(start_matrix,1) , " )"
            call abort()
        end if
        if( size(start_matrix,2) /= size(middle_matrix,1))then
            print *, "Start_matrix (column dimension: ", size(start_matrix,2),&
                    " ) and middle_matrix ( row dimension: ", size(middle_matrix,1), " ) have odd virtual dimension"
            call abort()
        end if
        if( size(middle_matrix,1) /= size(middle_matrix,2))then
            print *, " Non-square virtual middle_matrix ( row dimension: ", size(middle_matrix,1),&
                    " column dimension: ", size(middle_matrix,2), " )"
            call abort()
        end if
        if( size(middle_matrix,2) /= size(end_matrix,1))then
            print *, "middle_matrix (column dimension: ", size(start_matrix,2),&
                    " ) and end_matrix ( row dimension: ", size(middle_matrix,1), " ) have odd virtual dimension"
            call abort()
        end if
        if(size(end_matrix,2) /= 1) then
            print *, "end_matrix is not virtual column-vector (vDim2: ", size(end_matrix,2) , " )"
            call abort()
        end if

        if(allocated(self%nodes)) deallocate(self%nodes)

        allocate(self%nodes(length))

        self%nodes(1) = mpo_node(start_matrix)

        do i = 2, (length - 1)
            self%nodes(i) = mpo_node(middle_matrix)
        end do

        self%nodes(length) = mpo_node(end_matrix)

    end function dedicated_new_mpo


    pure function single_pDim_node(self,l,pInd1,pInd2) result(mat)

        class(mpo), intent(in) :: self !Matrix product operator to get the node from
        integer, intent(in) :: l        ! physical site in the one-dimensional system
        integer, intent(in) :: pInd1        ! physical row dimension on site l to return
        integer, intent(in) :: pInd2        ! physical column dimension on site l to return
        complex(R_P), allocatable, dimension(:,:) :: mat  ! node matrix on site l for physical dimension d

        mat = self%nodes(l)%node(pInd1,pInd2)

    end function


    pure function all_pDim_node(self,l) result(mat)

        class(mpo), intent(in) :: self !Matrix product operator to get the node from
        integer, intent(in) :: l        ! physical site in the one-dimensional system
        complex(R_P), allocatable, dimension(:,:,:,:) :: mat  ! node matrix on site l for physical dimension d

        mat = self%nodes(l)%node()

    end function


    pure function mpo_size(self,length_ind,dim_ind) result(sze)

        type(mpo), intent(in) :: self
        integer, optional, intent(in) :: length_ind
        integer, optional, intent(in) :: dim_ind
        integer :: sze

        sze = size(self%nodes)

        if(present(length_ind)) then
            sze = size(self%nodes(length_ind))
            if(present(dim_ind)) then
                sze = size(self%nodes(length_ind),dim_ind)
            end if
        end if

    end function mpo_size



end module mpo_module
