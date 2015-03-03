module mps_class
    use mps_parameters
    use mps_node
    use IR_Precision
    implicit none
    private


    type, public :: mps
        private
        type(mp_node), allocatable, dimension(:) :: node_array  !allocatable due to compiler issues
        integer :: ind_vMax
        integer :: left_normalization = 0
        integer :: right_normalization
    contains
        procedure, private :: all_pDim_node
        procedure, private :: single_pDim_node
        generic :: node => all_pDim_node, single_pDim_node
        procedure :: move_left
        procedure :: move_right
        procedure :: left_normalize
        procedure :: right_normalize
        procedure, private :: real_divide
        procedure, private :: complex_divide
        generic :: divide => complex_divide, real_divide
        procedure :: normalize
        procedure, private :: print_node
        procedure, private :: print_all
        generic :: print => print_node, print_all
        final :: delete
    end type mps

    private :: dedicated_constructor_mps, constant_override_mps, constant_mps

    public :: random_mps, scalar_product, norm

    interface mps
        module procedure constant_mps, constant_override_mps, dedicated_constructor_mps
    end interface mps


contains


    function dedicated_constructor_mps(fun,x,or_pDim,or_vMax,or_length) result(state)

        integer, intent(in) :: or_pDim     !override physical dimension
        integer, intent(in) :: or_vMax     !override maximal virtual dimension
        integer, intent(in) :: or_length   !override one-dimensional system length
        complex(R_P), intent(in) :: x
        interface
            function fun(vDim1,vDim2,pDim,x) result(z)
                use IR_Precision
                integer, intent(in) :: vDim1
                integer, intent(in) :: vDim2
                integer, intent(in) :: pDim
                complex(R_P), intent(in) :: x
                complex(R_P) :: z
             end function
        end interface
        type(mps) :: state ! output matrix product state
        integer :: i, j
        integer, dimension(or_length+1) :: dims

        dims = 1
        allocate(state%node_array(or_length))
! Implement min search tip of Oli
        if (mod(or_length,2) == 0) then
            do i = 2,(or_length/2 + 1)
                dims(i) = or_pDim**(i - 1)
            end do
            do i = (or_length/2 + 2),or_length
                dims(i) = or_pDim**(or_length - i + 1)
            end do
        else
            do i = 2,(or_length + 1)/2
                dims(i) = or_pDim**(i - 1)
            end do
            do i = (or_length + 1)/2 + 1, or_length
                dims(i) = or_pDim**(or_length- i + 1)
            end do
        end if

        where (dims > or_vMax)
            dims = or_vMax
        end where

        do i = 1,or_length
            state%node_array(i) = mp_node(reshape((/(fun(1,1,1,x), j=1,dims(i) * dims(i + 1) * or_pDim)/), &
            (/dims(i),dims(i+1),or_pDim/)))
        end do

        state%right_normalization = or_length + 1
        state%ind_vMax = or_vMax

    end function


    function constant_override_mps(x,or_pDim,or_vMax,or_length) result(state)

        integer, intent(in) :: or_pDim
        integer, intent(in) :: or_vMax
        integer, intent(in) :: or_length
        complex(R_P), intent(in) :: x
        type(mps) :: state ! output matrix product state

        state = dedicated_constructor_mps(constant_fun,x,or_pDim,or_vMax,or_length)

    end function


    function constant_mps(x) result(state)

        complex(R_P), intent(in) :: x
        type(mps) :: state ! output matrix product state

        state = constant_override_mps(x,pDim,vMax,length)

    end function


    subroutine random_mps(state,or_pDim,or_vMax,or_length,x)

        type(mps), intent(inout) :: state
        integer, optional, intent(in) :: or_pDim
        integer, optional, intent(in) :: or_vMax
        integer, optional, intent(in) :: or_length
        complex(R_P), intent(in), optional :: x

        complex(R_P) :: factor = (1.0,0.0_R_P)
        integer :: local_pDim = pDim
        integer :: local_vMax = vMax
        integer :: local_length = length

        if(present(x)) factor = x
        if(present(or_pDim)) local_pDim = or_pDim
        if(present(or_vMax)) local_vMax = or_vMax
        if(present(or_length)) local_length = or_length

        state = dedicated_constructor_mps(random_fun,factor,local_pDim,local_vMax,local_length)

    end subroutine


    pure function single_pDim_node(self,l,d) result(mat)

        class(mps), intent(in) :: self !Matrix product state to get the node from
        integer, intent(in) :: l        ! physical site in the one-dimensional system
        integer, intent(in) :: d        ! physical dimension on site l to return
        complex(R_P), allocatable, dimension(:,:) :: mat  ! node matrix on site l for physical dimension d

        mat = self%node_array(l)%node(d)

    end function


    pure function all_pDim_node(self,l) result(mat)

        class(mps), intent(in) :: self !Matrix product state to get the node from
        integer, intent(in) :: l        ! physical site in the one-dimensional system
        complex(R_P), allocatable, dimension(:,:,:) :: mat  ! node matrix on site l for physical dimension d

        mat = self%node_array(l)%node()

    end function


    subroutine move_left(self)

        class(mps), intent(inout) :: self

        if(self%left_normalization < size(self%node_array) - 1) then
            call self%node_array(self%left_normalization+2)%left_multiply(&
                    self%node_array(self%left_normalization+1)%left_normalize())
            self%left_normalization = self%left_normalization + 1
            if( self%right_normalization == self%left_normalization + 1) self%right_normalization = self%right_normalization + 1
        end if

    end subroutine


    subroutine move_right(self)

        class(mps), intent(inout) :: self

        if(self%right_normalization > 2) then
            call self%node_array(self%right_normalization - 2)%right_multiply(&
                    self%node_array(self%right_normalization - 1)%right_normalize())
            self%right_normalization = self%right_normalization - 1
            if( self%left_normalization == self%right_normalization - 1) self%left_normalization = self%left_normalization - 1
        end if

    end subroutine


    subroutine left_normalize(self)

        class(mps), intent(inout) :: self
        integer :: i

        do i = 1, (size(self%node_array) - 1 - self%left_normalization)
            call self%move_left
        end do

    end subroutine


    subroutine right_normalize(self)

        class(mps), intent(inout) :: self
        integer :: i

        do i = 1, (self%right_normalization - 2)
            call self%move_right
        end do

    end subroutine


    function scalar_product(self, other_self) result(overlap)

        class(mps), intent(in) :: self
        class(mps), intent(in) :: other_self
        complex(R_P) :: overlap
        complex(R_P), dimension(self%ind_vMax,other_self%ind_vMax) :: coupling_mat
        integer :: i, self_vDim1, self_vDim2, other_self_vDim1, other_self_vDim2

        call self%node_array(1)%dims(self_vDim1,self_vDim2)
        call other_self%node_array(1)%dims(other_self_vDim1,other_self_vDim2)

        coupling_mat(1:self_vDim2,1:other_self_vDim2) = contract(self%node_array(1), other_self%node_array(1))

        do i = 2, size(self%node_array)
            call self%node_array(i)%dims(self_vDim1,self_vDim2)
            call other_self%node_array(i)%dims(other_self_vDim1,other_self_vDim2)
            coupling_mat(1:self_vDim2,1:self_vDim2) &
                    = contract(self%node_array(i), other_self%node_array(i), coupling_mat(1:self_vDim1,1:other_self_vDim1))
        end do

        overlap = coupling_mat(1,1)

    end function scalar_product


    function norm(self) result(sze)

        class(mps), intent(in) :: self
        real(R_P) :: sze

        sze = sqrt(real(scalar_product(self,self),R_P))

    end function


    subroutine complex_divide(self, denominator)

        class(mps), intent(inout) :: self
        complex(R_P), intent(in) :: denominator

        call self%node_array(1)%complex_divide(denominator)

    end subroutine


    subroutine real_divide(self, denominator)

        class(mps), intent(inout) :: self
        real(R_P), intent(in) :: denominator

        call self%complex_divide(cmplx(denominator,kind=R_P))

    end subroutine real_divide


    subroutine normalize(self)

        class(mps), intent(inout) :: self

        call self%real_divide(norm(self))

    end subroutine


    subroutine print_node(self,l)

        class(mps), intent(in) :: self
        integer, intent(in) :: l
        type(mp_node) :: nde

        nde = self%node_array(l)
        call nde%print

    end subroutine


    subroutine print_all(self)

        class(mps), intent(in) :: self
        type(mp_node) :: nde
        integer :: i

        do i = 1,size(self%node_array)
            nde = self%node_array(i)
            print *,""
            print *,""
            print *,"Node:", i
            call nde%print
        end do

    end subroutine


    elemental subroutine delete(self)

        type(mps), intent(inout) :: self

        deallocate(self%node_array)

    end subroutine

! HELPER FUNCTIONS

    function constant_fun(vDim1,vDim2,pDim,x) result(z)

        integer, intent(in) :: vDim1
        integer, intent(in) :: vDim2
        integer, intent(in) :: pDim
        complex(R_P), intent(in) :: x
        complex(R_P) :: z

        z = x

     end function

     function random_fun(vDim1,vDim2,pDim,x) result(z)

        integer, intent(in) :: vDim1
        integer, intent(in) :: vDim2
        integer, intent(in) :: pDim
        complex(R_P), intent(in) :: x
        complex(R_P) :: z
        real(R_P) :: z1, z2

        call random_number(z1)
        call random_number(z2)

        z = cmplx(z1, z2, kind = R_P) * x

     end function

end module mps_class
