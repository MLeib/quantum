module mps_node
    use IR_Precision
    use lapack_wrap
    implicit none
    private

    type, public :: mp_node
        private
        complex(R_P), allocatable, dimension(:,:,:) :: matrix
        integer :: vDim1, vDim2, pDim
    contains
        procedure :: dims
        procedure, private :: all_node
        procedure, private :: single_node
        generic :: node => all_node, single_node
        procedure :: left_normalize
        procedure :: right_normalize
        procedure :: left_multiply
        procedure :: right_multiply
        procedure :: complex_divide
        procedure :: print
        final :: delete
    end type mp_node

    public :: contract

    interface mp_node
        module procedure  mat_mp_node
    end interface mp_node

contains

    function mat_mp_node(mat) result(node)

        complex(R_P), dimension(:,:,:), intent(in) :: mat
        type(mp_node) :: node
        integer, dimension(3) :: shp

        shp = shape(mat)
        node%vDim1 = shp(1)
        node%vDim2 = shp(2)
        node%pDim = shp(3)
        allocate(node%matrix(node%vDim1,node%vDim2,node%pDim), source=mat)

    end function

    pure function all_node(self) result(mat)

        class(mp_node), intent(in) :: self
        complex(R_P), dimension(self%vDim1,self%vDim2,self%pDim) :: mat

        mat = self%matrix

    end function


    pure function single_node(self,d) result(mat)

        class(mp_node), intent(in) :: self
        integer, intent(in) :: d                        !physical dimension
        complex(R_P), dimension(self%vDim1,self%vDim2) :: mat

        mat = self%matrix(:,:,d)

    end function


    subroutine dims(self, vDim1, vDim2)

        class(mp_node), intent(in) :: self
        integer, intent(out) :: vDim1
        integer, intent(out) :: vDim2

        vDim1 = self%vDim1
        vDim2 = self%vDim2

    end subroutine


    subroutine print(self)

        class(mp_node), intent(in) :: self
        integer :: i, p

        print '(a,i3)', "     |", self%pDim
        print '(i3,a,i3)', self%vDim1, "-node-", self%vDim2

        do p = 1,self%pDim
            print *, ""
            print '(a,i3)', "physical dimension:", p
            do i = 1,self%vDim1
                print *, self%matrix(i,:,p)
            end do
        end do

    end subroutine


    elemental subroutine delete(self)

        type(mp_node), intent(inout) :: self

        deallocate(self%matrix)

    end subroutine


    function left_normalize(self) result(mat)

        class(mp_node), intent(inout) :: self
        complex(R_P), dimension(self%vDim2,self%vDim2) :: mat
        complex(R_P), dimension(self%pDim * self%vDim1, self%vDim2) :: u
        complex(R_P), dimension(self%vDim2,self%vDim2) :: s
        complex(R_P), dimension(self%vDim2,self%vDim2) :: vh

        call svd(left_concat(self%matrix),u,s,vh)
        self%matrix = left_deConcat(u,self%pDim)
        mat = matmul(s,vh)

    end function left_normalize


    function right_normalize(self) result(mat)

        class(mp_node), intent(inout) :: self
        complex(R_P), dimension(self%vDim1,self%vDim1) :: mat
        complex(R_P), dimension(self%vDim1,self%vDim1) :: u
        complex(R_P), dimension(self%vDim1,self%vDim1) :: s
        complex(R_P), dimension(self%vDim1,self%vDim2 * self%pDim) :: vh

        call svd(right_concat(self%matrix),u,s,vh)
        self%matrix = right_deConcat(vh,self%pDim)
        mat = matmul(u,s)

    end function right_normalize


    subroutine left_multiply(self, mat)

        class(mp_node), intent(inout) :: self
        complex(R_P), dimension(self%vDim1,self%vDim1), intent(in) :: mat
        integer :: i

        do i = 1, self%pDim
            self%matrix(:,:,i) = multiply(mat,self%matrix(:,:,i))
        end do

    end subroutine left_multiply


    subroutine right_multiply(self, mat)

        class(mp_node), intent(inout) :: self
        complex(R_P), dimension(self%vDim2,self%vDim2), intent(in) :: mat
        integer :: i

        do i = 1, self%pDim
            self%matrix(:,:,i) = multiply(self%matrix(:,:,i),mat)
        end do

    end subroutine right_multiply


    subroutine complex_divide(self,denominator)

        class(mp_node), intent(inout) :: self
        complex(R_P), intent(in) :: denominator

        self%matrix = self%matrix / denominator

    end subroutine complex_divide


    function contract(node1,node2,old_coupling_mat) result(new_coupling_mat)

        class(mp_node), intent(in) :: node1
        class(mp_node), intent(in) :: node2
        complex(R_P), optional, dimension(node1%vDim1,node2%vDim1) :: old_coupling_mat
        complex(R_P), dimension(node1%vDim2,node2%vDim2) :: new_coupling_mat
        integer :: i

        new_coupling_mat = (0.0,0.0_R_P)

        if(present(old_coupling_mat)) then
            do i = 1,node1%pDim
            new_coupling_mat = new_coupling_mat &
                    + multiply(hc_multiply(node1%matrix(:,:,i),old_coupling_mat),node2%matrix(:,:,i))
            end do
        else
            do i = 1,node1%pDim
                new_coupling_mat = new_coupling_mat &
                        + hc_multiply(node1%matrix(:,:,i),node2%matrix(:,:,i))
            end do
        end if

    end function contract


! HELPER FUNCTIONS

    function left_concat(mat3d) result(mat2d)

        complex(R_P), dimension(:,:,:), intent(in) :: mat3d
        complex(R_P), allocatable, dimension(:,:) :: mat2d
        integer :: vDim1, vDim2, pDim, i, j

        vDim1 = size(mat3d,1)
        vDim2 = size(mat3d,2)
        pDim = size(mat3d,3)
        allocate(mat2d(vDim1 * pDim, vDim2))

        mat2d = reshape([((mat3d(:,i,j),j=1,pDim),i=1,vDim2)],[vDim1 * pDim, vDim2])

    end function


    function left_deConcat(mat2d,pDim) result(mat3d)

        complex(R_P), dimension(:,:), intent(in) :: mat2d
        integer, intent(in) :: pDim
        complex(R_P), allocatable, dimension(:,:,:) :: mat3d
        integer :: vDim1, vDim2,i

        vDim1 = size(mat2d,1) / pDim
        vDim2 = size(mat2d,2)

        allocate(mat3d(vDim1,vDim2,pDim))

        do i = 1,pDim
            mat3d(:,:,i) = mat2d((((i-1) * vDim1)+1):(i*vDim1),:)
        end do

    end function


    function right_concat(mat3d) result(mat2d)

        complex(R_P), dimension(:,:,:), intent(in) :: mat3d
        complex(R_P), dimension(size(mat3d,1),size(mat3d,2) * size(mat3d,3)) :: mat2d
        integer :: vDim1, vDim2, pDim, i, j

        vDim1 = size(mat3d,1)
        vDim2 = size(mat3d,2)
        pDim = size(mat3d,3)

        mat2d = reshape([((mat3d(:,i,j),i=1,vDim2),j=1,pDim)],[vDim1 , vDim2 * pDim])

    end function


    function right_deConcat(mat2d,pDim) result(mat3d)

        complex(R_P), dimension(:,:), intent(in) :: mat2d
        integer, intent(in) :: pDim
        complex(R_P), dimension(size(mat2d,1),size(mat2d,2) / pDim,pDim) :: mat3d
        integer :: vDim2,i

        vDim2 = size(mat2d,2) / pDim

        do i = 1,pDim
            mat3d(:,:,i) = mat2d(:,(((i-1) * vDim2)+1):(i*vDim2))
        end do

    end function


end module mps_node
