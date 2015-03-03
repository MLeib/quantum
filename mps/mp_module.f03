module mp_module
!    use node_class
!    use mps_class
!    use mpo_module
!    use mps_parameters
!    use IR_Precision
!    implicit none
!    private
!
!    public :: create_L
!    public :: create_R
!
!contains
!
!    function create_L(l,mp_operator,mp_state,old_L) result(new_L)
!
!        integer, intent(in) :: l
!        type(mpo), intent(in) :: mp_operator
!        type(mps), intent(in) :: mp_state
!        type(node), optional,  intent(in) :: old_L
!        type(node) :: new_L
!
!        type(node) :: mps_node_up, mps_node_down
!        type(node) :: mpo_node
!        complex(R_P), dimension(1,size(mp_state%node(1),2),pDim) :: boundary_state_tensor
!        complex(R_P), dimension(1,size(mp_operator%node(1),2),pDim,pDim) :: boundary_operator_tensor
!
!        integer, dimension(2,1), parameter :: first_contract_boundary = reshape([2,2],[2,1])
!        integer, dimension(2,1), parameter :: second_contract_boundary = reshape([3,2],[2,1])
!        integer, dimension(2,1), parameter :: first_contract_bulk = reshape([1,1],[2,1])
!        integer, dimension(2,2), parameter :: second_contract_bulk = reshape([2,3,3,1],[2,2])
!        integer, dimension(2,2), parameter :: third_contract_bulk = reshape([4,3,2,1],[2,2])
!
!        if(present(old_L)) then
!            mps_node_up = node(conjg(mp_state%node(l)))
!            mps_node_down = node(mp_state%node(l))
!            mpo_node = node(mp_operator%node(l))
!            new_L = merge(mps_node_up,old_L,first_contract_bulk)
!            new_L = merge(new_L,mpo_node,second_contract_bulk)
!            new_L = merge(new_L,mps_node_down,third_contract_bulk)
!        else
!            boundary_state_tensor = conjg(mp_state%node(1))
!            mps_node_up = node(boundary_state_tensor(1,:,:))
!            boundary_state_tensor = mp_state%node(1)
!            mps_node_down = node(boundary_state_tensor(1,:,:))
!            boundary_operator_tensor = mp_operator%node(1)
!            mpo_node = node(boundary_operator_tensor(1,:,:,:))
!            new_L = merge(mps_node_up,mpo_node,first_contract_boundary)
!            new_L = merge(new_L,mps_node_down,second_contract_boundary)
!        end if
!
!    end function create_L
!
!
!    function create_R(l,mp_operator,mp_state,old_R) result(new_R)
!
!        integer, intent(in) :: l
!        type(mpo), intent(in) :: mp_operator
!        type(mps), intent(in) :: mp_state
!        type(node), optional,  intent(in) :: old_R
!        type(node) :: new_R
!
!        type(node) :: mps_node_up, mps_node_down
!        type(node) :: mpo_node
!        complex(R_P), dimension(size(mp_state%node(length),1),1,pDim) :: boundary_state_tensor
!        complex(R_P), dimension(size(mp_operator%node(length),1),1,pDim,pDim) :: boundary_operator_tensor
!
!        integer, dimension(2,1), parameter :: first_contract_boundary = reshape([2,2],[2,1])
!        integer, dimension(2,1), parameter :: second_contract_boundary = reshape([3,2],[2,1])
!        integer, dimension(2,1), parameter :: first_contract_bulk = reshape([2,1],[2,1])
!        integer, dimension(2,2), parameter :: second_contract_bulk = reshape([2,3,3,2],[2,2])
!        integer, dimension(2,2), parameter :: third_contract_bulk = reshape([4,3,2,2],[2,2])
!
!        if(present(old_R)) then
!            mps_node_up = node(conjg(mp_state%node(l)))
!            mps_node_down = node(mp_state%node(l))
!            mpo_node = node(mp_operator%node(l))
!            new_R = merge(mps_node_up,old_R,first_contract_bulk)
!            new_R = merge(new_R,mpo_node,second_contract_bulk)
!            new_R = merge(new_R,mps_node_down,third_contract_bulk)
!        else
!            boundary_state_tensor = conjg(mp_state%node(length))
!            mps_node_up = node(boundary_state_tensor(:,1,:))
!            boundary_state_tensor = mp_state%node(length)
!            mps_node_down = node(boundary_state_tensor(:,1,:))
!            boundary_operator_tensor = mp_operator%node(length)
!            mpo_node = node(boundary_operator_tensor(:,1,:,:))
!            new_R = merge(mps_node_up,mpo_node,first_contract_boundary)
!            new_R = merge(new_R,mps_node_down,second_contract_boundary)
!        end if
!
!    end function create_R
!
!
!
!    subroutine variational_gs_search(hamilton,groundstate,vMax,gs_energy)
!
!        type(mpo), intent(in) :: hamilton
!        type(mps), intent(out) :: groundstate
!        integer, intent(in) :: vMax
!        real(R_P), intent(out) :: gs_energy
!
!        integer :: length = size(hamilton), pDim = size(hamilton,1,3)
!        type(node), dimension(1:(length-1)) :: l_array
!        type(node), dimension(2:length) :: r_array
!        type(node) :: opti_node, hamilton_node
!        integer :: i
!
!        call random_mps(groundstate,pDim,vMax,length)   ! take random initial guess for the groundstate
!        call right_normalize(groundstate)  ! right normalize for initial right sweep
!
!        ! Prepare and save R-matrices for first sweep
!        r_array(length) = create_R(length,hamilton,groundstate)
!        do i = length - 1,2,-1
!            r_array(i) = create_r(i,hamilton,groundstate,r_array(i+1))
!        end do
!
!        !First right-sweep
!            !left boundary
!        hamilton_node = node(hamilton%node(1))
!        opti_node = merge(hamilton_node,r_array(2),[1,2])
!
!
!        !bulk sweep
!        do i = 1,length
!        end do
!
!    end subroutine variational_gs_search

end module mp_module
