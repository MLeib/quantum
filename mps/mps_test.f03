module mps_test
use IR_Precision
use mps_class
use mpo_module
implicit none
private

public :: zero_constructor_test
public :: normalize_test
public :: hamilton

contains

    function zero_constructor_test() result(success)

        logical :: success
        type(mps) :: state

        success = .false.
        !L even
        !d^n smaller D
        state = mps((1.0,2.0_R_P),2,5,4)
        if(.not. all(shape(state%node(1)) == (/1,2,2/))) then
            print *, "L even, d^n smaller D : failed to correctly initialize node 1"
            print *, "node 1: ", state%node(1)
            return
        end if
        if(.not. all(shape(state%node(2)) == (/2,4,2/))) then
            print *, "L even, d^n smaller D : failed to correctly initialize node 2"
            print *, "node 2 shape: ", shape(state%node(2))
            print *, "node 2: ", state%node(2)
            return
        end if
        if(.not. all(shape(state%node(3)) == (/4,2,2/))) then
            print *, "L even, d^n smaller D : failed to correctly initialize node 3"
            print *, "node 3: ", state%node(3)
            return
        end if
        if(.not. all(shape(state%node(4)) == (/2,1,2/))) then
            print *, "L even, d^n smaller D : failed to correctly initialize node 4"
            print *, "node 4: ", state%node(4)
            return
        end if
        !d^n greater D
        state = mps((1.0,2.0_R_P),2,3,4)
        if(.not. all(shape(state%node(1)) == (/1,2,2/))) then
            print *, "L even, d^n greater D : failed to correctly initialize node 1"
            print *, "node 1: ", state%node(1)
            return
        end if
        if(.not. all(shape(state%node(2)) == (/2,3,2/))) then
            print *, "L even, d^n greater D : failed to correctly initialize node 2"
            print *, "node 2: ", state%node(2)
            return
        end if
        if(.not. all(shape(state%node(3)) == (/3,2,2/))) then
            print *, "L even, d^n greater D : failed to correctly initialize node 3"
            print *, "node 3: ", state%node(3)
            return
        end if
        if(.not. all(shape(state%node(4)) == (/2,1,2/))) then
            print *, "L even, d^n greater D : failed to correctly initialize node 4"
            print *, "node 4: ", state%node(4)
            return
        end if
        !L odd
        !d^n smaller D
        state = mps((1.0,2.0_R_P),3,4,3)
        if(.not. all(shape(state%node(1)) == (/1,3,3/))) then
            print *, "L odd, d^n smaller D : failed to correctly initialize node 1"
            print *, "node 1 shape: ", shape(state%node(1))
            print *, "node 1: ", state%node(1)
            return
        end if
        if(.not. all(shape(state%node(2)) == (/3,3,3/))) then
            print *, "L odd, d^n smaller D : failed to correctly initialize node 2"
            print *, "node 2: ", state%node(2)
            return
        end if
        if(.not. all(shape(state%node(3)) == (/3,1,3/))) then
            print *, "L odd, d^n smaller D : failed to correctly initialize node 3"
            print *, "node 3: ", state%node(3)
            return
        end if
        !d^n greater D
        state = mps((1.0,2.0_R_P),3,2,3)
        if(.not. all(shape(state%node(1)) == (/1,2,3/))) then
            print *, "L odd, d^n greater D : failed to correctly initialize node 1"
            print *, "node 1: ", state%node(1)
            return
        end if
        if(.not. all(shape(state%node(2)) == (/2,2,3/))) then
            print *, "L odd, d^n greater D : failed to correctly initialize node 2"
            print *, "node 2: ", state%node(2)
            return
        end if
        if(.not. all(shape(state%node(3)) == (/2,1,3/))) then
            print *, "L odd, d^n greater D : failed to correctly initialize node 3"
            print *, "node 3: ", state%node(3)
            return
        end if
        success = .true.
    end function zero_constructor_test


    subroutine normalize_test()

        real(R_P) :: left_prec, right_prec, acu_left_prec, acu_right_prec
        complex(R_P) :: rep_left_prec, rep_right_prec, acu_rep_left_prec, acu_rep_right_prec

        call individual_normalize_test(4,20,7, left_prec, right_prec,rep_left_prec,rep_right_prec)
            acu_left_prec = left_prec
            acu_right_prec = right_prec
            acu_rep_left_prec = rep_left_prec
            acu_rep_right_prec = rep_right_prec
        call individual_normalize_test(2,10,20, left_prec, right_prec,rep_left_prec,rep_right_prec)
            acu_left_prec = acu_left_prec + left_prec
            acu_right_prec = acu_right_prec + right_prec
            acu_rep_left_prec = acu_rep_left_prec + rep_left_prec
            acu_rep_right_prec = acu_rep_right_prec + rep_right_prec
        call individual_normalize_test(10,15,20, left_prec, right_prec,rep_left_prec,rep_right_prec)
            acu_left_prec = acu_left_prec + left_prec
            acu_right_prec = acu_right_prec + right_prec
            acu_rep_left_prec = acu_rep_left_prec + rep_left_prec
            acu_rep_right_prec = acu_rep_right_prec + rep_right_prec

        print *, "left normalization precision: ", acu_left_prec / 3_R_P
        print *, "right normalization precision: ", acu_right_prec / 3_R_P
        print *, "left reproduction precision: ", acu_rep_left_prec / 3_R_P
        print *, "right reproduction precision: ", acu_rep_right_prec / 3_R_P

    contains

        subroutine individual_normalize_test(local_pDim,local_vMax,local_length,left_prec,right_prec,rep_left_prec,rep_right_prec)

            integer, intent(in) :: local_pDim
            integer, intent(in) :: local_vMax
            integer, intent(in) :: local_length
            real(R_P), intent(out) :: left_prec
            real(R_P), intent(out) :: right_prec
            complex(R_P), intent(out) :: rep_left_prec
            complex(R_P), intent(out) :: rep_right_prec
            type(mps) :: left_state, right_state, left_state_copy, right_state_copy
            complex(R_P), allocatable, dimension(:,:) :: left_mat, right_mat
            integer :: i, j, vDim1, vDim2

            left_prec = 0.0_R_P
            right_prec = 0.0_R_P

            call random_mps(left_state,local_pDim,local_vMax,local_length)
            call left_state%normalize()
            left_state_copy = left_state
            call left_state%left_normalize
            rep_left_prec = scalar_product(left_state_copy,left_state) - (1.0,0.0_R_P)
            call random_mps(right_state,local_pDim,local_vMax,local_length)
            call right_state%normalize()
            right_state_copy = right_state
            rep_right_prec = scalar_product(right_state_copy,right_state) - (1.0,0.0_R_P)
            call right_state%right_normalize

            do i = 1, local_length

                vDim2 = size(left_state%node(i,1),2)
                allocate(left_mat(vDim2,vDim2))
                left_mat = (0.0, 0.0_R_P)
                vDim1 = size(right_state%node(i,1),1)
                allocate(right_mat(vDim1,vDim1))
                right_mat = (0.0, 0.0_R_P)

                do j = 1,local_pDim
                    if( i /= local_length) then
                        left_mat = left_mat + matmul(conjg(transpose(left_state%node(i,j))),left_state%node(i,j))
                    end if
                    if (i /= 1) then
                        right_mat = right_mat + matmul(right_state%node(i,j),conjg(transpose(right_state%node(i,j))))
                    end if
                end do

                if( i /= local_length) then
                    forall( j = 1:vDim2) left_mat(j,j) = left_mat(j,j) - (1.0 , 0.0_R_P)
                    left_prec = left_prec + frobNorm(left_mat) / sqrt(real(vDim2,R_P))
                end if
                if( i /= 1) then
                    forall( j = 1:vDim1) right_mat(j,j) = right_mat(j,j) - (1.0 , 0.0_R_P)
                    right_prec = right_prec + frobNorm(right_mat) / sqrt(real(vDim1,R_P))
                end if
                deallocate(left_mat)
                deallocate(right_mat)

            end do

        end subroutine

    end subroutine


    function hamilton(h,j,j_z) result(hamilton_op)

        real(R_P), intent(in) :: h
        real(R_P), intent(in) :: j
        real(R_P), intent(in) :: j_z
        type(mpo) :: hamilton_op

        complex(R_P), dimension(1,5,2,2) :: start_matrix
        complex(R_P), dimension(5,5,2,2) :: middle_matrix
        complex(R_P), dimension(5,1,2,2) :: end_matrix


            start_matrix(1,1,:,:) = -h * SIGMA_Z
            start_matrix(1,2,:,:) = (j / 2_R_P) * SIGMA_M
            start_matrix(1,3,:,:) = (j / 2_R_P) * SIGMA_P
            start_matrix(1,4,:,:) = j_z * SIGMA_Z
            start_matrix(1,5,:,:) = ID_2

            end_matrix(1,1,:,:) = ID_2
            end_matrix(2,1,:,:) = SIGMA_P
            end_matrix(3,1,:,:) = SIGMA_M
            end_matrix(4,1,:,:) = SIGMA_Z
            end_matrix(5,1,:,:) = -h * SIGMA_Z

            middle_matrix = (0.0,0.0_R_P)

            middle_matrix(:,1,:,:) = end_matrix(:,1,:,:)
            middle_matrix(5,:,:,:) = start_matrix(1,:,:,:)

            hamilton_op = mpo(start_matrix,middle_matrix,end_matrix)

    end function hamilton


! HELPER FUNCTIONS

    function frobNorm(mat) result(norm)

        complex(R_P), dimension(:,:), intent(in) :: mat
        real(R_P) :: norm

        norm = sqrt(sum(abs(mat)**2))

    end function


end module mps_test
