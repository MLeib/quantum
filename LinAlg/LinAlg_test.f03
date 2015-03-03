module LinAlg_test
    use IR_Precision
    use random
    use lapack_wrap
    use gnufor2
    use node_class
    implicit none
    private

    public :: svd_precision_speed_test
    public :: multiply_precision_speed_test
    public :: gs_search_precision_speed_test
    public :: transpose_precision_speed_test

contains

    subroutine svd_precision_speed_test()

        implicit none

        complex(R_P), allocatable, dimension(:,:) :: mat, matcopy, u, s, vh, uhu, vhv
        real(R_P) :: timer_start, timer_stop, duration
        real(R_P), parameter :: max_duration = 10_R_P
        integer, parameter :: average_cycles = 40 , mat_start_size = 100, mat_stride = 10
        real(R_P), dimension(200) :: avgTime = 0_R_P, avgRepPrec = 0_R_P, avgUPrec = 0_R_P, avgVPrec = 0_R_P
        real(R_P), allocatable, dimension(:) :: mat_size
        integer :: i, j, k, v1, v2

        j= 0
        call random_seed()
        duration = 0.0_R_P

        do while (duration < max_duration)

            j = j + 1

            do i = 1,average_cycles

            v1 = random_integer(mat_start_size + mat_stride * j - 2)
            v2 = mat_start_size + mat_stride * j - v1

            allocate(mat(v1,v2))
            allocate(matcopy(v1,v2))
            allocate(u(v1,min(v1,v2)))
            allocate(s(min(v1,v2),min(v1,v2)))
            allocate(vh(min(v1,v2),v2))
            allocate(uhu(min(v1,v2),min(v1,v2)))
            allocate(vhv(min(v1,v2),min(v1,v2)))

                call normal(mat)
                matcopy = mat

                call cpu_time(timer_start)
                call svd(mat,u,s,vh)
                call cpu_time(timer_stop)

                uhu = matmul(conjg(transpose(u)),u)
                vhv = matmul(vh,conjg(transpose(vh)))
                forall( k = 1:min(v1,v2)) uhu(k,k) = uhu(k,k) - 1_R_P
                forall( k = 1:min(v1,v2)) vhv(k,k) = vhv(k,k) - 1_R_P

                avgTime(j) = avgTime(j) + (timer_stop - timer_start)
                avgRepPrec(j) = avgRepPrec(j) + frobNorm(matmul(u,matmul(s,vh)) - matcopy) / frobNorm(matcopy)
                avgUPrec(j) = avgUPrec(j) + frobNorm(uhu) / sqrt(real(min(v1,v2),R_P))
                avgVPrec(j) = avgVPrec(j) + frobNorm(vhv) / sqrt(real(min(v2,v2),R_P))

            deallocate(mat)
            deallocate(matcopy)
            deallocate(u)
            deallocate(s)
            deallocate(vh)
            deallocate(uhu)
            deallocate(vhv)

            end do

            duration = avgTime(j)
            print *, duration
            avgTime(j) = avgTime(j) / average_cycles
            avgRepPrec(j) = avgRepPrec(j) / average_cycles
            avgUPrec(j) = avgUPrec(j) / average_cycles
            avgVPrec(j) = avgVPrec(j) / average_cycles

        end do

        allocate(mat_size(j))
        mat_size = [((mat_start_size + i * mat_stride), i = 1,j)]
        call plot(mat_size,avgTime(1:j),terminal='ps',filename='avgTime_svd.ps')
        call plot(mat_size,avgRepPrec(1:j),terminal='ps',filename='avgRepPrec_svd.ps')
        call plot(mat_size,avgUPrec(1:j),terminal='ps',filename='avgUPrec_svd.ps')
        call plot(mat_size,avgVPrec(1:j),terminal='ps',filename='avgVPrec_svd.ps')
        print *, mat_size(j), avgTime(j)

    end subroutine


    subroutine multiply_precision_speed_test()

        implicit none

        complex(R_P), allocatable, dimension(:,:) :: matA, matAcopy, matB, matBcopy, matC, matCcopy
        real(R_P) :: timer_start, timer_stop, duration
        real(R_P), parameter :: max_duration = 10_R_P
        integer, parameter :: average_cycles = 40 , mat_start_size = 100, mat_stride = 10
        real(R_P), dimension(200) :: avgTime_BI = 0_R_P, avgTime_LB = 0_R_P, avgRepPrec = 0_R_P
        real(R_P), allocatable, dimension(:) :: mat_size
        integer :: i, j

        j= 0
        call random_seed()
        duration = 0.0_R_P

        do while (duration < max_duration)

            j = j + 1

            allocate(matA(mat_start_size + j * mat_stride, mat_start_size + j * mat_stride))
            allocate(matAcopy(mat_start_size + j * mat_stride, mat_start_size + j * mat_stride))
            allocate(matB(mat_start_size + j * mat_stride, mat_start_size + j * mat_stride))
            allocate(matBcopy(mat_start_size + j * mat_stride, mat_start_size + j * mat_stride))
            allocate(matC(mat_start_size + j * mat_stride, mat_start_size + j * mat_stride))
            allocate(matCcopy(mat_start_size + j * mat_stride, mat_start_size + j * mat_stride))

            do i = 1,average_cycles

                call normal(matA)
                matAcopy = matA
                call normal(matB)
                matBcopy = matB

                call cpu_time(timer_start)
                matC = multiply(matA,matB)
                call cpu_time(timer_stop)
                avgTime_LB(j) = avgTime_LB(j) + (timer_stop - timer_start)


                call cpu_time(timer_start)
                matCcopy = matmul(matAcopy,matBcopy)
                call cpu_time(timer_stop)
                avgTime_BI(j) = avgTime_BI(j) + (timer_stop - timer_start)

                avgRepPrec(j) = avgRepPrec(j) + frobNorm(matC - matCcopy) / frobNorm(matCcopy)

            end do

            deallocate(matA)
            deallocate(matAcopy)
            deallocate(matB)
            deallocate(matBcopy)
            deallocate(matC)
            deallocate(matCcopy)
            duration = avgTime_LB(j) + avgTime_BI(j)

            print *, duration
            avgTime_LB(j) = avgTime_LB(j) / average_cycles
            avgTime_BI(j) = avgTime_BI(j) / average_cycles
            avgRepPrec(j) = avgRepPrec(j) / average_cycles

        end do

        allocate(mat_size(j))
        mat_size = [((mat_start_size + i * mat_stride), i = 1,j)]
        call plot(mat_size,avgTime_LB(1:j),mat_size,avgTime_BI(1:j),terminal='ps',filename='avgTime_multiply.ps')
        call plot(mat_size,avgRepPrec(1:j),terminal='ps',filename='avgRepPrec_multiply.ps')

    end subroutine


    subroutine gs_search_precision_speed_test()

        complex(R_P), allocatable, dimension(:,:) :: mat
        complex(R_P), allocatable, dimension(:) :: evec
        real(R_P) :: eval
        integer :: i, j
        real(R_P) :: timer_start, timer_stop, duration
        real(R_P), parameter :: max_duration = 10_R_P
        integer, parameter :: average_cycles = 40 , mat_start_size = 100, mat_stride = 10
        real(R_P), dimension(200) :: avgTime = 0_R_P, avgPrec = 0_R_P
        real(R_P), allocatable, dimension(:) :: mat_size

        j = 0
        call random_seed()
        duration = 0.0_R_P

        do while (duration < max_duration)

            j = j + 1
            allocate(mat(mat_start_size + j * mat_stride,mat_start_size + j * mat_stride))
            allocate(evec(mat_start_size + j * mat_stride))

            do i = 1,average_cycles
                call normal(mat)

                mat = mat + conjg(transpose(mat))

                call cpu_time(timer_start)
                call gs_search(mat,eval,evec)
                call cpu_time(timer_stop)

                avgTime(j) = avgTime(j) + (timer_stop - timer_start)
                avgPrec(j) = avgPrec(j) + real(sqrt(sum((matmul(mat,evec) - eval * evec)**2)),kind=R_P)

            end do

            duration = avgTime(j)
            print *, duration

            deallocate(mat)
            deallocate(evec)
            avgTime(j) = avgTime(j) / average_cycles
            avgPrec(j) = avgPrec(j) / average_cycles

        end do

        allocate(mat_size(j))
        mat_size = [((mat_start_size + i * mat_stride), i = 1,j)]
        call plot(mat_size,avgTime(1:j),terminal='ps',filename='avgTime_gs_search.ps')
        call plot(mat_size,avgPrec(1:j),terminal='ps',filename='avgPrec_gs_search.ps')

    end subroutine gs_search_precision_speed_test


    subroutine transpose_precision_speed_test()

        complex(R_P), dimension(10,20) :: mat, copy_mat
        complex(R_P), dimension(20,10) :: trans_mat, trans_copy_mat
        real(R_P) :: prec
        integer, dimension(2,1), parameter :: transpose_pair = reshape([1,2],[2,1])
        type(node) :: nde
        integer, dimension(2) :: shp

        call normal(mat)
        copy_mat = mat
        nde = node(mat)

        trans_copy_mat = transpose(copy_mat)
        nde = transpose(nde,transpose_pair)

        shp = nde%get_shape()
        trans_mat = reshape(nde%get_elements(),shp)

        prec = sqrt(sum(abs(trans_mat-trans_copy_mat)**2))

        print *, prec

    end subroutine transpose_precision_speed_test


! HELPER FUNCTIONS

    function frobNorm(mat) result(norm)

        complex(R_P), dimension(:,:), intent(in) :: mat
        real(R_P) :: norm

        norm = sqrt(sum(abs(mat)**2))

    end function


end module LinAlg_test
