module random
use IR_Precision
implicit none
private

public :: random_cmplx, random_integer
public :: normal

interface random_cmplx
    module procedure random_cmplx_sparse, random_cmplx_dense
end interface random_cmplx

interface normal
    module procedure normal_rank2, normal_rank1, normal_rank0
end interface normal

contains

    subroutine normal_rank2(cplxMat,sparse)

        complex(R_P), dimension(:,:), intent(inout) :: cplxMat
        real(R_P), optional, intent(in) :: sparse
        integer(I_P) :: i, j

        if(present(sparse)) then
            do i = 1,size(cplxMat,2)
                do j = 1,size(cplxMat,1)
                    cplxMat(j,i) = random_cmplx(sparse)
                end do
            end do
        else
            do i = 1,size(cplxMat,2)
                do j = 1,size(cplxMat,1)
                    cplxMat(j,i) = random_cmplx()
                end do
            end do

        end if

    end subroutine


    subroutine normal_rank1(cplxVec,sparse)

        complex(R_P), dimension(:), intent(inout) :: cplxVec
        real(R_P), optional, intent(in) :: sparse
        complex(R_P), dimension(size(cplxVec),1) :: cplxMat

        if(present(sparse))then
            call normal_rank2(cplxMat,sparse)
        else
            call normal_rank2(cplxMat)
        end if

        cplxVec = cplxMat(:,1)

    end subroutine normal_rank1


    subroutine normal_rank0(cplx,sparse)

        complex(R_P),  intent(inout) :: cplx
        real(R_P), optional, intent(in) :: sparse
        complex(R_P), dimension(1,1) :: cplxMat

        if(present(sparse))then
            call normal_rank2(cplxMat,sparse)
        else
            call normal_rank2(cplxMat)
        end if

        cplx = cplxMat(1,1)

    end subroutine normal_rank0


    function random_cmplx_sparse(sparse) result(cplx)

        real(R_P), intent(in) :: sparse
        complex(R_P) :: cplx
        real(R_P) :: u1, u2, x, y

        call random_number(u1)

        if(u1 < sparse)then
            call random_number(u2)
            x = sqrt(- 2_R_P * log(u1)) * sin(2_R_P * PI * u2)
            y = sqrt(- 2_R_P * log(u1)) * cos(2_R_P * PI * u2)
            cplx = cmplx(x,y, kind=R_P)
        else
            cplx = (0.0,0.0_R_P)
        end if

    end function random_cmplx_sparse


    function random_cmplx_dense() result(cplx)

        complex(R_P) :: cplx

        cplx = random_cmplx_sparse(1.0_R_P)

    end function random_cmplx_dense


    function random_integer(max) result(nat)

        integer, intent(in) :: max
        integer :: nat
        real(R_P) :: r

        call random_number(r)

        nat = floor(r * max) + 1

     end function random_integer

end module random

