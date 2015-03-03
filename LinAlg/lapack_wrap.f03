module lapack_wrap
    use IR_Precision
    use random
    implicit none
    private

    public :: svd
    public :: multiply, hc_multiply
    public :: gs_search
    public :: zndrv1

contains

subroutine svd(mat,u,s,vh)

    complex(R_P), dimension(:,:), intent(in) :: mat
    complex(R_P), dimension(size(mat,1),min(size(mat,1),size(mat,2))), intent(out) :: u
    complex(R_P), dimension(min(size(mat,1),size(mat,2)),min(size(mat,1),size(mat,2))), intent(out) :: s
    complex(R_P), dimension(min(size(mat,1),size(mat,2)),size(mat,2)) :: vh

    real(R_P), dimension(min(size(mat,1),size(mat,2))) :: s_vec
    complex(R_P), dimension(1) :: work
    complex(R_P), allocatable, dimension(:) :: work_array
    real(R_P), dimension(5*min(size(mat,1),size(mat,2))) :: rwork
    integer :: info, lwork, i


    call zgesvd('S','S',size(mat,1),size(mat,2),mat,size(mat,1),s_vec,u,size(u,1),vh,size(vh,1),work,-1,rwork,info)
    lwork = int(work(1))
    allocate(work_array(lwork))
    call zgesvd('S','S',size(mat,1),size(mat,2),mat,size(mat,1),s_vec,u,size(u,1),vh,size(vh,1),work_array,lwork,rwork,info)

    s = (0.0, 0.0_R_P)
    do i = 1,min(size(mat,1),size(mat,2))
        s(i,i) = cmplx(s_vec(i),kind=R_P)
    end do

end subroutine


function multiply(matA,matB) result(matC)

    complex(R_P), dimension(:,:), intent(in) :: matA
    complex(R_P), dimension(:,:), intent(in) :: matB
    complex(R_P), dimension(size(matA,1),size(matB,2)) :: matC

    call zgemm('n','n',size(matA,1),size(matB,2), size(matA,2),(1.0,0.0_R_P),matA,size(matA,1),matB,size(matB,1),&
             (0.0,0.0_R_P), matC, size(matA,1))

end function


function hc_multiply(matA,matB) result(matC)

    complex(R_P), dimension(:,:), intent(in) :: matA
    complex(R_P), dimension(:,:), intent(in) :: matB
    complex(R_P), dimension(size(matA,2),size(matB,2)) :: matC

    call zgemm('c','n',size(matA,2),size(matB,2), size(matA,1),(1.0,0.0_R_P),matA,size(matA,1),matB,size(matB,1),&
             (0.0,0.0_R_P), matC, size(matA,2))

end function


subroutine gs_search(mat,e_Val,e_Vec)

    complex(R_P), dimension(:,:), intent(in) :: mat
    real(R_P), intent(out) :: e_Val
    complex(R_P), dimension(:), intent(out) :: e_Vec

    complex(R_P), dimension(size(mat,1),size(mat,1)) :: matcopy
    integer :: m ! number of eigenvalues found
    integer, dimension(2) :: isuppz ! internal mumbo jumbo
    complex(R_P) :: opt_work
    complex(R_P), allocatable, dimension(:) :: work
    real(R_P) :: opt_rwork
    real(R_P), allocatable, dimension(:) :: rwork
    integer :: opt_iwork
    integer, allocatable, dimension(:) :: iwork
    integer :: info

    matcopy = mat ! make copy of mat to enforce intent(in) attribute

    ! optimal workspace querry
    call zheevr('V','I','U',size(mat,1),matcopy,size(mat,1),0.0_R_P,0.0_R_P,1,1,0.0_R_P,m,e_Val,e_Vec,size(mat,1),isuppz,&
            opt_work,-1,opt_rwork,1,opt_iwork,1,info)

    if(info == 0) then
        allocate(work(int(opt_work)))
        allocate(rwork(int(opt_rwork)))
        allocate(iwork(opt_iwork))
        call zheevr('V','I','U',size(mat,1),matcopy,size(mat,1),0.0_R_P,0.0_R_P,1,1,0.0_R_P,m,e_Val,e_Vec,size(mat,1),isuppz,&
            work,size(work),rwork,size(rwork),iwork,size(iwork),info)
    end if

end subroutine gs_search


subroutine zndrv1(mat,gs_vec,gs_val)

    complex(R_P), dimension(:,:), intent(in) :: mat
    complex(R_P), dimension(:), intent(inout) :: gs_vec
    real(R_P), intent(out) :: gs_val

    integer, parameter :: nev = 1, ncv = 10
    integer :: iparam(11), ipntr(14)
    logical :: select(ncv)
    complex(R_P) :: d(nev+1), v(size(mat,1),ncv), workd(3*size(mat,1)), workev(3*ncv), resid(size(mat,1))
    complex(R_P) :: workl(3*ncv*ncv+5*ncv)
    real(R_P) :: rwork(ncv)

    character(len=1) :: bmat = 'I'
    character(len=2) :: which = 'SR'
    integer, parameter :: ishfts = 1, maxitr = 300, mode = 1
    integer :: ido, lworkl, info, n, ldv, ierr
    complex(R_P) :: sigma
    real(R_P) :: tol
    logical, parameter :: rvec = .true.

    ierr = 0
    info = 1
    sigma = (0.0,0.0_R_P)
    tol = 0.0_R_P
    iparam = 0
    ipntr = 0
    select  = .true.
    d = (0.0,0.0_R_P)
    v = (0.0,0.0_R_P)
    workd = (0.0,0.0_R_P)
    workev = (0.0,0.0_R_P)
    workl = (0.0,0.0_R_P)
    lworkl = 3 * ncv**2+5*ncv
    rwork = 0.0_R_P
    resid = gs_vec
    n = size(mat,1)
    ldv = n
    ido = 0
    iparam(1) = ishfts
    iparam(3) = maxitr
    iparam(7) = mode

    do while (ido == 0 .or. ido == 1 .or. ido == -1)

         call znaupd( ido, bmat, n, which, nev, tol, resid, ncv,&
                v, ldv, iparam, ipntr, workd, workl, lworkl,&
                rwork,info )

         if (ido .eq. -1 .or. ido .eq. 1) then

            call av(workd(ipntr(1):(ipntr(1)+n-1)), workd(ipntr(2):(ipntr(2)+n-1)),mat)

         end if
    end do

    if ( info .lt. 0 ) then

        print *, ' '
        print *, ' Error with _naupd, info = ', info
        print *, ' Check the documentation of _naupd'
        print *, ' '

    else

        call zneupd (rvec, 'A', select, d, v, ldv, sigma,&
             workev, bmat, n, which, nev, tol, resid, ncv,&
             v, ldv, iparam, ipntr, workd, workl, lworkl, &
             rwork, ierr)

        if ( ierr .ne. 0) then

            print *, ' '
            print *, ' Error with _neupd, info = ', ierr
            print *, ' Check the documentation of _neupd. '
            print *, ' '

        end if

        gs_vec = v(:,1)
        gs_val = realpart(d(1))

!c        %-------------------------------------------%
!c        | Print additional convergence information. |
!c        %-------------------------------------------%

         if ( info .eq. 1) then
             print *, ' '
             print *, ' Maximum number of iterations reached.'
             print *, ' '
         else if ( info .eq. 3) then
             print *, ' '
             print *, ' No shifts could be applied during implicit Arnoldi update, try increasing NCV.'
             print *, ' '
         end if

         print *, ' '
         print *, '_NDRV1'
         print *, '====== '
         print *, ' '
         print *, ' Size of the matrix is ', n
         print *, ' The number of Ritz values requested is ', nev
         print *, ' The number of Arnoldi vectors generated (NCV) is ', ncv
         print *, ' What portion of the spectrum: ', which
         print *, ' The number of converged Ritz values is ', iparam(5)
         print *, ' The number of Implicit Arnoldi update iterations taken is ', iparam(3)
         print *, ' The number of OP*x is ', iparam(9)
         print *, ' The convergence criterion is ', tol
         print *, ' '

      end if



    contains

        subroutine av ( v, w, mat)

            complex(R_P), dimension(:), intent(in) :: v
            complex(R_P), dimension(:), intent(out) :: w
            complex(R_P), dimension(:,:), intent(in) :: mat

            w = matmul(mat,v)

        end subroutine av

    end subroutine zndrv1




end module lapack_wrap
