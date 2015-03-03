module abstract_test_class
    use gnufor2
    use IR_Precision
    implicit none
    private

    type, abstract, public :: test_case
    private
    contains
        procedure(abs_init), deferred :: init
        procedure(abs_run_test), deferred :: run_test
        procedure(abs_evaluate), deferred :: evaluate
        procedure, non_overridable :: test
    end type test_case

    abstract interface
        subroutine abs_init(self,test_size)
            import :: test_case
            class(test_case), intent(inout) :: self
            integer, intent(in) :: test_size
        end subroutine abs_init
    end interface

    abstract interface
        subroutine abs_run_test(self)
            import :: test_case
            class(test_case), intent(inout) :: self
        end subroutine abs_run_test
    end interface

    abstract interface
        function abs_evaluate(self) result(test_res)
            import :: test_case, R_P
            class(test_case), intent(in) :: self
            real(R_P) :: test_res
        end function abs_evaluate
    end interface

contains

    subroutine test(self,max_duration,avg,sze)

        class(test_case), intent(inout) :: self
        real(R_P), intent(in) :: max_duration
        integer, intent(in) :: avg
        integer, dimension(2), intent(in) :: sze

        integer :: i, j, test_size
        real(R_P) :: start, stop, accuracy, duration
        real(R_P), dimension(100) :: accu_array, time_array, sze_array

        test_size = sze(1)
        j = 0

        do while((duration < max_duration) .and. (j < 100))
            accuracy = 0.0_R_P
            duration = 0.0_R_P
            j = j + 1
            do i = 1,avg

                !Initialize Test
                call self%init(test_size)

                !Perform Test and Measure Time
                call cpu_time(start)
                call self%run_test()
                call cpu_time(stop)


                !Evaluate Accuracy and Speed of Test
                accuracy = accuracy + self%evaluate()
                duration = duration + stop - start

            end do

            accu_array(j) = accuracy / avg
            time_array(j) = duration / avg
            sze_array(j) = real(test_size,R_P)
            test_size = test_size + sze(2)

        end do

        call plot(sze_array(1:j),accu_array(1:j),terminal='ps',filename='accuracy.ps')
        call plot(sze_array(1:j),time_array(1:j),terminal='ps',filename='time.ps')

    end subroutine

end module abstract_test_class
