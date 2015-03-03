program testrun
    use IR_Precision
    use reference_node_test_class
    implicit none

    type(reference_node_test) :: denseNode

    call denseNode%test(60.0_R_P,10,[50,10])

! CORRECTNESS TESTS
!    if(zero_constructor_test()) then
!        print *, "Zero constructor test successful"
!    end if

! PRECISION AND SPEED TESTS
 !   call svd_precision_speed_test()
 !   call multiply_precision_speed_test()
 !   call normalize_test
 !   call gs_search_precision_speed_test()
 !   call transpose_precision_speed_test

! SPEED TESTS

end program testrun
