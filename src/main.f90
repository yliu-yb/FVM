!> WRF SCM Extension Program
!!
!! @author Yan Liu
!! @date August,18,2022
!! 

program SCMExtension
    use customTypes
    use FVMFor1DimConvectionDiffusion
    implicit none

    type(FVM_1DimType) :: para
    integer :: mk

    print *, "case1"
    para = FVM_1DimType(1.0, 0.0, 1.0, 0.1, 1.0, 0.1, 5)
    mk = FVMCalc1DimConvectionDiffusion(para)

    print *, "case2"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 5)
    mk = FVMCalc1DimConvectionDiffusion(para)

    print *, "case3"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 20)
    mk = FVMCalc1DimConvectionDiffusion(para)

end program SCMExtension
    
    

