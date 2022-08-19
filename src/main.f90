!> WRF SCM Extension Program
!!
!! @author Yan Liu
!! @date August,18,2022
!! 

program SCMExtension
    use customTypes
    use CentralDifferencingCalc1DimConvectionDiffusion
    use UpwindDifferencingCalc1DimConvectionDiffusion
    implicit none

    type(FVM_1DimType) :: para
    integer :: mk

    print *, "central differencing"
    print *, "case1"
    para = FVM_1DimType(1.0, 0.0, 1.0, 0.1, 1.0, 0.1, 5)
    mk = CentralDiffCalc1DimConvectionDiffusion(para)

    print *, "case2"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 5)
    mk = CentralDiffCalc1DimConvectionDiffusion(para)

    print *, "case3"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 20)
    mk = CentralDiffCalc1DimConvectionDiffusion(para)

    print *, "upwind differencing"
    print *, "case1 - positive u"
    para = FVM_1DimType(1.0, 0.0, 1.0, 0.1, 1.0, 0.1, 5)
    mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    print *, "case1 - negtive u"
    para = FVM_1DimType(0.0, 1.0, 1.0, -0.1, 1.0, 0.1, 5)
    mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    print *, "case2 - positive u"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 5)
    mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    print *, "case2 - negtive u"
    para = FVM_1DimType(0.0, 1.0, 1.0, -2.5, 1.0, 0.1, 5)
    mk = UpwindDiffCalc1DimConvectionDiffusion(para)

end program SCMExtension
    
    

