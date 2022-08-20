!> WRF SCM Extension Program
!!
!! @author Yan Liu
!! @date August,18,2022
!! 

program SCMExtension
    use customTypes
    use CentralDifferencingCalc1DimConvectionDiffusion
    use UpwindDifferencingCalc1DimConvectionDiffusion
    use HybridDifferencingCalc1DimConvectionDiffusion
    implicit none

    type(FVM_1DimType) :: para
    integer :: mk

    print *, "central differencing scheme"
    ! print *, "case1"
    ! para = FVM_1DimType(1.0, 0.0, 1.0, 0.1, 1.0, 0.1, 5)
    ! mk = CentralDiffCalc1DimConvectionDiffusion(para)

    ! print *, "case2"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 25)
    mk = CentralDiffCalc1DimConvectionDiffusion(para)

    ! print *, "case3"
    ! para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 20)
    ! mk = CentralDiffCalc1DimConvectionDiffusion(para)

    ! print *, "upwind differencing scheme"
    ! print *, "case1 - positive u"
    ! para = FVM_1DimType(1.0, 0.0, 1.0, 0.1, 1.0, 0.1, 5)
    ! mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    ! print *, "case1 - negtive u"
    ! para = FVM_1DimType(0.0, 1.0, 1.0, -0.1, 1.0, 0.1, 5)
    ! mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    ! print *, "case2 - positive u"
    ! para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 5)
    ! mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    ! print *, "case2 - negtive u"
    ! para = FVM_1DimType(0.0, 1.0, 1.0, -2.5, 1.0, 0.1, 5)
    ! mk = UpwindDiffCalc1DimConvectionDiffusion(para)

    print *, "hybrid differencing scheme"
    print *, "case2 - positive u, node num: 5"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 5)
    mk = HybirdDiffCalc1DimConvectionDiffusion(para)

    print *, "case2 - negtive u, node num: 5"
    para = FVM_1DimType(0.0, 1.0, 1.0, -2.5, 1.0, 0.1, 5)
    mk = HybirdDiffCalc1DimConvectionDiffusion(para)

    print *, "case2 - positive u, node num: 25"
    para = FVM_1DimType(1.0, 0.0, 1.0, 2.5, 1.0, 0.1, 25)
    mk = HybirdDiffCalc1DimConvectionDiffusion(para)

end program SCMExtension
    
    

