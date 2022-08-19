
module customTypes
    implicit none
    type FVM_1DimType
        real :: t_boundary_w
        real :: t_boundary_e
        real :: distance_we
        real :: u
        real :: density
        real :: diffusion_coefficient
        integer :: grid_num_we
    end type FVM_1DimType
contains
    
end module customTypes