
module customTypes
    implicit none
    ! boundary condition variable (/ w, west / e, east / n, north / s, south / b, bottom / t, top/)
    type FVM_1DimType
        real :: t_boundary_w
        real :: t_boundary_e
        real :: distance_we
        real :: u
        real :: density
        real :: diffusion_coefficient
        integer :: grid_num_we
    end type FVM_1DimType

    type FVM_2DimType
        real :: t_boundary_w
        real :: t_boundary_e
        real :: distance_we
        real :: u
        real :: diffusion_coefficient_we
        integer :: grid_num_we
        real :: t_boundary_s
        real :: t_boundary_n
        real :: distance_sn
        real :: v
        real :: diffusion_coefficient_sn
        integer :: grid_num_sn
        real :: density
    end type FVM_2DimType

    type FVM_3DimType
        real :: t_boundary_w
        real :: t_boundary_e
        real :: t_boundary_s
        real :: t_boundary_n
        real :: t_boundary_b
        real :: t_boundary_t
        real :: distance_we
        real :: distance_sn
        real :: distance_bt
        real :: u
        real :: v
        real :: w
        real :: density
        real :: diffusion_coefficient_we
        real :: diffusion_coefficient_sn
        real :: diffusion_coefficient_tp
        integer :: grid_num_we
        integer :: grid_num_sn
        integer :: grid_num_bt
    end type FVM_3DimType

contains
    
end module customTypes