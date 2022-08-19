!> Finite Volume Method(FVM) to calculate one-dimension steady convection diffusion case 
!!
!! @author Yan Liu
!! @date August,19,2022

module FVMFor1DimConvectionDiffusion
    implicit none
    
contains
    function FVMCalc1DimConvectionDiffusion(para) result(mk)
        ! Function for one-dimensional steady convection-diffusion case
            use matrixFunctions
            use customTypes
            implicit none

            integer :: mk, i

            ! declare boundary condition variable and so on (/ w, west / e, east / n, north / s, south /)
            type(FVM_1DimType) :: para
            
            ! declare convective mass flux per unit area, F and diffusion conductance at cell faces, D 
            real :: F, D
            
            ! declare euqation nodes coefficient
            real, dimension(:, :), allocatable :: nodes_coeff_we
            real, dimension(:), allocatable :: nodes_t, nodes_Su
            real :: node_aW, node_aP, node_aE, node_Sp

            F = para%density * para%u
            D = para%diffusion_coefficient / (para%distance_we / para%grid_num_we)
        
            allocate(nodes_coeff_we(para%grid_num_we, para%grid_num_we))
            allocate(nodes_t(para%grid_num_we))
            allocate(nodes_Su(para%grid_num_we))

            ! set nodes data to zero
            nodes_coeff_we = 0.0
            nodes_t = 0.0
            nodes_Su = 0.0
        
            ! calculate the nodes coefficient
            ! the first node
            node_aW = 0
            node_aE = D - F / 2
            node_Sp = -(2 * D + F)
            node_aP = node_aW + node_aE - node_Sp
            nodes_coeff_we(1, 1) = node_aP
            nodes_coeff_we(1, 2) = -node_aE
            nodes_Su(1) = (2 * D + F) * para%t_boundary_w
        
            ! the last node
            node_aW = D + F / 2
            node_aE = 0
            node_Sp = -(2 * D - F)
            node_aP = node_aW + node_aE - node_Sp
            nodes_coeff_we(para%grid_num_we, para%grid_num_we) = node_aP
            nodes_coeff_we(para%grid_num_we, para%grid_num_we - 1) = -node_aW
            nodes_Su(para%grid_num_we) = (2 * D - F) * para%t_boundary_e
        
            ! other nodes
            do i = 2, para%grid_num_we - 1
                node_aW = D + F / 2
                node_aE = D - F / 2
                node_aP = node_aW + node_aE
                nodes_coeff_we(i, i - 1) = -node_aW
                nodes_coeff_we(i, i) = node_aP
                nodes_coeff_we(i, i + 1) = -node_aE 
            end do
            
            ! calculate t
            nodes_t = matmul(inv(nodes_coeff_we), nodes_Su)

            ! print result to screen
            ! print *, 'nodes_coeff_we'
            ! do i = 1, ubound(nodes_coeff_we, 1)
            !     print *, nodes_coeff_we(i, :)
            ! end do
            ! print *, 'nodes_Su'
            ! print *, nodes_Su
            print *, 'nodes_t'
            print *, nodes_t
            
            ! deallocate data
            deallocate(nodes_coeff_we)
            deallocate(nodes_Su)
            deallocate(nodes_t)
            
            mk = 1
        end function FVMCalc1DimConvectionDiffusion
end module FVMFor1DimConvectionDiffusion