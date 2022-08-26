module HybridDifferencingCalc3DimConvectionDiffusion
    implicit none
    
contains
function HybirdDiffCalc3DimConvectionDiffusion(para) result(mk)
    ! Function for one-dimensional steady convection-diffusion case
        use matrixFunctions
        use customTypes
        implicit none

        integer :: mk, i

        ! declare 3dimension boundary condition variable, velocity, distance...
        type(FVM_3DimType) :: para
        
        ! declare convective mass flux per unit area, F and diffusion conductance at cell faces, D
        ! we : west - east, sn : south -north, bt : buttom - top 
        real :: F_we, F_sn, F_bt, D_we, D_sn, D_bt
        
        ! declare euqation nodes coefficient
        real, dimension(:, :), allocatable :: nodes_coeff_we
        real, dimension(:), allocatable :: nodes_t, nodes_Su
        real :: node_aW, node_aP, node_aE, node_Sp

        ! Peclet number
        real :: pew

        ! F = para%density * para%u
        ! D = para%diffusion_coefficient / (para%distance_we / para%grid_num_we)

        ! pew = F / D

        ! allocate(nodes_coeff_we(para%grid_num_we, para%grid_num_we))
        ! allocate(nodes_t(para%grid_num_we))
        ! allocate(nodes_Su(para%grid_num_we))

        ! ! set nodes data to zero
        ! nodes_coeff_we = 0.0
        ! nodes_t = 0.0
        ! nodes_Su = 0.0
    
        ! ! calculate the nodes coefficient
        ! ! the first node
        ! node_aW = 0
        ! node_aE = max(-F, D - F / 2, 0.0)
        ! if (abs(pew) < 2) then
        !     node_Sp = -(2 * D + F)
        ! else
        !     node_Sp = -(2 * D + max(F, 0.0))
        ! end if
        ! node_aP = node_aW + node_aE - node_Sp
        ! nodes_coeff_we(1, 1) = node_aP
        ! nodes_coeff_we(1, 2) = -node_aE
        ! nodes_Su(1) = (-node_Sp) * para%t_boundary_w
    
        ! ! the last node
        ! node_aW = max(F, D + F / 2, 0.0)
        ! node_aE = 0
        ! if (abs(pew) < 2) then
        !     node_Sp = -(2 * D - F)    
        ! else
        !     node_Sp = -(2 * D + max(-F, 0.0))
        ! end if  
        ! node_aP = node_aW + node_aE - node_Sp
        ! nodes_coeff_we(para%grid_num_we, para%grid_num_we) = node_aP
        ! nodes_coeff_we(para%grid_num_we, para%grid_num_we - 1) = -node_aW
        ! nodes_Su(para%grid_num_we) = (-node_Sp) * para%t_boundary_e
    
        ! ! other nodes
        ! do i = 2, para%grid_num_we - 1
        !     node_aW = max(F, D + F / 2, 0.0)
        !     node_aE = max(-F, D - F / 2, 0.0)
        !     node_aP = node_aW + node_aE
        !     nodes_coeff_we(i, i - 1) = -node_aW
        !     nodes_coeff_we(i, i) = node_aP
        !     nodes_coeff_we(i, i + 1) = -node_aE 
        ! end do
        
        ! ! calculate t
        ! nodes_t = matmul(inv(nodes_coeff_we), nodes_Su)

        ! ! print result to screen
        ! print *, 'nodes_t'
        ! print *, nodes_t
        
        ! ! deallocate data
        ! deallocate(nodes_coeff_we)
        ! deallocate(nodes_Su)
        ! deallocate(nodes_t)
        
        mk = 1
    end function HybirdDiffCalc3DimConvectionDiffusion
    
end module HybridDifferencingCalc3DimConvectionDiffusion