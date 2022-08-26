module HybridDifferencingCalc2DimConvectionDiffusion
    implicit none
    
contains
function HybirdDiffCalc2DimConvectionDiffusion(para) result(mk)
    ! Function for two-dimension steady convection-diffusion case
        use matrixFunctions
        use customTypes
        implicit none

        integer :: mk, i, j

        ! declare 2-dimension boundary condition variable, velocity, distance...
        type(FVM_2DimType) :: para
        
        ! declare convective mass flux per unit area, F and diffusion conductance at cell faces, D
        ! we : west - east, sn : south -north, bt : buttom - top 
        real :: F_we, F_sn, D_we, D_sn
        
        ! declare euqation nodes coefficient
        real, dimension(:, :), allocatable :: nodes_coeff

        ! declare result and source value 
        real, dimension(:), allocatable :: nodes_t, nodes_Su

        ! declare coefficient for every single node and it's adjacent nodes  
        real :: node_aW, node_aE, node_aS, node_aN, node_aP, node_Sp_we, node_Sp_sn

        ! declare Peclet number
        real :: pew_we, pew_sn

        ! declare node nums
        integer :: grid_sum

        F_we = para%density * para%u
        F_sn = para%density * para%v

        D_we = para%diffusion_coefficient_we / (para%distance_we / para%grid_num_we)
        D_sn = para%diffusion_coefficient_sn / (para%distance_sn / para%grid_num_sn)

        pew_we = F_we / D_we
        pew_sn = F_sn / D_sn

        grid_sum = para%grid_num_we * para%grid_num_sn

        allocate(nodes_coeff(grid_sum, grid_sum))
        allocate(nodes_t(grid_sum))
        allocate(nodes_Su(grid_sum))

        ! set nodes data to zero
        nodes_coeff = 0.0
        nodes_t = 0.0
        nodes_Su = 0.0
    
        ! calculate the nodes coefficient
        ! the first node (bottom left)
        node_aW = 0
        node_aE = max(-F_we, D_we - F_we / 2.0, 0.0)
        node_aS = 0
        node_aN = max(-F_sn, D_sn - F_sn / 2.0, 0.0)
        node_Sp_we = CalcFirstNodeSoure(pew_we, D_we, F_we)
        node_Sp_sn = CalcFirstNodeSoure(pew_sn, D_sn, F_sn)
        node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
        
        nodes_coeff(1, 1) = node_aP
        nodes_coeff(1, 2) = -node_aE
        nodes_coeff(1, 1 + para%grid_num_we) = -node_aN
        nodes_Su(1) = (-node_Sp_we) * para%t_boundary_w + (-node_Sp_sn) * para%t_boundary_s
    
        ! the south line nodes (except bottom left and bottom right)
        do i = 2, para%grid_num_we - 1
            node_aW = max(F_we, D_we + F_we / 2.0, 0.0)
            node_aE = max(-F_we, D_we - F_we / 2.0, 0.0)
            node_aS = 0
            node_aN = max(-F_sn, D_sn - F_sn / 2.0, 0.0)
            node_Sp_we = 0
            node_Sp_sn = CalcFirstNodeSoure(pew_sn, D_sn, F_sn)
            node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
            nodes_coeff(i, i - 1) = -node_aW
            nodes_coeff(i, i) = node_aP
            nodes_coeff(i, i + 1) = -node_aE
            nodes_coeff(i, i + para%grid_num_we) = -node_aN
            nodes_Su(i) = (-node_Sp_sn) * para%t_boundary_s
        end do
        
        ! the bottom right node
        node_aW = max(F_we, D_we + F_we / 2.0, 0.0)
        node_aE = 0
        node_aS = 0
        node_aN = max(-F_sn, D_sn - F_sn / 2.0, 0.0)
        node_Sp_we = CalcLastNodeSoure(pew_we, D_we, F_we)
        node_Sp_sn = CalcFirstNodeSoure(pew_sn, D_sn, F_sn)
        node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
        
        nodes_coeff(para%grid_num_we, para%grid_num_we) = node_aP
        nodes_coeff(para%grid_num_we, para%grid_num_we - 1) = -node_aW
        nodes_coeff(para%grid_num_we, para%grid_num_we + para%grid_num_we) = -node_aN
        nodes_Su(para%grid_num_we) = (-node_Sp_we) * para%t_boundary_e + (-node_Sp_sn) * para%t_boundary_s
        
        ! the east line nodes (except top right and bottom right)
        do i = 2, para%grid_num_sn - 1
            node_aW = max(F_we, D_we + F_we / 2.0, 0.0)
            node_aE = 0
            node_aS = max(F_sn, D_sn + F_sn / 2.0, 0.0)
            node_aN = max(-F_sn, D_sn - F_sn / 2.0, 0.0)
            node_Sp_we = CalcLastNodeSoure(pew_we, D_we, F_we)
            node_Sp_sn = 0
            node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
            nodes_coeff(i * para%grid_num_we, i * para%grid_num_we - 1) = -node_aW
            nodes_coeff(i * para%grid_num_we, i * para%grid_num_we) = node_aP
            nodes_coeff(i * para%grid_num_we, (i - 1) * para%grid_num_we) = -node_aS
            nodes_coeff(i * para%grid_num_we, (i + 1) * para%grid_num_we) = -node_aN
            nodes_Su(i * para%grid_num_we) = (-node_Sp_we) * para%t_boundary_e
        end do

        ! the last node (top right)
        node_aW = max(F_we, D_we + F_we / 2.0, 0.0)
        node_aE = 0
        node_aS = max(F_sn, D_sn + F_sn / 2.0, 0.0)
        node_aN = 0
        node_Sp_we = CalcLastNodeSoure(pew_we, D_we, F_we)
        node_Sp_sn = CalcLastNodeSoure(pew_sn, D_sn, F_sn)
        node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
        
        nodes_coeff(grid_sum, grid_sum) = node_aP
        nodes_coeff(grid_sum, grid_sum - 1) = -node_aW
        nodes_coeff(grid_sum, grid_sum - para%grid_num_we) = -node_aS
        nodes_Su(grid_sum) = (-node_Sp_we) * para%t_boundary_e + (-node_Sp_sn) * para%t_boundary_n
    
        ! the north line nodes (except top left and top right)
        do i = 2, para%grid_num_we - 1
            node_aW = max(F_we, D_we + F_we / 2.0, 0.0)
            node_aE = max(-F_we, D_we - F_we / 2.0, 0.0)
            node_aS = max(F_sn, D_sn + F_sn / 2.0, 0.0)
            node_aN = 0
            node_Sp_we = 0
            node_Sp_sn = CalcLastNodeSoure(pew_sn, D_sn, F_sn)
            node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
            nodes_coeff(grid_sum - para%grid_num_we + i, grid_sum - para%grid_num_we + i - 1) = -node_aW
            nodes_coeff(grid_sum - para%grid_num_we + i, grid_sum - para%grid_num_we + i) = node_aP
            nodes_coeff(grid_sum - para%grid_num_we + i, grid_sum - para%grid_num_we + i + 1) = -node_aE
            nodes_coeff(grid_sum - para%grid_num_we + i, grid_sum - para%grid_num_we + i - para%grid_num_we) = -node_aS
            nodes_Su(grid_sum - para%grid_num_we + i) = (-node_Sp_sn) * para%t_boundary_n
        end do

        ! the top left node
        node_aW = 0
        node_aE = max(-F_we, D_we - F_we / 2.0, 0.0)
        node_aS = max(F_sn, D_sn + F_sn / 2.0, 0.0)
        node_aN = 0
        node_Sp_we = CalcFirstNodeSoure(pew_we, D_we, F_we)
        node_Sp_sn = CalcLastNodeSoure(pew_sn, D_sn, F_sn)
        node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
        
        nodes_coeff(grid_sum - para%grid_num_we + 1, grid_sum - para%grid_num_we + 1) = node_aP
        nodes_coeff(grid_sum - para%grid_num_we + 1, grid_sum - para%grid_num_we + 1 + 1) = -node_aE
        nodes_coeff(grid_sum - para%grid_num_we + 1, grid_sum - para%grid_num_we + 1 - para%grid_num_we) = -node_aS
        nodes_Su(grid_sum - para%grid_num_we + 1) = (-node_Sp_we) * para%t_boundary_w + (-node_Sp_sn) * para%t_boundary_n
        
        ! the west line nodes (except bottom left and top left)
        do i = 2, para%grid_num_sn - 1
            node_aW = 0
            node_aE = max(-F_we, D_we - F_we / 2.0, 0.0)
            node_aS = max(F_sn, D_sn + F_sn / 2.0, 0.0)
            node_aN = max(-F_sn, D_sn - F_sn / 2.0, 0.0)
            node_Sp_we = CalcFirstNodeSoure(pew_we, D_we, F_we)
            node_Sp_sn = 0
            node_aP = node_aW + node_aE + node_aS + node_aN - node_Sp_we - node_Sp_sn
            nodes_coeff(1 + (i - 1) * para%grid_num_we, 1 + (i - 1) * para%grid_num_we + 1) = -node_aE
            nodes_coeff(1 + (i - 1) * para%grid_num_we, 1 + (i - 1) * para%grid_num_we) = node_aP
            nodes_coeff(1 + (i - 1) * para%grid_num_we, 1 + (i - 1) * para%grid_num_we - para%grid_num_we) = -node_aS
            nodes_coeff(1 + (i - 1) * para%grid_num_we, 1 + (i - 1) * para%grid_num_we + para%grid_num_we) = -node_aN
            nodes_Su(1 + (i - 1) * para%grid_num_we) = (-node_Sp_we) * para%t_boundary_w
        end do

        ! boundary nodes
    
        ! other nodes
        do i = 2, para%grid_num_sn - 1
            do j = 2, para%grid_num_we - 1
                node_aW = max(F_we, D_we + F_we / 2.0, 0.0)
                node_aE = max(-F_we, D_we - F_we / 2.0, 0.0)
                node_aS = max(F_sn, D_sn + F_sn / 2.0, 0.0)
                node_aN = max(-F_sn, D_sn - F_sn / 2.0, 0.0)
                node_aP = node_aW + node_aE + node_aS + node_aN
                nodes_coeff(j + (i - 1) * para%grid_num_we, j + (i - 1) * para%grid_num_we - 1) = -node_aW
                nodes_coeff(j + (i - 1) * para%grid_num_we, j + (i - 1) * para%grid_num_we) = node_aP
                nodes_coeff(j + (i - 1) * para%grid_num_we, j + (i - 1) * para%grid_num_we + 1) = -node_aE
                nodes_coeff(j + (i - 1) * para%grid_num_we, j + (i - 1) * para%grid_num_we - para%grid_num_we) = -node_aS
                nodes_coeff(j + (i - 1) * para%grid_num_we, j + (i - 1) * para%grid_num_we + para%grid_num_we) = -node_aN
            end do
        end do
        
        ! print *, ubound(nodes_coeff, 1)

        ! print nodes coefficents
        ! do i = 1, ubound(nodes_coeff, 1)
        !     print *, nodes_coeff(i,:)
        ! end do  

        ! calculate t
        nodes_t = matmul(inv(nodes_coeff), nodes_Su)

        ! print result to screen
        do i = 1, para%grid_num_sn
            print *, nodes_t(grid_sum - i * para%grid_num_we + 1: grid_sum - (i - 1) * para%grid_num_we)
        end do
        
        ! print *, nodes_Su

        print *, 'F:', F_we , 'D:', D_we, pew_we
        ! print *, F_sn , D_sn, pew_sn

        ! deallocate data
        deallocate(nodes_coeff)
        deallocate(nodes_Su)
        deallocate(nodes_t)
        
        mk = 1
    contains
        function CalcFirstNodeSoure(pew, D, F) result(Sp)
            implicit none
            real :: pew, D, F, Sp

            if (abs(pew) < 2) then
                Sp = -(2 * D + F)    
            else
                Sp = -(2 * D + max(F, 0.0))
            end if 
        end function CalcFirstNodeSoure

        function CalcLastNodeSoure(pew, D, F) result(Sp)
            implicit none
            real :: pew, D, F, Sp

            if (abs(pew) < 2) then
                Sp = -(2 * D - F)    
            else
                Sp = -(2 * D + max(-F, 0.0))
            end if 
        end function CalcLastNodeSoure
end function HybirdDiffCalc2DimConvectionDiffusion
    
end module HybridDifferencingCalc2DimConvectionDiffusion