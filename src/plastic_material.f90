module plastic_material

  !> plastic_material module defines initial plastic strains and strain rates  on a block
  
  use common, only : wp
  use datatypes, only : block_plastic,block_grid_t,block_indices
  implicit none
  
contains
  
  
  !> initialize plastic strains and strain-rates
  
  subroutine init_plastic_material(P, G, I, problem, mu_beta_eta)
    
    use mpi3dcomm
    
    implicit none
    
    type(block_plastic),intent(out) :: P
    type(block_grid_t),intent(in) :: G
    type(block_indices),intent(in) :: I
    real(kind = wp), intent(in) :: mu_beta_eta(:)
    character(*),intent(in) :: problem
    integer :: l,j,k
    
    
    call allocate_array_body(P%P,G%C,2, ghost_nodes= .false.)
    
    select case(problem)
       
    case default
       
       P%P(:,:,:,1) = 0.0_wp                  ! plastic strain 
       P%P(:,:,:,2) = 0.0_wp                  ! plastic strain rate
       
       
    case('rate-weakening')
       
       P%mu_beta_eta(1) = mu_beta_eta(1)                                                                                        
       P%mu_beta_eta(2) = mu_beta_eta(2) 
       P%mu_beta_eta(3) = mu_beta_eta(3) 
       P%P(:,:,:,1) = 0.0_wp                  ! plastic strain                                                                 
       P%P(:,:,:,2) = 0.0_wp                  ! plastic strain rate
       
    case('TPV31')

       !P%mu_beta_eta(1) = mu_beta_eta(1)
       !P%mu_beta_eta(2) = mu_beta_eta(2)
       !P%mu_beta_eta(3) = mu_beta_eta(3)

       P%P(:,:,:,1) = 0.0_wp                  ! plastic strain       
       P%P(:,:,:,2) = 0.0_wp                  ! plastic strain rate  
       
    end select
    
    
    
  end subroutine init_plastic_material

end module plastic_material

