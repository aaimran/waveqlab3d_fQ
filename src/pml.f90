module pml

  !> pml module allocates PML blocks and initializes auxiliary fieldls and rates inside the PML

  use common, only : wp
  use datatypes, only : block_type, block_grid_t,block_indices, block_pml
  implicit none

contains

  subroutine init_pml(B, G, id, npml)
    
    use mpi3dcomm
    implicit none
    
    integer, intent(in) :: id, npml
    !logical, intent(in) :: pml
    type(block_pml), intent(inout) :: B
    type(block_grid_t),intent(in) :: G

    
    if(B%pml .EQV. .TRUE.) then
       
    
       B%N_pml = npml           ! number of grid points inside the PML

    
       ! left q-boundary
       if (id == 1) then
          B%C = G%C
          B%C%pq = npml
          
          ! allocate fields and rates only if points live inside the PML
          if(G%C%mq .ge. npml+1) return 
          
          if (G%C%pq .ge. npml+1) then 
             allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             
          else
             allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             
          end if
          B%Q(:,:,:,:) = 0.0_wp
          B%DQ(:,:,:,:) = 1.0e40_wp
       end if
       
       ! right q-boundary
       if (id == 2) then
          B%C = G%C
          B%C%mq = G%C%nq-npml+1
          
          ! allocate fields and rates only if points live inside the PML
          if(G%C%pq .le. (G%C%nq-npml)) return
          
          if (G%C%mq .le. G%C%nq-npml+1) then
             allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             
          else
             allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             
          end if
          
          !allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:6))
          !allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:6))
          !allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:6))
          !allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:6))
          
          B%Q(:,:,:,:) = 0.0_wp
          B%DQ(:,:,:,:) = 1.0e40_wp
       end if
       
       ! bottom r-boundary
       if (id == 3) then
          B%C = G%C
          B%C%pr = npml
          
          ! allocate fields and rates only if points live inside the PML
          if(G%C%mr .ge. npml+1) return
          
          if (G%C%pr .ge. npml+1) then
             allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             
          else
             allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             
          end if
          
           B%Q(:,:,:,:) = 0.0_wp
           B%DQ(:,:,:,:) = 1.0e40_wp
        end if
        
        ! top r-boundary
        if (id == 4) then
           B%C = G%C
           B%C%mr = G%C%nr-npml+1
           
           ! allocate fields and rates only if points live inside the PML
           if(G%C%pr .le. (G%C%nr-npml)) return
           
           if (G%C%mr .le. G%C%nr-npml+1) then
              allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
              allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
              
           else
              allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
              allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
              
           end if
           
           B%Q(:,:,:,:) = 0.0_wp
           B%DQ(:,:,:,:) = 1.0e40_wp
        end if
       
        
        ! front s-boundary
        if (id == 5) then
           B%C = G%C
           B%C%ps = npml
           
           ! allocate fields and rates only if points live inside the PML
           if(G%C%ms .ge. npml+1) return
           
           if (G%C%ps .ge. npml+1) then
              allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
              allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
              
           else
              allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
              allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
              
           end if
           
           B%Q(:,:,:,:) = 0.0_wp
           B%DQ(:,:,:,:) = 1.0e40_wp
        end if
        
        ! back s-boundary
       if (id == 6) then
          B%C = G%C
          B%C%ms = G%C%ns-npml+1
          
          
          ! allocate fields and rates only if points live inside the PML
          if(G%C%ps .le. G%C%ns-npml) return
          
          if (G%C%ms .le. G%C%ns-npml+1) then
             allocate(B%Q(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             allocate(B%DQ(B%C%mq:B%C%pq,B%C%mr:B%C%pr,B%C%ms:B%C%ps,1:9))
             
          else
             allocate(B%Q(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
             allocate(B%DQ(G%C%mq:G%C%pq,G%C%mr:G%C%pr,G%C%ms:G%C%ps,1:9))
                          
          end if
          B%Q(:,:,:,:) = 0.0_wp
          B%DQ(:,:,:,:) = 1.0e40_wp
       end if
       ! set fields to zero and rates to a large number
       
    end if
    
  end subroutine init_pml
  
end module pml
