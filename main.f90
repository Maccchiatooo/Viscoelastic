program main
  
  use lbm


 
  

  call DataInput
  call MPI_INIT (ierr)
  call MPI_COMM_RANK (MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE (MPI_COMM_WORLD, num_procs, ierr)
  call MPISetup

  
  
  call AllocateArrays
  
  call Initilization
  
  call test 
  iter=0
  call MPIoutput

  do iter=1,max_step
    call Collision

    call BoundaryCondition
    
    call Propagation
    
    call PostProcessing
    
    if(mod(iter,interv) == 0)then
    call MPIoutput
    endif

  enddo

  call DeAllocateArrays
  call MPI_FINALIZE (ierr)

endprogram main
