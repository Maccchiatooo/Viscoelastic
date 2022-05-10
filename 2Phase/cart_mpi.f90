module cart_mpi
  use mpi
  implicit none
  
  integer,parameter:: dim=2,nq=9,ghost=2,num_neighbors=8
  integer nx,ny
  
  integer ierr, num_procs, rank, fh

  integer CART_COMM
  integer,dimension(dim)::cart_coords,cart_num_procs
  logical,dimension(dim)::cart_periods
  integer,dimension(dim)::global_length,local_length,local_start,local_end
  
  integer,dimension(0:num_neighbors)::neighbor_procs  
  integer,dimension(0:num_neighbors)::INT_RECV_TYPE,F_RECV_TYPE,DPR_RECV_TYPE,LOG_RECV_TYPE,U_RECV_TYPE
  integer,dimension(0:num_neighbors)::INT_SEND_TYPE,F_SEND_TYPE,DPR_SEND_TYPE,LOG_SEND_TYPE,U_SEND_TYPE

  integer,dimension(0:num_neighbors)::neighbor_procsv  
  integer,dimension(0:num_neighbors)::V_RECV_TYPE
  integer,dimension(0:num_neighbors)::V_SEND_TYPE

contains
  subroutine MPISetup
     
    call SetCartesian
    call SetLocalLength
    call SetNeighbors

  endsubroutine MPISetup
  
!------------------------------------------
  subroutine SetCartesian   
    double precision max_e, config_e
    integer px,py
    logical reorder_flag
    integer i
    
    global_length(1)=nx+1
    global_length(2)=ny+1
    
    max_e=0.d0
    do i=1,num_procs
       if(mod(num_procs,i).eq.0)then
          px=i
          py=num_procs/i
          config_e=1.d0/(2*(global_length(2)*(px-1)/py&
               &+global_length(1)*(py-1)/px))
          if(config_e.ge.max_e)then
             max_e=config_e
             cart_num_procs(1)=px
             cart_num_procs(2)=py
          endif
       endif
    enddo

    do i=1,dim
       cart_periods(i)=.true.
    enddo

    reorder_flag=.true.
    call MPI_CART_CREATE(MPI_COMM_WORLD,dim,cart_num_procs,cart_periods,reorder_flag,CART_COMM,ierr)
    call MPI_COMM_RANK(CART_COMM,rank,ierr)
    call MPI_CART_COORDS(CART_COMM,rank,dim,cart_coords,ierr)
    call MPI_BARRIER(CART_COMM,ierr)    
  endsubroutine SetCartesian
  
!-------------------------------------------------------------  
  subroutine SetLocalLength
    integer i,type1
    integer,dimension(dim,num_procs)::length_public
    integer,dimension(dim)::this_cart_coords
    integer,dimension(:,:,:),allocatable::local_length_db
    do i=1,dim
       local_length(i)=global_length(i)/cart_num_procs(i)+(1-sign(1,(cart_coords(i)-mod(global_length(i),cart_num_procs(i)))))/2
       if(local_length(i)<ghost)then
          print*, 'Dimension = ',i,' in Cartesian proc id =',rank,' smaller than GhostCells'
          call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
       endif
    enddo

    allocate(local_length_db(dim,0:cart_num_procs(1)-1,0:cart_num_procs(2)-1))
    
    call MPI_TYPE_VECTOR(dim,1,1,MPI_INTEGER,type1,ierr)
    call MPI_TYPE_COMMIT(type1,ierr)
    call MPI_ALLGATHER(local_length(:),1,type1,length_public,1,type1,CART_COMM,ierr)
    call MPI_TYPE_FREE(type1,ierr)
    do i=1,num_procs
       call MPI_CART_COORDS(CART_COMM,i-1,dim,this_cart_coords,ierr)
       local_length_db(:,this_cart_coords(1),this_cart_coords(2))=length_public(:,i)
    enddo
    
    local_start(1)=sum(local_length_db(1,0:cart_coords(1)-1,cart_coords(2)))
    local_start(2)=sum(local_length_db(2,cart_coords(1),0:cart_coords(2)-1))

    local_end(1)=sum(local_length_db(1,0:cart_coords(1),cart_coords(2)))-1
    local_end(2)=sum(local_length_db(2,cart_coords(1),0:cart_coords(2)))-1

    deallocate(local_length_db)
    
  endsubroutine SetLocalLength

 !--------------------------------------------
  subroutine SetNeighbors
    integer i,j,ind, neighbor_local_rank,neighbor_pos_p
    integer,dimension(dim)::neighbor_pos
    integer,dimension(dim)::neighbor_cart_coords
    integer,dimension(dim+1)::fsizes,fsubsize,frecv_start,fsend_start,usizes,usubsize,urecv_start,usend_start
    integer,dimension(dim+1)::vsizes,vsubsize,vrecv_start,vsend_start
    integer,dimension(dim)::sizes,subsize,recv_start,send_start
    
    do i=-1,1
       do j=-1,1
          neighbor_pos(1)=i
          neighbor_pos(2)=j
          neighbor_local_rank=(i+1)*3+(j+1)

          do ind=1,dim
             neighbor_cart_coords(ind)=cart_coords(ind)+neighbor_pos(ind)
          enddo

          call MPI_CART_RANK(CART_COMM,neighbor_cart_coords,neighbor_procs(neighbor_local_rank),ierr)

          fsizes(1)=nq
          fsubsize(1)=nq
          frecv_start(1)=0
          fsend_start(1)=0

          vsizes(1)=5
          vsubsize(1)=5
          vrecv_start(1)=0
          vsend_start(1)=0
          
          do ind=1,dim
             neighbor_pos_p=neighbor_pos(ind)
             
             sizes(ind)=local_length(ind)+2*ghost             
             subsize(ind)=local_length(ind)*(1-abs(neighbor_pos_p))+ghost*abs(neighbor_pos_p)
             recv_start(ind)=neighbor_pos_p*(1-neighbor_pos_p)/2*ghost+neighbor_pos_p*(1+neighbor_pos_p)/2*local_length(ind)+ghost
             send_start(ind)=neighbor_pos_p*(1+neighbor_pos_p)/2*(local_length(ind)-ghost)+ghost

             fsizes(ind+1)=sizes(ind)
             fsubsize(ind+1)=subsize(ind)
             frecv_start(ind+1)=recv_start(ind)
             fsend_start(ind+1)=send_start(ind)

             vsizes(ind+1)=sizes(ind)
             vsubsize(ind+1)=subsize(ind)
             vrecv_start(ind+1)=recv_start(ind)
             vsend_start(ind+1)=send_start(ind)

             usizes(ind)=sizes(ind)
             usubsize(ind)=subsize(ind)
             urecv_start(ind)=recv_start(ind)
             usend_start(ind)=send_start(ind)
          enddo
          usizes(dim+1)=dim
          usubsize(dim+1)=dim
          urecv_start(dim+1)=0!1
          usend_start(dim+1)=0!1          

          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,recv_start,MPI_ORDER_FORTRAN,&
               &MPI_INTEGER,INT_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(INT_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,send_start,MPI_ORDER_FORTRAN,&
               &MPI_INTEGER,INT_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(INT_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim+1,vsizes,vsubsize,vrecv_start,MPI_ORDER_FORTRAN,&
          &MPI_DOUBLE_PRECISION,V_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(V_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,vsizes,vsubsize,vsend_start,MPI_ORDER_FORTRAN,&
          &MPI_DOUBLE_PRECISION,V_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(V_SEND_TYPE(neighbor_local_rank),ierr)
          
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,fsizes,fsubsize,frecv_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,F_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(F_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,fsizes,fsubsize,fsend_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,F_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(F_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,recv_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,DPR_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(DPR_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,send_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,DPR_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(DPR_SEND_TYPE(neighbor_local_rank),ierr)

          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,recv_start,MPI_ORDER_FORTRAN,&
               &MPI_LOGICAL,LOG_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(LOG_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim,sizes,subsize,send_start,MPI_ORDER_FORTRAN,&
               &MPI_LOGICAL,LOG_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(LOG_SEND_TYPE(neighbor_local_rank),ierr)
          
          
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,usizes,usubsize,urecv_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,U_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(U_RECV_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_CREATE_SUBARRAY(dim+1,usizes,usubsize,usend_start,MPI_ORDER_FORTRAN,&
               &MPI_DOUBLE_PRECISION,U_SEND_TYPE(neighbor_local_rank),ierr)
          call MPI_TYPE_COMMIT(U_SEND_TYPE(neighbor_local_rank),ierr)


          
       enddo
    enddo
    
    call MPI_BARRIER(MPI_COMM_WORLD,ierr)
    
  endsubroutine SetNeighbors
  
!-------------------------------------------------
subroutine PassF(Array)

   double precision,dimension(0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::array
   integer neighbor_local_rank,neighbor_oppos_rank
   
   integer,dimension(0:num_neighbors,2):: req
   integer,dimension(MPI_STATUS_SIZE,0:num_neighbors,2) :: communication_status

   call MPI_BARRIER(CART_COMM,ierr)

   do neighbor_local_rank=0,num_neighbors

      neighbor_oppos_rank=num_neighbors-neighbor_local_rank

      call MPI_ISEND(Array,1,F_SEND_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,1),ierr)
      call MPI_IRECV(Array,1,F_RECV_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,2),communication_status(:,neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,1),communication_status(:,neighbor_local_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

   enddo
endsubroutine PassF

!-------------------------------------------------
subroutine PassD(Array)

   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::array
   integer neighbor_local_rank,neighbor_oppos_rank
   
   integer,dimension(0:num_neighbors,2):: req
   integer,dimension(MPI_STATUS_SIZE,0:num_neighbors,2) :: communication_status

   call MPI_BARRIER(CART_COMM,ierr)

   do neighbor_local_rank=0,num_neighbors

      neighbor_oppos_rank=num_neighbors-neighbor_local_rank

      call MPI_ISEND(Array,1,DPR_SEND_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,1),ierr)
      call MPI_IRECV(Array,1,DPR_RECV_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,2),communication_status(:,neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,1),communication_status(:,neighbor_local_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

   enddo
 endsubroutine PassD

 !-------------------------------------------------
subroutine PassU(Array)

   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost,dim)::array
   integer neighbor_local_rank,neighbor_oppos_rank
   
   integer,dimension(0:num_neighbors,2):: req
   integer,dimension(MPI_STATUS_SIZE,0:num_neighbors,2) :: communication_status

   call MPI_BARRIER(CART_COMM,ierr)

   do neighbor_local_rank=0,num_neighbors

      neighbor_oppos_rank=num_neighbors-neighbor_local_rank

      call MPI_ISEND(Array,1,U_SEND_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,1),ierr)
      call MPI_IRECV(Array,1,U_RECV_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
           &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,2),communication_status(:,neighbor_local_rank,2),ierr)
      call MPI_WAIT(req(neighbor_local_rank,1),communication_status(:,neighbor_local_rank,1),ierr)

      call MPI_BARRIER(CART_COMM,ierr)

   enddo
 endsubroutine PassU

 subroutine PassV(Array)

     double precision,dimension(0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::array
     integer neighbor_local_rank,neighbor_oppos_rank
     
     integer,dimension(0:num_neighbors,2):: req
     integer,dimension(MPI_STATUS_SIZE,0:num_neighbors,2) :: communication_status
  
     call MPI_BARRIER(CART_COMM,ierr)
  
     do neighbor_local_rank=0,num_neighbors
  
        neighbor_oppos_rank=num_neighbors-neighbor_local_rank
  
        call MPI_ISEND(Array,1,V_SEND_TYPE(neighbor_local_rank),neighbor_procs(neighbor_local_rank),&
             &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,1),ierr)
        call MPI_IRECV(Array,1,V_RECV_TYPE(neighbor_oppos_rank),neighbor_procs(neighbor_oppos_rank),&
             &neighbor_local_rank,CART_COMM,req(neighbor_local_rank,2),ierr)

        call MPI_WAIT(req(neighbor_local_rank,2),communication_status(:,neighbor_local_rank,2),ierr)
        call MPI_WAIT(req(neighbor_local_rank,1),communication_status(:,neighbor_local_rank,1),ierr)
  
        call MPI_BARRIER(CART_COMM,ierr)
  
     enddo
  endsubroutine PassV
 
endmodule cart_mpi
