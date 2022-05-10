module lbm
   use cart_mpi
   implicit none

   !mpioutput
   double precision,dimension(:,:),allocatable:: xx,yy,ppp,ppp2

   !double precision
   double precision v1,c1,v2,c2
   !control parameter
   double precision Oh,Ar,Bo
   double precision eps
   double precision cs2
   !position parameter
   double precision x0,y0,x1,y1,r0,r1
   !density 
   double precision rho1,rho2,rho3
   !relaxation time
   double precision t1,t2,t3,M,tk,ratio
   !surface tension
   double precision sigma12,sigma13,sigma23,sigma1,sigma2,sigma3
   !interface thickness
   double precision delta
   
   double precision forsx1,forsy1,forsx2,forsy2,forsx3,forsy3,forspx,forspy,forsppx,forsppy,forsnx,forsny
   double precision forsa11,forsa12,forsa21,forsa22
   double precision forsax,forsay,fors

   double precision gamma,feq,geq1,geq2,gamma5,aaeq11,aaeq12,aaeq21,aaeq22
   double precision forphi1,forphix1,forphiy1,sqad1,sqad2,sqad3,forphi2,forphix2,forphiy2


  
   double precision psq


   !double precision dimension
   !distribution function
   double precision,dimension(:,:,:),allocatable:: f,g1,g2,fb,gb1,gb2,b11,b12,b21,b22,bb11,bb12,bb21,bb22
   !density
   double precision,dimension(:,:),allocatable::rho,aa11,aa12,aa21,aa22,aaa11,aaa12,aaa21,aaa22
   !pressure
   double precision,dimension(:,:),allocatable:: p,pp
   !viscosity
   double precision,dimension(:,:),allocatable:: nu,mu,dissip
   !velocity
   double precision,dimension(:,:),allocatable:: ua,va,rhou,rhov
   !order parameter
   double precision,dimension(:,:),allocatable:: phi1,phi2,phi3,con
   !relaxation time
   double precision,dimension(:,:),allocatable:: tphi,tau


   !derivative
   !order parameter gradient
   double precision,dimension(:,:),allocatable:: dx1,dy1,dx2,dy2,dx3,dy3
   !rho gradient
   double precision,dimension(:,:),allocatable:: drx,dry
   !pressure gradient
   double precision,dimension(:,:),allocatable:: dpx,dpy,dppx,dppy
   !divergence
   double precision,dimension(:,:),allocatable:: divx1,divx2,divx3,divy1,divy2,divy3,div1,div2,div3,sqd1,sqd2,sqd3
   double precision,dimension(:,:),allocatable:: divx11,divx12,divx21,divx22,divy11,divy12,divy21,divy22
   double precision,dimension(:,:),allocatable:: dux,duy,dvx,dvy 
   


   integer max_step,iter,interv

   double precision,dimension(0:nq-1)::t
   double precision,dimension(0:4)::tv
   integer,dimension(0:nq-1,dim)::e
   integer,dimension(0:4,dim)::ev
   !MRT
   double precision,dimension(0:nq-1,0:nq-1):: MM,MMT,Stau,SStau,stemp,tmp
   double precision,dimension(0:nq-1)::temp
contains

   !---------------------------------------------
   subroutine AllocateArrays
   !mpiio
   allocate(xx (local_start(1):local_end(1),local_start(2):local_end(2)))
   allocate(yy (local_start(1):local_end(1),local_start(2):local_end(2)))
   allocate(ppp (local_start(1):local_end(1),local_start(2):local_end(2)))
   allocate(ppp2 (local_start(1):local_end(1),local_start(2):local_end(2)))
   !distribution function
   allocate(f    (0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(fb   (0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(g1   (0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(gb1  (0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(g2   (0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(gb2  (0:nq-1,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))

   allocate(b11    (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(b12  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(b21  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(b22  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(bb11  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(bb12  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(bb21  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(bb22  (0:4,local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))

   !pressure
   allocate(p     (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(pp    (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dppx  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dppy  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dpx   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dpy   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))

   allocate(aa11  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(aa12  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(aa21   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(aa22   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))

   allocate(aaa11  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(aaa12  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(aaa21   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(aaa22   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))


   !relaxation time
   allocate(tau   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(tphi  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(nu    (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(mu    (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(rho   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dissip   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   


   !derivative
   allocate(divx1 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divx2 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divx3 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy1 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy2 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy3 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))

   allocate(divx11 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divx12 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divx21 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divx22 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy11 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy12 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy21 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(divy22 (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))


   allocate(sqd1  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(sqd2  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(sqd3  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(div1  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(div2  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(div3  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dx1   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dx2   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dx3   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dy1   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dy2   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dy3   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(drx   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dry   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(duy   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dvy   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dux   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(dvx   (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   

   
   !order parameter
   allocate(phi1  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(phi2  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(phi3  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(con  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   !velocity
   allocate(ua    (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(va    (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(rhou  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))
   allocate(rhov  (local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost))


   endsubroutine AllocateArrays

!-----------------------------------------------------
   subroutine DeAllocateArrays
   deallocate(xx)
   deallocate(yy)
   deallocate(ppp)
   deallocate(f)
   deallocate(fb)
   deallocate(g1)
   deallocate(gb1)
   deallocate(g2)
   deallocate(gb2)
   deallocate(b11)
   deallocate(b12)
   deallocate(b21)
   deallocate(b22)
   deallocate(bb11)
   deallocate(bb12)
   deallocate(bb21)
   deallocate(bb22)

   deallocate(p)
   deallocate(pp)
   deallocate(dpx)
   deallocate(dpy)
   deallocate(dppx)
   deallocate(dppy)

   deallocate(tau)
   deallocate(tphi)
   deallocate(nu)
   deallocate(mu)
   deallocate(rho)
   deallocate(dissip)
   deallocate(divx1)
   deallocate(divx2)
   deallocate(divx3)
   deallocate(divy1)
   deallocate(divy2)
   deallocate(divy3)
   deallocate(sqd1)
   deallocate(sqd2)
   deallocate(sqd3)
   deallocate(div1)
   deallocate(div2)
   deallocate(div3)
   deallocate(dx1)
   deallocate(dx2)
   deallocate(dx3)
   deallocate(dy1)
   deallocate(dy2)
   deallocate(dy3)
   deallocate(drx)
   deallocate(dry)
   deallocate(duy)
   deallocate(dvy)
   deallocate(dux)
   deallocate(dvx)


   deallocate(phi1)
   deallocate(phi2)
   deallocate(phi3)

   deallocate(ua)
   deallocate(va)
   deallocate(rhou)
   deallocate(rhov)


   endsubroutine DeAllocateArrays

!------------------------------------------------------  
   subroutine DataInput
   integer i,j,k
   open(100,file='input.in',status='old')
   !length and scaling
   read(100,*)nx,ny,delta
   !read(100,*)max_step,interv
   !fluid properties
   

   M=0.1d0
   cs2=1.d0/3.d0
   !density
   rho1=1.d0
   rho2=1.d0
   rho3=0.1d0

   !initialization
   r1=0.2d0*nx
   r0=0.2d0*nx

   x1=0.d0*nx
   y1=r1*0.2d0
   x0=0.d0*nx
   y0=y1+r1+delta



   t1=0.05d0
   t2=0.01
   t3=0.05d0
   tk=1e-6
   ratio=1e-6


   sigma12=1e-7
   sigma13=sigma12
   sigma23=sigma12/10

   sigma1=(sigma12+sigma13-sigma23)/2.d0
   sigma2=(sigma12+sigma23-sigma13)/2.d0
   sigma3=(sigma23+sigma13-sigma12)/2.d0

   Oh=rho2*t2*cs2/dsqrt(sigma12*rho2*r1)
   max_step=100


   interv=max_step/10
   endsubroutine DataInput

!-------------------------------------------------------------------------
!Test
   subroutine test
   if (rank.eq.0)then
      print*,'00000000'
      write(*,'(A24,F7.4)')   "1 Relaxation Time      :",t1
      write(*,'(A24,F7.4)')   "2 Relaxation Time      :",t2
      write(*,'(A24,F7.4)')   "3 Relaxation Time      :",t3
      write(*,'(A24,F10.8)')  "Surface tension 12     :",sigma12
      write(*,'(A24,F10.8)')  "Surface tension 13     :",sigma13
      write(*,'(A24,F10.8)')  "Surface tension 23     :",sigma23
      write(*,'(A24,F10.6)')  "OH                     :",Oh
      write(*,'(A24,F10.4)')  "Radius                 :",r1
      write(*,'(A24,I9.2)')   "Maximum # of Steps     :",max_step
      write(*,*)"==========================================="
   endif

   endsubroutine test

!---------------------------------------------
   subroutine Initilization

   integer i,j,iq,k



   !Lattice velocities and weights(2D)
   t(0)=4.d0/9.d0
   t(1)=1.d0/9.d0
   t(2)=1.d0/9.d0
   t(3)=1.d0/9.d0
   t(4)=1.d0/9.d0
   t(5)=1.d0/36.d0
   t(6)=1.d0/36.d0
   t(7)=1.d0/36.d0
   t(8)=1.d0/36.d0

   tv(0)=1.d0/3.d0
   tv(1)=1.d0/6.d0
   tv(2)=1.d0/6.d0
   tv(3)=1.d0/6.d0
   tv(4)=1.d0/6.d0

   ev(0,1)=0
   ev(1,1)=1
   ev(2,1)=0
   ev(3,1)=-1
   ev(4,1)=0

   ev(0,2)=0
   ev(1,2)=0
   ev(2,2)=1
   ev(3,2)=0
   ev(4,2)=-1

   e(0,1)=0
   e(1,1)=1
   e(2,1)=0
   e(3,1)=-1
   e(4,1)=0
   e(5,1)=1
   e(6,1)=-1
   e(7,1)=-1
   e(8,1)=1

   e(0,2)=0
   e(1,2)=0
   e(2,2)=1
   e(3,2)=0
   e(4,2)=-1
   e(5,2)=1
   e(6,2)=1
   e(7,2)=-1
   e(8,2)=-1


   call SetPhase

   if(rank == 0) Then
      open(200, file  =  'height.dat', status  =  'unknown')
      write(200,*),'"center"'
      close(200)
   endif



   !Initilize distributions
   call PassD(phi1)
   call Gradient(phi1,dx1,dy1)

   call PassD(phi2)
   call Gradient(phi2,dx2,dy2)
   
   call PassD(phi3)
   call Gradient(phi3,dx3,dy3)

   call PassD(p)
   call Gradient(p,dpx,dpy)

   call PassD(pp)
   call Gradient(pp,dppx,dppy)
   
   call PassD(ua)
   call Gradient(ua,dux,duy)

   call PassD(va)
   call Gradient(va,dvx,dvy)


      do j=local_start(2),local_end(2)
         do i=local_start(1),local_end(1)
            !calculation of curvature

            sqd1(i,j) = dsqrt(dx1(i,j)**2+dy1(i,j)**2)+1e-6

            sqd2(i,j) = dsqrt(dx2(i,j)**2+dy2(i,j)**2)+1e-6
   
            sqd3(i,j) = dsqrt(dx3(i,j)**2+dy3(i,j)**2)+1e-6

            divx1(i,j)=dx1(i,j)/sqd1(i,j)
            divy1(i,j)=dy1(i,j)/sqd1(i,j)
   

            divx2(i,j)=dx2(i,j)/sqd2(i,j)
            divy2(i,j)=dy2(i,j)/sqd2(i,j)
   

            divx3(i,j)=dx3(i,j)/sqd3(i,j)
            divy3(i,j)=dy3(i,j)/sqd3(i,j)

            drx(i,j)=dx1(i,j)*rho1+dx2(i,j)*rho2+dx3(i,j)*rho3
            dry(i,j)=dy1(i,j)*rho1+dy2(i,j)*rho2+dy3(i,j)*rho3
   
         enddo
      enddo




   call PassD(divx1)
   call PassD(divy1)
   call div(div1,divx1,divy1)

   call PassD(divx2)
   call PassD(divy2)
   call div(div2,divx2,divy2)

   call PassD(divx3)
   call PassD(divy3)
   call div(div3,divx3,divy3)



   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)


         do iq=0,4
            gamma5 =  3.d0 * ev(iq,1) * ua(i,j)+3.d0*ev(iq,2)*va(i,j)

            b11(iq,i,j)=tv(iq) * aa11(i,j) * (1.d0+ gamma5)
            b12(iq,i,j)=tv(iq) * aa12(i,j) * (1.d0+ gamma5)
            b21(iq,i,j)=tv(iq) * aa21(i,j) * (1.d0+ gamma5)
            b22(iq,i,j)=tv(iq) * aa22(i,j) * (1.d0+ gamma5)
         enddo

      enddo
   enddo

   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         aaa11(i,j)=phi2(i,j)*(aa11(i,j)-1.d0)
         aaa12(i,j)=phi2(i,j)*aa12(i,j)
         aaa21(i,j)=phi2(i,j)*aa21(i,j)
         aaa22(i,j)=phi2(i,j)*(aa22(i,j)-1.d0)

         divx11(i,j)=0.d0
         divx12(i,j)=0.d0
         divx21(i,j)=0.d0
         divx22(i,j)=0.d0

         divy11(i,j)=0.d0
         divy12(i,j)=0.d0
         divy21(i,j)=0.d0
         divy22(i,j)=0.d0

      enddo
   enddo



   
      do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         psq=phi1(i,j)**2+phi2(i,j)**2+phi3(i,j)**2

         forsx1=-3.d0/2.d0*delta*sigma1*div1(i,j)*sqd1(i,j)*dx1(i,j)  
         forsy1=-3.d0/2.d0*delta*sigma1*div1(i,j)*sqd1(i,j)*dy1(i,j)
   
         forsx2=-3.d0/2.d0*delta*sigma2*div2(i,j)*sqd2(i,j)*dx2(i,j)  
         forsy2=-3.d0/2.d0*delta*sigma2*div2(i,j)*sqd2(i,j)*dy2(i,j)
   
         forsx3=-3.d0/2.d0*delta*sigma3*div3(i,j)*sqd3(i,j)*dx3(i,j)  
         forsy3=-3.d0/2.d0*delta*sigma3*div3(i,j)*sqd3(i,j)*dy3(i,j)                 
         

         forspx=-dpx(i,j)
         forspy=-dpy(i,j)

         forsppx=dppx(i,j)
         forsppy=dppy(i,j)
         
         forsnx=nu(i,j)*(2*dux(i,j)*drx(i,j)+(duy(i,j)+dvx(i,j))*dry(i,j))
         forsny=nu(i,j)*(2*dvy(i,j)*dry(i,j)+(duy(i,j)+dvx(i,j))*drx(i,j))

         forsax=0.d0!ratio*(divx11(i,j)+divy21(i,j))
         forsay=0.d0!ratio*(divx12(i,j)+divy22(i,j))


         sqad1 = 4.d0*phi1(i,j)*(1.d0-phi1(i,j))/delta
         sqad2 = 4.d0*phi2(i,j)*(1.d0-phi2(i,j))/delta
         sqad3 = 4.d0*phi3(i,j)*(1.d0-phi3(i,j))/delta

         forphix1 =  (1.d0-phi1(i,j)**2/psq)*sqad1*divx1(i,j)-&
         phi1(i,j)**2/psq*sqad2*divx2(i,j)-&
         phi1(i,j)**2/psq*sqad3*divx3(i,j)
 
         forphiy1 =  (1.d0-phi1(i,j)**2/psq)*sqad1*divy1(i,j)-&
         phi1(i,j)**2/psq*sqad2*divy2(i,j)-&
         phi1(i,j)**2/psq*sqad3*divy3(i,j)
                                         
 
         forphix2 =  (1.d0-phi2(i,j)**2/psq)*sqad2*divx2(i,j)-&
         phi2(i,j)**2/psq*sqad1*divx1(i,j)-&
         phi2(i,j)**2/psq*sqad3*divx3(i,j)
 
         forphiy2 =  (1.d0-phi2(i,j)**2/psq)*sqad2*divy2(i,j)-&
         phi2(i,j)**2/psq*sqad1*divy1(i,j)-&
         phi2(i,j)**2/psq*sqad3*divy3(i,j)


         do iq=0,nq-1
            
            gamma =  3.d0*(e(iq,1) * ua(i,j)+     &
                           e(iq,2) * va(i,j))+    &
                     4.5d0*(e(iq,1) * ua(i,j)+     &
                            e(iq,2) * va(i,j))**2- &
                     1.5d0*(ua(i,j)**2+va(i,j)**2)


            fors =   (e(iq,1)-ua(i,j))*(1.d0+gamma)*t(iq)/rho(i,j)*&
                     (forspx+forsnx+forsx1+forsx2+forsx3)+&
                     (e(iq,2)-va(i,j))*(1.d0+Gamma)*t(iq)/rho(i,j)*&
                     (forspy+forsny+forsy1+forsy2+forsy3)+&
                     (e(iq,1)-ua(i,j))*t(iq)*forsppx+&
                     (e(iq,2)-va(i,j))*t(iq)*forsppy
         
            forphi1 = t(iq)* (1.d0+gamma)*  &
                     ((e(iq,1)-ua(i,j))*forphix1+ &
                     (e(iq,2)-va(i,j))*forphiy1)

            forphi2 = t(iq)* (1.d0+gamma)*  &
                     ((e(iq,1)-ua(i,j))*forphix2+ &
                     (e(iq,2)-va(i,j))*forphiy2)


             g1(iq,i,j) = t(iq) * phi1(i,j) * (1.d0 + gamma) - 0.5d0 * forphi1
             g2(iq,i,j) = t(iq) * phi2(i,j) * (1.d0 + gamma) - 0.5d0 * forphi2
             f(iq,i,j)  = t(iq) * pp(i,j) + t(iq) * gamma * cs2 - 0.5d0 * fors

         enddo
      enddo
   enddo

   endsubroutine Initilization

!-------------------------------------------------------------------


!-----------------------------------------------------------
   subroutine SetPhase
   double precision dist,dist1
   integer i,j,iq

   do j=local_start(2)-ghost,local_end(2)+ghost
      do i=local_start(1)-ghost,local_end(1)+ghost
         do iq=0,4
            b11(iq,i,j)=0.d0
            bb11(iq,i,j)=0.d0

            b12(iq,i,j)=0.d0
            bb12(iq,i,j)=0.d0

            b21(iq,i,j)=0.d0
            bb21(iq,i,j)=0.d0

            b22(iq,i,j)=0.d0
            bb22(iq,i,j)=0.d0
         enddo
         aa11(i,j)=1.d0
         aa12(i,j)=0.d0
         aa21(i,j)=0.d0
         aa22(i,j)=1.d0

         aaa11(i,j)=0.d0
         aaa12(i,j)=0.d0
         aaa21(i,j)=0.d0
         aaa22(i,j)=0.d0
      enddo
   enddo
   do j=local_start(2)-ghost,local_end(2)+ghost
      do i=local_start(1)-ghost,local_end(1)+ghost
         do iq=0,nq-1
            f(iq,i,j)=0.d0
            fb(iq,i,j)=0.d0

            g1(iq,i,j)=0.d0
            gb1(iq,i,j)=0.d0

            g2(iq,i,j)=0.d0
            gb2(iq,i,j)=0.d0

            
         enddo
         dist=2.d0*(j-y1)/delta

         dist1=2.d0*(sqrt((i-x0)**2+(j-y0)**2+0.d0)-r0)/delta

         phi1(i,j)=0.5d0-0.5d0*tanh(dist1)

         phi2(i,j)=0.5d0-0.5d0*tanh(dist)

         phi3(i,j)=1.d0-phi1(i,j)-phi2(i,j)

         con(i,j)=-phi1(i,j)+phi3(i,j)

         tau(i,j)=t1*phi1(i,j)+t3*phi3(i,j)+t2*phi2(i,j)
         tphi(i,j)=M/cs2
         rho(i,j)=rho1*phi1(i,j)+rho2*phi2(i,j)+rho3*phi3(i,j)
         nu(i,j)=tau(i,j)*cs2
         mu(i,j)=nu(i,j)*rho(i,j)

         p(i,j)=0.d0
         pp(i,j)=0.d0
         ua(i,j)=0.d0
         va(i,j)=0.d0
         rhou(i,j)=0.d0
         rhov(i,j)=0.d0

      enddo
   enddo

   endsubroutine SetPhase


!---------------------------------------------
   subroutine Collision
   integer i,j,iq,k

   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         forsa11=2.d0*aa11(i,j)*dux(i,j)+aa21(i,j)*dvx(i,j)+aa12(i,j)*dvx(i,j)
         forsa12=aa12(i,j)*dux(i,j)+aa11(i,j)*duy(i,j)+aa22(i,j)*dvx(i,j)+aa12(i,j)*dvy(i,j)
         forsa21=aa21(i,j)*dux(i,j)+aa11(i,j)*duy(i,j)+aa22(i,j)*dvx(i,j)+aa21(i,j)*dvy(i,j)
         forsa22=2.d0*aa22(i,j)*dvy(i,j)+aa21(i,j)*duy(i,j)+aa12(i,j)*duy(i,j)
         do iq=0,4
            gamma5 =  3.d0*(ev(iq,1) * ua(i,j)+     &
                           ev(iq,2) * va(i,j))

            aaeq11=tv(iq) * aa11(i,j) * (1.d0+ gamma5)- 0.5d0 * forsa11*tv(iq)* (1.d0+gamma5)
            aaeq12=tv(iq) * aa12(i,j) * (1.d0+ gamma5)- 0.5d0 * forsa12*tv(iq)* (1.d0+gamma5)
            aaeq21=tv(iq) * aa21(i,j) * (1.d0+ gamma5)- 0.5d0 * forsa21*tv(iq)* (1.d0+gamma5)
            aaeq22=tv(iq) * aa22(i,j) * (1.d0+ gamma5)- 0.5d0 * forsa22*tv(iq)* (1.d0+gamma5)


            bb11(iq,i,j)=b11(iq,i,j)-1.d0/(0.5d0+tk)*(b11(iq,i,j)-aaeq11)+forsa11*tv(iq)* (1.d0+gamma5)
            bb12(iq,i,j)=b12(iq,i,j)-1.d0/(0.5d0+tk)*(b12(iq,i,j)-aaeq12)+forsa12*tv(iq)* (1.d0+gamma5)
            bb21(iq,i,j)=b21(iq,i,j)-1.d0/(0.5d0+tk)*(b21(iq,i,j)-aaeq21)+forsa21*tv(iq)* (1.d0+gamma5)
            bb22(iq,i,j)=b22(iq,i,j)-1.d0/(0.5d0+tk)*(b22(iq,i,j)-aaeq22)+forsa22*tv(iq)* (1.d0+gamma5)
         enddo

      enddo
   enddo
   
   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         psq=phi1(i,j)**2+phi2(i,j)**2+phi3(i,j)**2

         forsx1=-3.d0/2.d0*delta*sigma1*div1(i,j)*sqd1(i,j)*dx1(i,j)  
         forsy1=-3.d0/2.d0*delta*sigma1*div1(i,j)*sqd1(i,j)*dy1(i,j)
   
         forsx2=-3.d0/2.d0*delta*sigma2*div2(i,j)*sqd2(i,j)*dx2(i,j)  
         forsy2=-3.d0/2.d0*delta*sigma2*div2(i,j)*sqd2(i,j)*dy2(i,j)
   
         forsx3=-3.d0/2.d0*delta*sigma3*div3(i,j)*sqd3(i,j)*dx3(i,j)  
         forsy3=-3.d0/2.d0*delta*sigma3*div3(i,j)*sqd3(i,j)*dy3(i,j)            
            


         forspx=-dpx(i,j)
         forspy=-dpy(i,j)

         forsppx=dppx(i,j)
         forsppy=dppy(i,j)
         
         forsnx=nu(i,j)*(2*dux(i,j)*drx(i,j)+(duy(i,j)+dvx(i,j))*dry(i,j))
         forsny=nu(i,j)*(2*dvy(i,j)*dry(i,j)+(duy(i,j)+dvx(i,j))*drx(i,j))

         forsax=ratio*(divx11(i,j)+divy21(i,j))
         forsay=ratio*(divx12(i,j)+divy22(i,j))

         sqad1 = 4.d0*phi1(i,j)*(1.d0-phi1(i,j))/delta
         sqad2 = 4.d0*phi2(i,j)*(1.d0-phi2(i,j))/delta
         sqad3 = 4.d0*phi3(i,j)*(1.d0-phi3(i,j))/delta

         forphix1 =  (1.d0-phi1(i,j)**2/psq)*sqad1*divx1(i,j)-&
         phi1(i,j)**2/psq*sqad2*divx2(i,j)-&
         phi1(i,j)**2/psq*sqad3*divx3(i,j)

         forphiy1 =  (1.d0-phi1(i,j)**2/psq)*sqad1*divy1(i,j)-&
         phi1(i,j)**2/psq*sqad2*divy2(i,j)-&
         phi1(i,j)**2/psq*sqad3*divy3(i,j)


         forphix2 =  (1.d0-phi2(i,j)**2/psq)*sqad2*divx2(i,j)-&
         phi2(i,j)**2/psq*sqad1*divx1(i,j)-&
         phi2(i,j)**2/psq*sqad3*divx3(i,j)

         forphiy2 =  (1.d0-phi2(i,j)**2/psq)*sqad2*divy2(i,j)-&
         phi2(i,j)**2/psq*sqad1*divy1(i,j)-&
         phi2(i,j)**2/psq*sqad3*divy3(i,j)


         do iq=0,nq-1

            gamma =  3.d0*(e(iq,1) * ua(i,j)+     &
                           e(iq,2) * va(i,j))+    &
                     4.5d0*(e(iq,1) * ua(i,j)+     &
                            e(iq,2) * va(i,j))**2- &
                     1.5d0*(ua(i,j)**2+va(i,j)**2)


            fors =   (e(iq,1)-ua(i,j))*(1.d0+gamma)*t(iq)/rho(i,j)*&
                     (forspx+forsnx+forsx1+forsx2+forsx3+forsax)+&
                     (e(iq,2)-va(i,j))*(1.d0+Gamma)*t(iq)/rho(i,j)*&
                     (forspy+forsny+forsy1+forsy2+forsy3+forsay)+&
                     (e(iq,1)-ua(i,j))*t(iq)*forsppx+&
                     (e(iq,2)-va(i,j))*t(iq)*forsppy
         
            forphi1 = t(iq)* (1.d0+gamma)*  &
                     ((e(iq,1)-ua(i,j))*forphix1+ &
                     (e(iq,2)-va(i,j))*forphiy1)

            forphi2 = t(iq)* (1.d0+gamma)*  &
                     ((e(iq,1)-ua(i,j))*forphix2+ &
                     (e(iq,2)-va(i,j))*forphiy2)      
      
               
            geq1 = t(iq) * phi1(i,j) * (1.d0 + gamma) - 0.5d0 * forphi1
            geq2 = t(iq) * phi2(i,j) * (1.d0 + gamma) - 0.5d0 * forphi2
            feq  = t(iq) * pp(i,j) + t(iq) * gamma * cs2 - 0.5d0 * fors


            gb1(iq,i,j)=g1(iq,i,j)-1.d0/(0.5d0+tphi(i,j))*&
            (g1(iq,i,j)-geq1)+forphi1   

            gb2(iq,i,j)=g2(iq,i,j)-1.d0/(0.5d0+tphi(i,j))*&
            (g2(iq,i,j)-geq2)+forphi2

            fb(iq,i,j)=f(iq,i,j)-1.d0/(0.5d0+tau(i,j))*(f(iq,i,j)-feq)+fors
   
         enddo
      enddo
   enddo
   endsubroutine Collision
!---------------------------------------------

   subroutine BoundaryCondition
   integer i,j,iq,io
   call PassF(fb)
   call PassF(gb1)
   call PassF(gb2)

   call PassV(bb11)
   call PassV(bb12)
   call PassV(bb21)
   call PassV(bb22)


   if(nx.eq.local_end(1))then
      do j=local_start(2)-1,local_end(2)+1
         fb (3,nx+1,j) = fb (1,nx-1,j)
         gb1(3,nx+1,j) = gb1(1,nx-1,j)
         gb2(3,nx+1,j) = gb2(1,nx-1,j)

         fb (6,nx+1,j) = fb (5,nx-1,j)
         gb1(6,nx+1,j) = gb1(5,nx-1,j)
         gb2(6,nx+1,j) = gb2(5,nx-1,j)
      
         fb (7,nx+1,j) = fb (8,nx-1,j)
         gb1(7,nx+1,j) = gb1(8,nx-1,j)
         gb2(7,nx+1,j) = gb2(8,nx-1,j)

         !bb11(3,nx+1,j)= -bb11(0,nx,j)*(1.d0-1.d0/tv(0))-bb11(1,nx-1,j)-bb11(2,nx,j-1)-bb11(4,nx,j+1)
         !bb12(3,nx+1,j)= -bb12(0,nx,j)*(1.d0-1.d0/tv(0))-bb12(1,nx-1,j)-bb12(2,nx,j-1)-bb12(4,nx,j+1)
         !bb21(3,nx+1,j)= -bb21(0,nx,j)*(1.d0-1.d0/tv(0))-bb21(1,nx-1,j)-bb21(2,nx,j-1)-bb21(4,nx,j+1)
         !bb22(3,nx+1,j)= -bb22(0,nx,j)*(1.d0-1.d0/tv(0))-bb22(1,nx-1,j)-bb22(2,nx,j-1)-bb22(4,nx,j+1)
      

      enddo

   elseif(0.eq.local_start(1))then
      do j=local_start(2)-1,local_end(2)+1
         fb (1,-1,j) = fb (3,1,j)
         gb1(1,-1,j) = gb1(3,1,j)
         gb2(1,-1,j) = gb2(3,1,j)

         fb (5,-1,j) = fb (6,1,j)
         gb1(5,-1,j) = gb1(6,1,j)
         gb2(5,-1,j) = gb2(6,1,j)
      
         fb (8,-1,j) = fb (7,1,j)
         gb1(8,-1,j) = gb1(7,1,j)
         gb2(8,-1,j) = gb2(7,1,j)

         !bb11(1,-1,j)= -bb11(0,0,j)*(1.d0-1.d0/tv(0))-bb11(3,1,j)-bb11(2,0,j-1)-bb11(4,0,j+1)
         !bb12(1,-1,j)= -bb12(0,0,j)*(1.d0-1.d0/tv(0))-bb12(3,1,j)-bb12(2,0,j-1)-bb12(4,0,j+1)
         !bb21(1,-1,j)= -bb21(0,0,j)*(1.d0-1.d0/tv(0))-bb21(3,1,j)-bb21(2,0,j-1)-bb21(4,0,j+1)
         !bb22(1,-1,j)= -bb22(0,0,j)*(1.d0-1.d0/tv(0))-bb22(3,1,j)-bb22(2,0,j-1)-bb22(4,0,j+1)
      

      
      enddo
   endif

   if(ny.eq.local_end(2))then
      do i=local_start(1)-1,local_end(1)+1
         fb (4,i,ny+1) = fb (2,i,ny-1)
         gb1(4,i,ny+1) = gb1(2,i,ny-1)
         gb2(4,i,ny+1) = gb2(2,i,ny-1)
      
         fb (7,i,ny+1) = fb (5,i-2,ny-1)
         gb1(7,i,ny+1) = gb1(5,i-2,ny-1)
         gb2(7,i,ny+1) = gb2(5,i-2,ny-1)
      
         fb (8,i,ny+1) = fb (6,i+2,ny-1)
         gb1(8,i,ny+1) = gb1(6,i+2,ny-1)
         gb2(8,i,ny+1) = gb2(6,i+2,ny-1)

         !bb11(4,i,ny+1)=  -bb11(0,i,ny)*(1.d0-1.d0/tv(0))-bb11(1,i-1,ny)-bb11(2,i,ny-1)-bb11(3,i+1,ny)
         !bb12(4,i,ny+1)=  -bb12(0,i,ny)*(1.d0-1.d0/tv(0))-bb12(1,i-1,ny)-bb12(2,i,ny-1)-bb12(3,i+1,ny)
         !bb21(4,i,ny+1)=  -bb21(0,i,ny)*(1.d0-1.d0/tv(0))-bb21(1,i-1,ny)-bb21(2,i,ny-1)-bb21(3,i+1,ny)
         !bb22(4,i,ny+1)=  -bb22(0,i,ny)*(1.d0-1.d0/tv(0))-bb22(1,i-1,ny)-bb22(2,i,ny-1)-bb22(3,i+1,ny)
      enddo

   elseif(0.eq.local_start(2))then
      do i=local_start(1)-1,local_end(1)+1
         fb (2,i,-1) = fb (4,i,1)
         gb1(2,i,-1) = gb1(4,i,1)
         gb2(2,i,-1) = gb2(4,i,1)
      
         fb (5,i,-1) = fb (7,i+2,1)
         gb1(5,i,-1) = gb1(7,i+2,1)
         gb2(5,i,-1) = gb2(7,i+2,1)
      
         fb (6,i,-1) = fb (8,i-2,1)
         gb1(6,i,-1) = gb1(8,i-2,1)
         gb2(6,i,-1) = gb2(8,i-2,1)

         !bb11(2,i,-1)= -bb11(0,i,0)*(1.d0-1.d0/tv(0))-bb11(1,i-1,0)-bb11(3,i+1,0)-bb11(4,i,1)
         !bb12(2,i,-1)= -bb12(0,i,0)*(1.d0-1.d0/tv(0))-bb12(1,i-1,0)-bb12(3,i+1,0)-bb12(4,i,1)
         bb21(2,i,-1)= -bb21(0,i,0)*(1.d0-1.d0/tv(0))-bb21(1,i-1,0)-bb21(3,i+1,0)-bb21(4,i,1)
         bb22(2,i,-1)= -bb22(0,i,0)*(1.d0-1.d0/tv(0))-bb22(1,i-1,0)-bb22(3,i+1,0)-bb22(4,i,1)
      
      enddo
   endif


   


   call MPI_BARRIER(CART_COMM,ierr)
   endsubroutine BoundaryCondition

!------------------------------------------------
   subroutine Propagation

   integer i,j,iq

   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)      
         do iq=0,nq-1            
            f(iq,i,j)=fb(iq,i-e(iq,1),j-e(iq,2))  
            g1(iq,i,j)=gb1(iq,i-e(iq,1),j-e(iq,2))         
            g2(iq,i,j)=gb2(iq,i-e(iq,1),j-e(iq,2))   
         enddo
      enddo
   enddo

   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)      
         do iq=0,4          
            b11(iq,i,j)=bb11(iq,i-ev(iq,1),j-ev(iq,2))  
            b12(iq,i,j)=bb12(iq,i-ev(iq,1),j-ev(iq,2))         
            b21(iq,i,j)=bb21(iq,i-ev(iq,1),j-ev(iq,2))   
            b22(iq,i,j)=bb22(iq,i-ev(iq,1),j-ev(iq,2))   
         enddo
      enddo
   enddo
   endsubroutine Propagation

!---------------------------------------------
   subroutine PostProcessing
      integer i,j,iq

   !PHI and RHO


   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)
         phi1(i,j)=0.d0
         phi2(i,j)=0.d0
         do iq=0,nq-1
            phi1(i,j)=phi1(i,j)+g1(iq,i,j)
            phi2(i,j)=phi2(i,j)+g2(iq,i,j)
         enddo
         if(local_start(2).eq.0)then
            phi1(i,0)=0.d0
            phi2(i,0)=1.d0
         endif
         phi3(i,j)=1.d0-phi1(i,j)-phi2(i,j)
         con(i,j)=-phi1(i,j)+phi3(i,j)

         rho(i,j)=rho1*phi1(i,j)+rho2*phi2(i,j)+rho3*phi3(i,j)

         tau(i,j)=t1*phi1(i,j)+t3*phi3(i,j)+t2*phi2(i,j)
         nu(i,j)=tau(i,j)*cs2
         mu(i,j)=nu(i,j)*rho(i,j)

      enddo
   enddo

   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         if(rho(i,j)<rho3)then
            rho(i,j)=rho3
         else if(rho(i,j)>rho1)then
            rho(i,j)=rho1
         endif

      enddo
   enddo




   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         pp(i,j)=0.d0
         do iq=0,nq-1
            pp(i,j)=pp(i,j)+f(iq,i,j)
         enddo
         pp(i,j)=pp(i,j)-ua(i,j)/2.d0*dppx(i,j)-va(i,j)/2.d0*dppy(i,j)
         p(i,j)=pp(i,j)*rho(i,j)
      enddo
   enddo



   call PassD(phi1)
   call Gradient(phi1,dx1,dy1)
   
   call PassD(phi2)
   call Gradient(phi2,dx2,dy2)
   
   call PassD(phi3)
   call Gradient(phi3,dx3,dy3)

   call PassD(p)
   call Gradient(p,dpx,dpy)

   call PassD(pp)
   call Gradient(pp,dppx,dppy)
   
   call PassD(ua)
   call Gradient(ua,dux,duy)

   call PassD(va)
   call Gradient(va,dvx,dvy)



   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         sqd1(i,j) =dsqrt(dx1(i,j)**2+dy1(i,j)**2)+1e-6
         divx1(i,j)=dx1(i,j)/sqd1(i,j)
         divy1(i,j)=dy1(i,j)/sqd1(i,j)

         sqd2(i,j) = dsqrt(dx2(i,j)**2+dy2(i,j)**2)+1e-6
         divx2(i,j)=dx2(i,j)/sqd2(i,j)
         divy2(i,j)=dy2(i,j)/sqd2(i,j)

         sqd3(i,j) = dsqrt(dx3(i,j)**2+dy3(i,j)**2)+1e-6
         divx3(i,j)=dx3(i,j)/sqd3(i,j)
         divy3(i,j)=dy3(i,j)/sqd3(i,j)

         drx(i,j)=dx1(i,j)*rho1+dx2(i,j)*rho2+dx3(i,j)*rho3
         dry(i,j)=dy1(i,j)*rho1+dy2(i,j)*rho2+dy3(i,j)*rho3

      enddo
   enddo
   call PassD(divx1)
   call PassD(divy1)
   call div(div1,divx1,divy1)

   call PassD(divx2)
   call PassD(divy2)
   call div(div2,divx2,divy2)

   call PassD(divx3)
   call PassD(divy3)
   call div(div3,divx3,divy3)

      do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)
         rhou(i,j)=0.d0
         rhov(i,j)=0.d0

         do iq=0,nq-1
            rhou(i,j)=rhou(i,j)+f(iq,i,j)*e(iq,1)
            rhov(i,j)=rhov(i,j)+f(iq,i,j)*e(iq,2)
         enddo

      enddo
   enddo

   do j=local_start(2),local_end(2)
   do i=local_start(1),local_end(1)

         forsx1=-3.d0/2.d0*delta*sigma1*div1(i,j)*sqd1(i,j)*dx1(i,j)  
         forsy1=-3.d0/2.d0*delta*sigma1*div1(i,j)*sqd1(i,j)*dy1(i,j)
   
         forsx2=-3.d0/2.d0*delta*sigma2*div2(i,j)*sqd2(i,j)*dx2(i,j)  
         forsy2=-3.d0/2.d0*delta*sigma2*div2(i,j)*sqd2(i,j)*dy2(i,j)
   
         forsx3=-3.d0/2.d0*delta*sigma3*div3(i,j)*sqd3(i,j)*dx3(i,j)  
         forsy3=-3.d0/2.d0*delta*sigma3*div3(i,j)*sqd3(i,j)*dy3(i,j)            
            

         forspx=-dpx(i,j)
         forspy=-dpy(i,j)

         forsppx=dppx(i,j)*rho(i,j)
         forsppy=dppy(i,j)*rho(i,j)
         
         forsnx=nu(i,j)*(2*dux(i,j)*drx(i,j)+(duy(i,j)+dvx(i,j))*dry(i,j))
         forsny=nu(i,j)*(2*dvy(i,j)*dry(i,j)+(duy(i,j)+dvx(i,j))*drx(i,j)) 
         forsax=ratio*(divx11(i,j)+divy21(i,j))
         forsay=ratio*(divx12(i,j)+divy22(i,j))

         ua(i,j)=rhou(i,j)/cs2+(forspx+forsppx+forsnx+forsx1+forsx2+forsx3+forsax)/2.d0/rho(i,j)
         va(i,j)=rhov(i,j)/cs2+(forspy+forsppy+forsny+forsy1+forsy2+forsy3+forsay)/2.d0/rho(i,j)
      enddo          
   enddo
   call PassD(ua)
   call Gradient(ua,dux,duy)

   call PassD(va)
   call Gradient(va,dvx,dvy)
   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)
         aa11(i,j)=0.d0
         aa12(i,j)=0.d0
         aa21(i,j)=0.d0
         aa22(i,j)=0.d0

         do iq=0,4
            aa11(i,j)=aa11(i,j)+b11(iq,i,j)
            aa12(i,j)=aa12(i,j)+b12(iq,i,j)
            aa21(i,j)=aa21(i,j)+b21(iq,i,j)
            aa22(i,j)=aa22(i,j)+b22(iq,i,j)
         enddo
         forsa11=2.d0*aa11(i,j)*dux(i,j)+aa21(i,j)*dvx(i,j)+aa12(i,j)*dvx(i,j)
         forsa12=aa12(i,j)*dux(i,j)+aa11(i,j)*duy(i,j)+aa22(i,j)*dvx(i,j)+aa12(i,j)*dvy(i,j)
         forsa21=aa21(i,j)*dux(i,j)+aa11(i,j)*duy(i,j)+aa22(i,j)*dvx(i,j)+aa21(i,j)*dvy(i,j)
         forsa22=2.d0*aa22(i,j)*dvy(i,j)+aa21(i,j)*duy(i,j)+aa12(i,j)*duy(i,j)

         aa11(i,j)=aa11(i,j)+forsa11/2.d0
         aa12(i,j)=aa12(i,j)+forsa12/2.d0
         aa21(i,j)=aa21(i,j)+forsa21/2.d0
         aa22(i,j)=aa22(i,j)+forsa22/2.d0

      enddo
   enddo

   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         aaa11(i,j)=phi2(i,j)*(aa11(i,j)-1.d0)
         aaa12(i,j)=phi2(i,j)*aa12(i,j)
         aaa21(i,j)=phi2(i,j)*aa21(i,j)
         aaa22(i,j)=phi2(i,j)*(aa22(i,j)-1.d0)
      enddo
   enddo

   call PassD(aaa11)
   call PassD(aaa12)
   call PassD(aaa21)
   call PassD(aaa22)

   call Gradient(aaa11,divx11,divy11)
   call Gradient(aaa12,divx12,divy12)
   call Gradient(aaa21,divx21,divy21)
   call Gradient(aaa22,divx22,divy22)


   endsubroutine PostProcessing 
!----------------------------------------




!------------------------------------

   subroutine MPIoutput
      
      INTEGER::i,j
      INTEGER::FILETYPE,DATATYPE
      INTEGER(KIND=MPI_OFFSET_KIND):: offset
      INTEGER:: Status(MPI_STATUS_SIZE)
      INTEGER::ierr,fh
      DOUBLE PRECISION::umin,umax,vmin,vmax,rhomin,rhomax,pmin,pmax,ppmin,ppmax,disum,dismin,dismax
      DOUBLE PRECISION::phi1min,phi1max,phi2min,phi2max,phi3min,phi3max,conmin,conmax,psum,psum2,dissum
      INTEGER::tp
      REAL::ttp
      Double precision::fp
      INTEGER,DIMENSION(dim):: start
      CHARACTER:: version*8,NULCHAR*1,filename*50,num*4
      ppp=0.d0
      do j=local_start(2),local_end(2)
         do i=local_start(1),local_end(1)
            xx(i,j)=0.5*(i+0.d0)/nx
            yy(i,j)=0.6*(j+0.d0)/ny
            if(i.ne.0)then
            ppp(i,j) = j*phi1(i,j)*2
            ppp2(i,j)= j*phi2(i,j)*2
            else
               ppp(i,j) = j*phi1(i,j)
               ppp2(i,j)= j*phi2(i,j)
            endif
         enddo
      enddo
      psum=sum(ppp)
      psum2=sum(ppp2)
      dissum=sum(dissip)

      
      call MPI_BARRIER(CART_COMM,ierr)
      call mpi_reduce(minval(phi1), phi1min, 1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(phi1), phi1max, 1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(aa11), phi2min, 1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(aa11), phi2max, 1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(aa12), phi3min, 1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(aa12), phi3max, 1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(con),  conmin,  1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(con),  conmax,  1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(rho),  rhomin,  1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(rho),  rhomax,  1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(ua),   umin,    1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(ua),   umax,    1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(va),   vmin,    1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(va),   vmax,    1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(p),    pmin,    1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(p),    pmax,    1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(aaa11),   ppmin,   1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(aaa11),   ppmax,   1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      call mpi_reduce(minval(dissip),   dismin,   1,MPI_double_precision,mpi_min,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(maxval(dissip),   dismax,   1,MPI_double_precision,mpi_max,0,MPI_COMM_WORLD,ierr)

      start=(/ghost,ghost/)
      offset=0




      call MPI_TYPE_CREATE_SUBARRAY(2,local_length+2*ghost,local_length,start,MPI_ORDER_FORTRAN,&
      MPI_double_precision,datatype,ierr)
      call MPI_TYPE_COMMIT(datatype,ierr)
      call MPI_TYPE_CREATE_SUBARRAY(2,global_length,local_length,local_start,MPI_ORDER_FORTRAN,&
      MPI_double_precision,filetype,ierr)
      call MPI_TYPE_COMMIT(filetype,ierr)

      write(num,'(i4.4)')iter/interv
      filename='rho'//num(1:4)//'.plt'
      call MPI_File_open(MPI_COMM_WORLD,filename, &
      MPI_MODE_WRONLY + MPI_MODE_CREATE,MPI_INFO_NULL,fh,ierr)

      if(rank == 0) Then

         call MPI_File_seek (fh,offset,MPI_seek_set,ierr)
         !header
         !version number
         version='#!TDV112'
         call MPI_File_write(fh,version,8,MPI_char,status,ierr)
         !INTEGER 1
         tp=1
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !TITLE AND VARIABLE NAME
         !FULL FILE TYPE
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !TITLE
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !3*4+8=20
         !VARIABLE NAME
         tp=11
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('x'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('y'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('u'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('v'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('p'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('h'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('i'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('1'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('p'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('h'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('i'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('2'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('p'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('h'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('i'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('3'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('c'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('o'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('n'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('p'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('r'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('h'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('o'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('d'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('s'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         !20+4*37=168
         !ZONE MARKER
         ttp=299.0
         call MPI_File_write(fh,ttp,            1,MPI_REAL,status,ierr)
         !ZONE NAME
         call MPI_File_write(fh,ichar('Z'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('O'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('N'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('E'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(' '),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('0'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('0'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar('1'),     1,MPI_integer,status,ierr)
         call MPI_File_write(fh,ichar(nulchar), 1,MPI_integer,status,ierr)
         !168+10*4=208
         !PARAENTS
         tp=-1
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !STRENDID
         tp=-2
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !SOLUTION TIME
         fp=iter/interv
         call MPI_File_write(fh,fp,             1,MPI_double_precision,status,ierr)
         !ZONE COLOR
         tp=-1
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !ZONE TYPE
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !SPECIFY VAR LOCATION
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !ARE RAW LOCAL
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !NUMBER OF MISCELLANEOUS
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !ORDERED ZONE
         tp=nx+1
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         tp=ny+1
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         tp=1
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !AUXILIARY
         tp=0
         call MPI_File_write(fh,tp,             1,MPI_integer,status,ierr)
         !208+13*4=260
         !EOHMARKER
         ttp=357.0
         call MPI_File_write(fh,ttp,             1,MPI_REAL,status,ierr)

         !DATA SECTION
         ttp=299.0
         call MPI_File_write(fh,ttp,             1,MPI_REAL,status,ierr)
         !VARIABLE DATA FORMAT
         tp=2
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         call MPI_File_write(fh,tp,              1,MPI_integer,status,ierr)
         
         !PASSIVE VARIABLE
         tp=0
         call MPI_File_write(fh,tp,                1,MPI_integer,status,ierr)
         !SHARING VARIABLE
         call MPI_File_write(fh,tp,                1,MPI_integer,status,ierr)
         !ZONE NUMBER
         tp=-1
         call MPI_File_write(fh,tp,                1,MPI_integer,status,ierr)
         !260+16*4=324
         !MIN AND MAX VALUE FLOAT 64
         fp=0.d0
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=0.5d0
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=0.d0
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=0.6d0
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=umin
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=umax
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=vmin
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=vmax
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=phi1min
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=phi1max
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=phi2min
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=phi2max
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=phi3min
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=phi3max
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=conmin
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=conmax
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=pmin
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=pmax
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=rhomin
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=rhomax
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=dismin
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         fp=dismax
         call MPI_File_write(fh,fp, 1,MPI_double_precision,status,ierr)
         !324+22*8=500
      endif

      offset=500

      call MPI_File_set_view(fh,offset,MPI_double_precision,filetype,&
      "native",MPI_INFO_NULL,ierr)
      !DATA
      CALL MPI_FILE_WRITE_all(fh, xx, local_length(1)*local_length(2), MPI_double_precision, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, yy, local_length(1)*local_length(2), MPI_double_precision, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, ua   ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, va   ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, phi1 ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, aa11 ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, aa12 ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, con  ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, p    ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, rho  ,1, datatype, MPI_STATUS_ignore, ierr)
      CALL MPI_FILE_WRITE_all(fh, dissip   ,1, datatype, MPI_STATUS_ignore, ierr)

      call MPI_File_close(fh,ierr)
      call mpi_reduce(dissum  ,   disum,   1,MPI_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      call mpi_reduce(psum/(3.1415926*r1**2)   ,   psum,   1,MPI_double_precision,mpi_sum,0,MPI_COMM_WORLD,ierr)
      psum=(psum-y1)/2/r1
      psum2=(psum2-y1)/y1
      if(rank == 0) Then
         open(200, file = 'height.dat', access = 'append', action = 'write')
         write(200,*)psum,psum2,disum
         close(200)
         print*
         print*, 'The result '//num(1:4)//' is written'
         print*
         print*,'========================================'

      endif
   
      call MPI_BARRIER(CART_COMM,ierr)
   endsubroutine MPIoutput
   
   !------------------------------------
   subroutine Gradient(c,dcx,dcy)
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::c
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::dcx
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::dcy
   integer i,j,iq,ind

   
   call MPI_BARRIER(CART_COMM,ierr)
   if(nx.eq.local_end(1))then
      do j=local_start(2)-ghost,local_end(2)+ghost
         c(nx+1,j)=c(nx-1,j)
      enddo
   elseif(0.eq.local_start(1))then
      do j=local_start(2)-ghost,local_end(2)+ghost
         c(-1,j)=c(1,j)

      enddo
   endif

   if(ny.eq.local_end(2))then
      do i=local_start(1)-ghost,local_end(1)+ghost
         c(i,ny+1)=c(i,ny-1)
      enddo
   elseif(0.eq.local_start(2))then
      do i=local_start(1)-ghost,local_end(1)+ghost
         c(i,-1)=c(i,1)
      enddo
   endif
   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)
         
         dcx(i,j)=(c(i+1,j)-c(i-1,j))/3.d0+&
         (c(i+1,j+1)-c(i-1,j-1))/12.d0+&
         (c(i+1,j-1)-c(i-1,j+1))/12.d0

         dcy(i,j)=(c(i,j+1)-c(i,j-1))/3.d0+&
         (c(i+1,j+1)-c(i-1,j-1))/12.d0+&
         (c(i-1,j+1)-c(i+1,j-1))/12.d0 
      
      enddo
   enddo
      call MPI_BARRIER(CART_COMM,ierr)
end subroutine Gradient

subroutine div(cc,uu,vv)
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::cc
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::uu
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::vv
   integer i,j,iq,ind

   
   call MPI_BARRIER(CART_COMM,ierr)

   if(nx.eq.local_end(1))then
      do j=local_start(2)-ghost,local_end(2)+ghost
         uu(nx+1,j)=-uu(nx-1,j)
         vv(nx+1,j)=vv(nx-1,j)
      enddo
   elseif(0.eq.local_start(1))then
      do j=local_start(2)-ghost,local_end(2)+ghost
         uu(-1,j)=-uu(1,j)
         vv(-1,j)=vv(1,j)
      enddo
   endif

   if(ny.eq.local_end(2))then
      do i=local_start(1)-ghost,local_end(1)+ghost
         uu(i,ny+1)=uu(i,ny-1)
         vv(i,ny+1)=vv(i,ny-1)

      enddo
   elseif(0.eq.local_start(2))then
      do i=local_start(1)-ghost,local_end(1)+ghost
         uu(i,-1)=uu(i,1)
         vv(i,-1)=vv(i,1)
      enddo
   endif
   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)

         cc(i,j)=(uu(i+1,j)-uu(i-1,j))/3.d0+&
               (uu(i+1,j+1)-uu(i-1,j-1))/12.d0+&
               (uu(i+1,j-1)-uu(i-1,j+1))/12.d0+&
               (vv(i,j+1)-vv(i,j-1))/3.d0+&
               (vv(i+1,j+1)-vv(i-1,j-1))/12.d0+&
               (vv(i-1,j+1)-vv(i+1,j-1))/12.d0 


      enddo
   enddo

   call MPI_BARRIER(CART_COMM,ierr)
endsubroutine div

subroutine multi(m1,m2,m3)
      double precision,dimension(0:nq-1,0:nq-1)::m1,m2,m3,tmp
      integer::i,j,k

      do j=0,8
         do i=0,8
            tmp(i,j)=0.d0
            do k=0,8
                  tmp(i,j)=tmp(i,j)+m1(k,j)*m2(i,k)
            enddo
            m3(i,j)=tmp(i,j)
         enddo
      enddo

endsubroutine

subroutine Gradient5(c,dcx,dcy)
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::c
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::dcx
   double precision,dimension(local_start(1)-ghost:local_end(1)+ghost,local_start(2)-ghost:local_end(2)+ghost)::dcy
   integer i,j,iq,ind

   
   call MPI_BARRIER(CART_COMM,ierr)
   if(nx.eq.local_end(1))then
      do j=local_start(2)-ghost,local_end(2)+ghost
         c(nx+1,j)=c(nx-1,j)
      enddo
   elseif(0.eq.local_start(1))then
      do j=local_start(2)-ghost,local_end(2)+ghost
         c(-1,j)=c(1,j)

      enddo
   endif

   if(ny.eq.local_end(2))then
      do i=local_start(1)-ghost,local_end(1)+ghost
         c(i,ny+1)=c(i,ny-1)
      enddo
   elseif(0.eq.local_start(2))then
      do i=local_start(1)-ghost,local_end(1)+ghost
         c(i,-1)=c(i,1)
      enddo
   endif
   do j=local_start(2),local_end(2)
      do i=local_start(1),local_end(1)
         
         dcx(i,j)=(c(i+1,j)-c(i-1,j))/2.d0

         dcy(i,j)=(c(i,j+1)-c(i,j-1))/2.d0
      
      enddo
   enddo
      call MPI_BARRIER(CART_COMM,ierr)
end subroutine Gradient5

endmodule lbm
