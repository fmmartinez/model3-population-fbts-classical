program popmodel3
use ifport
implicit none

integer,parameter :: nmap = 3

real(8),parameter :: pi=3.1415926535d0, oc=37.7d0, kc=sqrt(10d0)*oc

complex(8) :: coeff,fact1,fact2,fact3
complex(8),dimension(:),allocatable :: pop,pop1,pop2,pop3

integer :: i,j,ng,nb,nd,basispc,stp
integer :: np,nosc,nmcs,nmds,seed_dimension,bath,init,mcs,it,is,ib
!integer,dimension(:),allocatable :: seed

real(8) :: delta,ome_max,dt,lumda_d,eg,eb,ed,mu,e0,beta,time_j,taw_j,omega_j,check
real(8) :: dt2,uj,qbeta,lambdacheck,a1,a2,et,gaussian
real(8) :: pc,qc,av1,av2,fc
real(8),dimension(1:3) :: rm,pm,rn,pn
real(8),dimension(:),allocatable :: ome,c2,kosc,x,p,fx
real(8),dimension(:,:),allocatable :: hm

call iniconc()

call srand(seed_dimension)

allocate(ome(1:nosc),c2(1:nosc),kosc(1:nosc))
allocate(hm(1:nmap,1:nmap))

call iniconq_d()

allocate(pop(1:nmds+1),pop1(1:nmds+1),pop2(1:nmds+1),pop3(1:nmds+1))
allocate(x(1:nosc),p(1:nosc),fx(1:nosc))

pop  = 0d0
pop1 = 0d0
pop2 = 0d0
pop3 = 0d0

dt  = 2d0*pi*dt
dt2 = 0.5d0*dt

MC: do mcs = 1, nmcs
   do is=1,nosc
      uj = 0.5d0*beta*dsqrt(kosc(is))
      
      qbeta = beta/(uj/tanh(uj))
      
      p(is) = gauss_noise2()/dsqrt(qbeta)
      x(is) = gauss_noise2()/dsqrt(qbeta*kosc(is))
      
      if(bath == 1) x(is)=x(is)+c2(is)/kosc(is)
   end do
   
   uj = 0.5d0*beta*dsqrt(oc**2)
   qbeta = beta/(uj/tanh(uj))
   pc = gauss_noise2()/dsqrt(qbeta)
   qc = gauss_noise2()/dsqrt(qbeta*oc**2)

   if (init == 3) then
      do i = 1, nmap
         rm(i) = gauss_noise2()/sqrt(2d0)
         pm(i) = gauss_noise2()/sqrt(2d0)
         rn(i) = gauss_noise2()/sqrt(2d0)
         pn(i) = gauss_noise2()/sqrt(2d0)
      end do
   else
      rm = 0d0
      pm = 0d0
      rn = 0d0
      pn = 0d0
   end if
   
   coeff = 0.25d0*(cmplx(rm(1),-pm(1))*cmplx(rn(1),pn(1)))
   
   call get_force_bath(nmap,kosc,x,c2,rm,pm,rn,pn,fx)
   call get_force_coupledosc(nmap,oc,qc,kc,rm,pm,rn,pn,fc)

   ib = 1
   
   call get_facts_pop(ng,nb,nd,coeff,rm,pm,rn,pn,fact1,fact2,fact3)
   
   pop(ib)  = pop(ib)  + (fact1+fact2+fact3)
   pop1(ib) = pop1(ib) + (fact1)
   pop2(ib) = pop2(ib) + (fact2)
   pop3(ib) = pop3(ib) + (fact3)
   
   a1 = 0d0
   a2 = 0d0
   do is=1,nosc
      a1 = a1 + 2.d0*c2(is)**2/ome(is)**2
      a2 = a2 + 2.d0*c2(is)*x(is)
   end do

   av1 = 2d0*kc**2/oc**2
   av2 = 2d0*kc*qc
   
   MD: do it = 1, nmds
      gaussian=sqrt(4.d0*log(2.d0)/(pi*taw_j**2))*exp(-4.d0*log(2.d0)*((it-0.5d0)*dt-time_j)**2/(taw_j**2))
      et = gaussian*e0*cos(omega_j*((it-0.5d0)*dt-time_j))
   
      do is = 1, nosc
         p(is) = p(is) + dt2*fx(is)
      end do
      pc = pc + dt2*fc
      
      call get_hm(delta,mu,et,a1,a2,av1,av2,pc,oc,qc,hm)

      call evolve_pm(nmap,dt2,hm,rm,pm)
      call evolve_pm(nmap,dt2,hm,rn,pn)

      do is = 1, nosc
         x(is) = x(is) + dt*p(is)
      end do
      qc = qc + dt*pc

      a2=0.d0
      do is = 1, nosc
          a2 = a2 + 2.d0*c2(is)*x(is)
      end do
      av2 = 2d0*kc*qc

      call update_hm(a1,a2,av1,av2,pc,oc,qc,hm)

      call evolve_rm(nmap,dt,hm,pm,rm)
      call evolve_rm(nmap,dt,hm,pn,rn)

      call evolve_pm(nmap,dt2,hm,rm,pm)
      call evolve_pm(nmap,dt2,hm,rn,pn)

      call get_force_bath(nmap,kosc,x,c2,rm,pm,rn,pn,fx)
      call get_force_coupledosc(nmap,oc,qc,kc,rm,pm,rn,pn,fc)
      
      do is = 1, nosc
         p(is) = p(is) + dt2*fx(is)
      end do
      pc = pc + dt2*fc

      ib = it + 1
      
      call get_facts_pop(ng,nb,nd,coeff,rm,pm,rn,pn,fact1,fact2,fact3)
      
      pop(ib)  = pop(ib)  + (fact1+fact2+fact3)
      pop1(ib) = pop1(ib) + (fact1)
      pop2(ib) = pop2(ib) + (fact2)
      pop3(ib) = pop3(ib) + (fact3)
   end do MD

   if (mod(mcs,1000) == 0) then
      open(444,file='temp.out')
      do i = 1, nmds+1
         write(444,'(i10,4f20.9)') i-1,dble(pop1(i)/pop(i)),dble(pop2(i)/pop(i)),dble(pop3(i)/pop(i)),dble(pop(i))
      end do
      close(444)
   end if
end do MC

do ib = 1, nmds+1
   pop1(ib) = pop1(ib)/pop(ib)
   pop2(ib) = pop2(ib)/pop(ib)
   pop3(ib) = pop3(ib)/pop(ib)
   write(333,'(i10,4f15.9)') ib-1, dble(pop1(ib)),dble(pop2(ib)),dble(pop3(ib)),dble(pop(ib))!/dnmcs
end do

deallocate(ome,c2,kosc)
deallocate(pop,pop1,pop2,pop3)
deallocate(x,p,fx)

contains

subroutine evolve_rm(nmap,dt,hm,pm,rm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt
real(8),dimension(:),intent(in) :: pm
real(8),dimension(:),intent(inout) :: rm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      rm(i) = rm(i) + dt*hm(i,j)*pm(j)
   end do
end do

end subroutine evolve_rm

subroutine evolve_pm(nmap,dt2,hm,rm,pm)
implicit none

integer :: i, j
integer,intent(in) :: nmap

real(8),intent(in) :: dt2
real(8),dimension(:),intent(in) :: rm
real(8),dimension(:),intent(inout) :: pm
real(8),dimension(:,:),intent(in) :: hm

do i = 1, nmap
   do j = 1, nmap
      pm(i) = pm(i) - dt2*hm(i,j)*rm(j)
   end do
end do

end subroutine evolve_pm

subroutine get_facts_pop(vg,vb,vd,coeff,rm,pm,rn,pn,fact1,fact2,fact3)
implicit none

integer :: a,b,i
integer,intent(in) :: vg,vb,vd

complex(8),intent(in) :: coeff
complex(8),intent(out) :: fact1,fact2,fact3

real(8) :: temp
real(8),dimension(:),intent(in) :: rm,pm,rn,pn

fact1 = 0d0
do a = 1, vg
   fact1 = fact1 + coeff*(cmplx(rm(a),pm(a))*cmplx(rn(a),-pn(a)))
end do

fact2 = 0d0
do a = vg+1, vg+vb
   fact2 = fact2 + coeff*(cmplx(rm(a),pm(a))*cmplx(rn(a),-pn(a)))
end do

fact3 = 0d0
do a = vg+vb+1, vg+vb+vd
   fact3 = fact3 + coeff*(cmplx(rm(a),pm(a))*cmplx(rn(a),-pn(a)))
end do
end subroutine get_facts_pop

function kronecker_delta(i,j) result (d)
implicit none

integer,intent(in) :: i, j

real(8) :: d

if (i == j) then
   d = 1d0
else
   d = 0d0
end if

end function kronecker_delta

subroutine get_force_bath(nmap,kosc,x,c2,rm,pm,rn,pn,f)
implicit none

integer :: a,b,i,j,n
integer,intent(in) :: nmap

real(8) :: trace
real(8),dimension(:),intent(in) :: kosc,x,c2,rm,pm,rn,pn
real(8),dimension(:),intent(out) :: f
real(8),dimension(1:3,1:3) :: dh

n = size(x)

f = 0d0
do j = 1, n
   f(j) = -kosc(j)*x(j)
   
   !original dh/dR 
!   dh = 0d0 
!   dh(3,3)= c2(j)
   
   trace = -2d0*c2(j)
!   do a = 1, nmap
!      trace = trace + dh(a,a)
!   end do
   !traceless
!   do a = 1, nmap
!      dh(a,a) = dh(a,a) - trace/nmap
!   end do
   
!   do a = 1, nmap
!      b = a
      f(j) = f(j) + (-c2(j))* &
                  (rm(3)*rm(3) + pm(3)*pm(3) + rn(3)*rn(3) + pn(3)*pn(3) )
!   end do

   f(j) = f(j) - trace
end do

end subroutine get_force_bath

subroutine get_force_coupledosc(nmap,oc,qc,kc,rm,pm,rn,pn,f)
implicit none

integer :: a,b
integer,intent(in) :: nmap

real(8) :: trace
real(8),intent(in) :: oc,qc,kc
real(8),intent(out) :: f
real(8),dimension(:),intent(in) :: rm,pm,rn,pn
real(8),dimension(1:3,1:3) :: dh

f = -oc**2*qc

!original dh/dq
dh = 0d0
dh(2,2) = -kc
dh(3,3) = -0.25d0*kc

trace = 0d0
do a = 1, nmap
   trace = trace + dh(a,a)
end do

do a = 1, nmap
   dh(a,a) = dh(a,a) - trace/nmap
end do

do a = 1, nmap
   b = a
   f = f - 0.5d0*dh(a,a)*(rm(a)*rm(b) + pm(a)*pm(b) + rn(a)*rn(b) + pn(a)*pn(b))
end do

end subroutine get_force_coupledosc

subroutine update_hm(a1,a2,av1,av2,pc,oc,qc,hm)
implicit none

real(8),parameter :: eg=0, eb=240, ed=240

real(8) :: ev
real(8),intent(in) :: a1,a2,av1,av2,pc,oc,qc
real(8),dimension(:,:),intent(inout) :: hm

ev = 0.5d0*(pc**2 + (oc*qc)**2)

!only diagonal part is updated
!1 x 1
hm(1,1) = 0d0
hm(1,1) = eg + ev
!2 x 2
hm(2,2) = 0d0
hm(2,2) = eb + ev + av1 - av2
!3 x 3
hm(3,3) = 0d0
hm(3,3) = ed + ev + a1 + a2 + 0.25d0*av1 - 0.5d0*av2

end subroutine update_hm

subroutine make_hmtraceless(nmap,hm)
implicit none

integer :: i
integer,intent(in) :: nmap

real(8) :: trace
real(8),dimension(:,:),intent(inout) :: hm

trace = 0d0
do i = 1, nmap
   trace = trace + hm(i,i)
end do

do i = 1, nmap
   hm(i,i) = hm(i,i) - trace/nmap
end do

end subroutine make_hmtraceless


subroutine get_hm(delta,mu,et,a1,a2,av1,av2,pc,oc,qc,hm)
implicit none

real(8),parameter :: eg=0, eb=240, ed=240

real(8) :: ev
real(8),intent(in) :: delta,mu,et,a1,a2,pc,oc,qc,av1,av2
real(8),dimension(:,:),intent(out) :: hm

ev = 0.5d0*(pc**2 + (oc*qc)**2)

hm = 0d0
!1 x 1
hm(1,1) = eg + ev
!1 x 2
hm(1,2) = -mu*et
!1 x 3
!this part is zero
!2 x 1
hm(2,1) = hm(1,2)
!2 x 2
hm(2,2) = eb + ev + av1 - av2
!2 x 3
hm(2,3) = delta
!3 x 1
!this part is zero
!3 x 2
hm(3,2) = hm(2,3)
!3 x 3
hm(3,3) = ed + ev + a1 + a2 + 0.25d0*av1 - 0.5d0*av2

end subroutine get_hm

subroutine iniconc()
implicit none

integer :: i

open (666,file='md.in')
read(666,*)
read(666,*) np,delta,nosc,ome_max
read(666,*)
read(666,*) nmcs,nmds,seed_dimension,dt,lumda_d
read(666,*)
read(666,*) eg,eb,ed,mu,e0,beta
read(666,*)
read(666,*) time_j,taw_j,omega_j 
read(666,*)
read(666,*) bath,init
read(666,*)
read(666,*) basispc,ng,nb,nd
close(666)

!call random_seed(size=seed_dimension)
!allocate (seed(seed_dimension))
!do i=1,seed_dimension
!  seed(i) = 3*2**i-1
!enddo
!call random_seed(put=seed)

end subroutine iniconc

subroutine iniconq_d()
implicit none

integer :: i

check=0
do i=1,nosc
!  ome(i)=tan(pi*i/(4*nosc))
   ome(i)=tan(i*datan(ome_max)/nosc)
!  ome(i)=i**2*ome_max/nosc**2
!  rho(i)=sqrt(ome(i)*ome_max)/(1+ome(i)**2)
!  c2(i)=0.5d0*ome(i)*dsqrt(lumda_d/nosc)
   c2(i)=ome(i)*dsqrt(datan(ome_max)*lumda_d/(pi*nosc))
!  c2(i)=ome(i)*sqrt(lumda_d*rho(i)*2/(pi*nosc))
   check=check+c2(i)**2/ome(i)**2
   kosc(i)=ome(i)**2 
end do
write(6,*) check

end subroutine iniconq_d

function gauss_noise2() result(g)
implicit none

real(8),parameter :: pi2=2.0*3.141592654

real(8) :: g,z1,z2

!call random_number(z1)
!call random_number(z2)
z1 = rand()
z2 = rand()
g = sqrt(-2.d0*log(z1))*cos(pi2*z2)

end function gauss_noise2

end program popmodel3
