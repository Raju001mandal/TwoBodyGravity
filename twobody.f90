program gravity2body


real::f11,f21,g11,g21,x1,y1,x10,y10,to,tf,t,dt,omegaz,x1dot,y1dot,x1dot0,y1dot0,x2,x20,y2,y20,x2dot,x2dot0,y2dot,y2dot0 ,&
      k11,k21,k31,k41,l11,l21,l31,l41,k12,k22,k32,k42,l12,l22,l32,l42,k13,k23,k33,k43,l13,l23,l33,l43 ,&
      f12,f22,g12,g22,k14,k24,k34,k44,l14,l24,l34,l44
      
      
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0
integer(8)::n,i



print*,"give the initial time and final time"
read*,t0,tf
print*,"give the number of interval"
read*,n
!print*,"give the length of the pendulum "
!read*,R
print*,"give the initial value of x1,y1,x1dot,y1dot,x2,y2,x2dot,y2dot"
read*,x10,y10,x1dot0,y1dot0,x20,y20,x2dot0,y2dot0

!==================================================================================
!Initialization
!===================================================================================
x1=xor+x10
y1=yor+y10

x1dot=x1dot0
y1dot=y1dot0
!====================================================================================

x2=xor+x20
y2=yor+y20

x2dot=x2dot0
y2dot=y2dot0

!======================================================================================


open(30,file="twobody.dat")
!open(35,file="modeuler")
dt=(tf-t0)/n

do i=0,n


!=======================================================================================================================
k11=dt*f11(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
l11=dt*g11(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
k21=dt*f11(x1+k11*0.5,x2+k11*0.5,y1+k11*0.5,y2+k11*0.5,x1dot+l11*0.5,x2dot+l11*0.5,y1dot+l11*0.5,y2dot+l11*0.5,t+dt*0.5)
l21=dt*g11(x1+k11*0.5,x2+k11*0.5,y1+k11*0.5,y2+k11*0.5,x1dot+l11*0.5,x2dot+l11*0.5,y1dot+l11*0.5,y2dot+l11*0.5,t+dt*0.5)
k31=dt*f11(x1+k21*0.5,x2+k21*0.5,y1+k21*0.5,y2+k21*0.5,x1dot+l21*0.5,x2dot+l21*0.5,y1dot+l21*0.5,y2dot+l21*0.5,t+dt*0.5)
l31=dt*g11(x1+k21*0.5,x2+k21*0.5,y1+k21*0.5,y2+k21*0.5,x1dot+l21*0.5,x2dot+l21*0.5,y1dot+l21*0.5,y2dot+l21*0.5,t+dt*0.5)
k41=dt*f11(x1+k31*0.5,x2+k31*0.5,y1+k31*0.5,y2+k31*0.5,x1dot+l31*0.5,x2dot+l31*0.5,y1dot+l31*0.5,y2dot+l31*0.5,t+dt*0.5)
l41=dt*g11(x1+k31*0.5,x2+k31*0.5,y1+k31*0.5,y2+k31*0.5,x1dot+l31*0.5,x2dot+l31*0.5,y1dot+l31*0.5,y2dot+l31*0.5,t+dt*0.5)
!========================================================================================================================
k12=dt*f21(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
l12=dt*g21(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
k22=dt*f21(x1+k12*0.5,x2+k12*0.5,y1+k12*0.5,y2+k12*0.5,x1dot+l12*0.5,x2dot+l12*0.5,y1dot+l12*0.5,y2dot+l12*0.5,t+dt*0.5)
l22=dt*g21(x1+k12*0.5,x2+k12*0.5,y1+k12*0.5,y2+k12*0.5,x1dot+l12*0.5,x2dot+l12*0.5,y1dot+l12*0.5,y2dot+l12*0.5,t+dt*0.5)
k32=dt*f21(x1+k22*0.5,x2+k22*0.5,y1+k22*0.5,y2+k22*0.5,x1dot+l22*0.5,x2dot+l22*0.5,y1dot+l22*0.5,y2dot+l22*0.5,t+dt*0.5)
l32=dt*g21(x1+k22*0.5,x2+k22*0.5,y1+k22*0.5,y2+k22*0.5,x1dot+l22*0.5,x2dot+l22*0.5,y1dot+l22*0.5,y2dot+l22*0.5,t+dt*0.5)
k42=dt*f21(x1+k32*0.5,x2+k32*0.5,y1+k32*0.5,y2+k32*0.5,x1dot+l32*0.5,x2dot+l32*0.5,y1dot+l32*0.5,y2dot+l32*0.5,t+dt*0.5)
l42=dt*g21(x1+k32*0.5,x2+k32*0.5,y1+k32*0.5,y2+k32*0.5,x1dot+l32*0.5,x2dot+l32*0.5,y1dot+l32*0.5,y2dot+l32*0.5,t+dt*0.5)
!===================================================================================================

k13=dt*f12(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
l13=dt*g12(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
k23=dt*f12(x1+k13*0.5,x2+k13*0.5,y1+k13*0.5,y2+k13*0.5,x1dot+l13*0.5,x2dot+l13*0.5,y1dot+l13*0.5,y2dot+l13*0.5,t+dt*0.5)
l23=dt*g12(x1+k13*0.5,x2+k13*0.5,y1+k13*0.5,y2+k13*0.5,x1dot+l13*0.5,x2dot+l13*0.5,y1dot+l13*0.5,y2dot+l13*0.5,t+dt*0.5)
k33=dt*f12(x1+k23*0.5,x2+k23*0.5,y1+k23*0.5,y2+k23*0.5,x1dot+l23*0.5,x2dot+l23*0.5,y1dot+l23*0.5,y2dot+l23*0.5,t+dt*0.5)
l33=dt*g12(x1+k23*0.5,x2+k23*0.5,y1+k23*0.5,y2+k23*0.5,x1dot+l23*0.5,x2dot+l23*0.5,y1dot+l23*0.5,y2dot+l23*0.5,t+dt*0.5)
k43=dt*f12(x1+k33*0.5,x2+k33*0.5,y1+k33*0.5,y2+k33*0.5,x1dot+l33*0.5,x2dot+l33*0.5,y1dot+l33*0.5,y2dot+l33*0.5,t+dt*0.5)
l43=dt*g12(x1+k33*0.5,x2+k33*0.5,y1+k33*0.5,y2+k33*0.5,x1dot+l33*0.5,x2dot+l33*0.5,y1dot+l33*0.5,y2dot+l33*0.5,t+dt*0.5)
!=================================================================================================================

k14=dt*f22(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
l14=dt*g22(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
k24=dt*f22(x1+k14*0.5,x2+k14*0.5,y1+k14*0.5,y2+k14*0.5,x1dot+l14*0.5,x2dot+l14*0.5,y1dot+l14*0.5,y2dot+l14*0.5,t+dt*0.5)
l24=dt*g22(x1+k14*0.5,x2+k14*0.5,y1+k14*0.5,y2+k14*0.5,x1dot+l14*0.5,x2dot+l14*0.5,y1dot+l14*0.5,y2dot+l14*0.5,t+dt*0.5)
k34=dt*f22(x1+k24*0.5,x2+k24*0.5,y1+k24*0.5,y2+k24*0.5,x1dot+l24*0.5,x2dot+l24*0.5,y1dot+l24*0.5,y2dot+l24*0.5,t+dt*0.5)
l34=dt*g22(x1+k24*0.5,x2+k24*0.5,y1+k24*0.5,y2+k24*0.5,x1dot+l24*0.5,x2dot+l24*0.5,y1dot+l24*0.5,y2dot+l24*0.5,t+dt*0.5)
k44=dt*f22(x1+k34*0.5,x2+k34*0.5,y1+k34*0.5,y2+k34*0.5,x1dot+l34*0.5,x2dot+l34*0.5,y1dot+l34*0.5,y2dot+l34*0.5,t+dt*0.5)
l44=dt*g22(x1+k34*0.5,x2+k34*0.5,y1+k34*0.5,y2+k34*0.5,x1dot+l34*0.5,x2dot+l34*0.5,y1dot+l34*0.5,y2dot+l34*0.5,t+dt*0.5)

!========================================================================================================================
if(mod(i,5)==0)then
write(30,*)x1,y1,x2,y2
end if
!===================================================================================================
t=i*dt+t0
!===================================================================================================
x1=x1+(k11+2*k21+2*k31+k41)/6.0
x1dot=x1dot+(l11+2*l21+2*l31+l41)/6.0

y1=y1+(k12+2*k22+2*k32+k42)/6.0
y1dot=y1dot+(l12+2*l22+2*l32+l42)/6.0

x1=xor+x1
y1 =yor+y1
!=======================================================================================
x2=x2+(k13+2*k23+2*k33+k43)/6.0
x2dot=x2dot+(l13+2*l23+2*l33+l43)/6.0

y2=y2+(k14+2*k24+2*k34+k44)/6.0
y2dot=y2dot+(l14+2*l24+2*l34+l44)/6.0

x2=xor+x2
y2=yor+y2
!=======================================================================================
!print*,x,y
!write(30,*)x,y,t

end do

!close(35)
!call system ('gnuplot -p euler_plot.plt')
!call system ('gnuplot -p modeuler_plot.plt')
end program
!---------------------------------------------------------------------
real function f11(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

f11=x1dot

end function f11
!======================================================================

real function f21(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

f21=y1dot

end function f21

!===========================================================================
real function f12(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

f12=x2dot

end function f12
!===========================================================================

real function f22(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

f22=y2dot

end function f22

!===========================================================================
real function g11(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

g11=-G*m2*(x1-x2)/((x1-x2)**2.0+(y1-y2)**2.0)**(3/2)

end function g11
!=============================================================================

real function g21(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

g21=-G*m2*(y1-y2)/((x1-x2)**2+(y1-y2)**2)**(3/2)


end function g21

!==============================================================================
real function g12(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

g12=-G*m1*(x2-x1)/((x1-x2)**2.0+(y1-y2)**2.0)**(3/2)

end function g12
!==================================================================================

real function g22(x1,x2,y1,y2,x1dot,x2dot,y1dot,y2dot,t)
implicit none

real,intent(in)::x1,x2,y1,y2,x1dot,y1dot,x2dot,y2dot,t
real,parameter::pi=acos(-1.0),xor=0.0,yor=0.0,G=1.0,m1=10,m2=1.0

g22=-G*m1*(y2-y1)/((x1-x2)**2+(y1-y2)**2)**(3/2)

end function g22

!===================================================================================

