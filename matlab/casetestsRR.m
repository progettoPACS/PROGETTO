clc
close all
clear all

S=sym('S');
Bx=sym('Bx');
By=sym('By');
Bz=sym('Bz');
Mu=sym('Mu');
Chi=sym('Chi');
gamma=Chi/Mu;
x=sym('x');
y=sym('y');
z=sym('z');
Ly=sym('Ly');
Lx=sym('Lx');
Lz=sym('Lz');

g=sym('0');
b= 1.0/( 1.0 + x*z);
fstar=(-2.0/3.0/Ly*b*y^3+b*y^2+g)*70;
%ues=sym('z*(Lz-z)*exp ( -gamma / Ly * (y -  Ly/2 ) ^ 2 + fstar )');
ues=10^5*z*(Lz-z)*exp( - gamma/Ly * (y - Ly/2 ) ^ 2 + fstar  )*(Lx-x)^2;
pretty(ues)
dxues=diff(ues,'x');
dyues=diff(ues,'y');
dzues=diff(ues,'z');
dxxues=diff(dxues,'x');
dyyues=diff(dyues,'y');
dzzues=diff(dzues,'z');
LAPLues=dxxues+dyyues+dzzues;

neumann_out  = subs(dxues,'x',Lx)
dirichlet_in = subs(ues,'x',0)
robin_left   = -Mu*subs(dyues,'y',0)  +  Chi*subs(ues,'y',0)
robin_right = +Mu*subs(dyues,'y',Ly) +  Chi*subs(ues,'y',Ly)



f=(-Mu*LAPLues+Bx*dxues+By*dyues+Bz*dzues+S*ues);

s=2;
bx=5;
by=1;
bz=2;
mu=1;
chi=3;
ly=0.1;
lz=0.1;
lx=0.1;

f=subs(subs(subs(subs(subs(f,'S',s),'Bx',bx),'By',by),'Mu',mu),'Bz',bz)

[x,y,z] = meshgrid(0:lx/100:lx,0:ly/100:ly,0:lz/100:lz);
u=matlabFunction(ues);
v=u(chi,lx,ly,lz,mu,x,y,z);

xslice = linspace(0,lx,10); yslice = []; zslice = [];
contourslice(x,y,z,v,xslice,yslice,zslice);daspect([1,1,1])
camva(6); 
camproj perspective;
campos([0.6,-0.9,0.45])
set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
set(gca,'Color','black','XColor','white', ...
	'YColor','white','ZColor','white')
box on
