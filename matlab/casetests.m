Ly=0.1;
Lz=0.1;
Lx=0.2;

ues=sym('10^7*y*(Ly-y)*z*(Lz-z)*(x-Lx)^2*exp((x-Lx)^2*y*z*2)');

dxues=diff(ues,'x');
dyues=diff(ues,'y');
dzues=diff(ues,'z');
dxxues=diff(dxues,'x');
dyyues=diff(dyues,'y');
dzzues=diff(dzues,'z');
LAPLues=dxxues+dyyues+dzzues;

neumann_out=subs(subs(dxues,'Lx',Lx),'x',Lx)

dirichlet_in=subs(ues,'x',0)

S=sym('S');
Bx=sym('Bx');
By=sym('By');
Bz=sym('Bz');
Mu=sym('Mu');

f=(-Mu*LAPLues+Bx*dxues+By*dyues+Bz*dzues+S*ues);

s=3;
bx=5;
by=1;
bz=1;
mu=1;

f=subs(subs(subs(subs(subs(f,'S',s),'Bx',bx),'By',by),'Mu',mu),'Bz',bz)
collect(f,sym('exp(2*y*z*(Lx - x)^2)'))
[x,y,z] = meshgrid(0:Lx/100:Lx,0:Ly/100:Ly,0:Lz/100:Lz);
u=matlabFunction(ues);
v=u(Lx,Ly,Lz,x,y,z);

xslice = linspace(0,Lx,10); yslice = []; zslice = [];
contourslice(x,y,z,v,xslice,yslice,zslice);daspect([1,1,1])
camva(6); 
camproj perspective;
campos([0.6,-0.9,0.45])
set(gcf,'Color',[.5,.5,.5],'Renderer','zbuffer')
set(gca,'Color','black','XColor','white', ...
	'YColor','white','ZColor','white')
box on