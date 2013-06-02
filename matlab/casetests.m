Ly=2;
Lz=2;
Lx=3;

ues=sym('y*(Ly-y)*z*(Lz-z)*exp((x-Lx)^2*y*z*2)');
dxues=diff(ues,'x');
dyues=diff(ues,'y');
dzues=diff(ues,'z');
dxxues=diff(dxues,'x');
dyyues=diff(dyues,'y');
dzzues=diff(dzues,'z');
LAPLues=dxxues+dyyues+dzzues;

neumann_out=subs(subs(dxues,'Lx',Lx),'x',Lx)

dirichlet_in=subs(subs(dxues,'Lx',Lx),'x',0)

S=sym('S');
Bx=sym('Bx');
By=sym('By');
Bz=sym('Bz');
Mu=sym('Mu');

f=(-Mu*LAPLues+Bx*dxues+By*dyues+Bz*dzues+S*ues);

s=0;
bx=0;
by=0;
bz=0;
mu=1;

f=subs(subs(subs(subs(subs(f,'S',s),'Bx',bx),'By',by),'Mu',mu),'Bz',bz)