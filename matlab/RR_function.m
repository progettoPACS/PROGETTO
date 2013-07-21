L=0.1;
mu=1;
chi=0.1;
h=1000;

f=@(x) 2*mu*x/L+tan(x).*(chi-mu*mu*x.*x/L^2/chi);
x=linspace(0,3*pi,1e5);
plot(x,f(x))
ylim([-h,h])
grid on
hold on
