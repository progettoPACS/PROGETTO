clc
N=10;
m=8;

subplot(1,2,1);
B=ones(m,m*3);
A=[ones(m,m*2),zeros(m,m*(N-2))];
for k=2:N-1;
    A=[A;zeros(m,m*(k-2)),B,zeros(m,m*(N-k-1))];
end
A=[A;zeros(m,m*(N-2)),ones(m,m*2)];
T=diag(ones(N,1)) + diag(ones(N-1,1),1) + diag(ones(N-1,1),-1);
B=[];
v=[];
for k=1:m;
    v=[v,T];
end
for k=1:m
    B=[B;v];
end
spy(B);
subplot(1,2,2);
spy(A);


