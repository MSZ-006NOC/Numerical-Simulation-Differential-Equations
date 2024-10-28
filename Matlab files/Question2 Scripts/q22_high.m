function [r]=q22_h(u,D) %高阶精度
L=1; %取L=1
T=10; %取模拟的时间为1
N=10;
M=200;
h=L/N;
t=T/M;
a=(u*t)/(2*h);
b=(D*t)/(h^2);

Ck=zeros(N,1);
Ck1=zeros(N,1);
ukj=zeros(M+1,N);
B=zeros(N,N);
I=eye(N);
%% initiate the Ck and B
for i=0:N-1
    Ck(i+1,1)=sin(pi*i*h);
end
Ck(1,1)=0;
ukj(M+1,:)=Ck;
for i=1:N-1
    B(i,i+1)=a-b;
    B(i+1,i)=-a-b;
end
B=(1+2*b)*I+B;B(1,2)=-4;B(1,1)=3;B(1,3)=1;
B
%% calculate all layers
for i=1:M %T_a*Ck1=T_b*Ck
    Ck1=gs(B,Ck);
    Ck=Ck1;
    ukj(M+1-i,:)=Ck;
    Ck(1,1)=0;
end
%% display result
ukj=[ukj,zeros(M+1,1)];
ukj
X=linspace(0,1,N+1);
Y=linspace(0,10,M+1);
subplot(1,1,1)
mesh(X,Y,flipud(ukj));
xlabel('x');
ylabel('t');
zlabel('溢油浓度');
title(['高阶精度格式  ','h=',num2str(h),'  \tau=',num2str(t),'  u=',num2str(u),'  D=',num2str(D),'  a=',num2str(a),'  b=',num2str(b)]);