function [r]=q21_h(u,D) %最简隐格式
L=1; %取L=1
T=1; %取模拟的时间为1
N=5;
M=50;
h=L/N;
t=T/M;
a=(u*t)/(2*h);
b=(D*t)/(h^2);

Ck=zeros(N-1,1);
Ck1=zeros(N-1,1);
ukj=zeros(M+1,N-1);
B=zeros(N-1,N-1);
I=eye(N-1);
%% initiate the Ck and B
for i=1:N-1
    Ck(i,1)=sin(pi*i*h);
end
ukj(M+1,:)=Ck;
for i=1:N-2
    B(i,i+1)=a-b;
    B(i+1,i)=-a-b;
end
%% calculate all layers
T_a=((1+2*b)*I+B);
for i=1:M %T_a*Ck1=T_b*Ck
    Ck1=gs(T_a,Ck);
    Ck=Ck1;
    ukj(M+1-i,:)=Ck;
end
%% display result
ukj=[zeros(M+1,1),ukj,zeros(M+1,1)];
ukj
X=linspace(0,1,N+1);
Y=linspace(0,1,M+1);
subplot(1,1,1)
mesh(X,Y,flipud(ukj));
xlabel('x');
ylabel('t');
zlabel('溢油浓度');
axis([0,1,0,1,0,1]);
title(['最简隐格式  ','h=',num2str(h),'  \tau=',num2str(t),'  u=',num2str(u),'  D=',num2str(D),'  a=',num2str(a),'  b=',num2str(b)]);

