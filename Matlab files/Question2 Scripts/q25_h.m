function [r]=q25_h(v,w,u,D) %最简隐格式 暂时默认v=w
L=1; %取L=1
T=1; %取模拟的时间为1
N=40;
h=L/N;
t=h/w;
M=T/t;

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
ukj(M+1,:)=Ck;
for i=1:N-1
    B(i,i+1)=a-b;
    B(i+1,i)=-a-b;
end
B=(1+2*b)*I+B;B(1,2)=-2*b;
B
Ck
%% calculate all layers
for i=1:M 
    if(length(Ck)>1)
    Ck1=gs(B,Ck);
    Ck=Ck1(2:length(Ck1)-1);
    B=B(1:N-i-1,1:N-i-1);
    ukj(M+1-i,:)=[zeros(i,1);Ck;zeros(i,1)];
    end
end
%% display result
ukj=[ukj,zeros(M+1,1)];
ukj
set(gcf,'unit','centimeters','position',[3,3,25,25]);
X=linspace(0,1,N+1);
Y=linspace(0,1,M+1);
subplot(2,1,1)
mesh(X,Y,flipud(ukj));
xlabel('x');
ylabel('t');
zlabel('溢油浓度');
title(['最简隐格式  ','v=',num2str(v),'  w=',num2str(w),'  h=',num2str(h),'  \tau=',num2str(t),'  u=',num2str(u),'  D=',num2str(D),'  a=',num2str(a),'  b=',num2str(b)]);
for i=1:M+1
subplot(2,1,2)
plot(X,ukj(i,:));
hold on
end
xlabel('x');
ylabel('溢油浓度');
title(['最简隐格式  ','v=',num2str(v),'  w=',num2str(w),'  h=',num2str(h),'  \tau=',num2str(t),'  u=',num2str(u),'  D=',num2str(D),'  a=',num2str(a),'  b=',num2str(b)]);