function [r]=q1u(func_num,ia,ib)
if(func_num==1)
    x0=0;
    x1=pi;
    y0=0;
    y1=pi;
    mui=@(a,b) 0.5*(a+sqrtm(a^2-4*b));
    f=@(x,y,a,b) (4+2*a+b)*sin(x)*sin(y);
    
    u_xy=@(x,y) sin(x)*sin(y);
    g1=@(x,y) sin(x)*sin(y);
    g2=@(x,y) -2*sin(x)*sin(y);
    
    v=@(x,y,a,b) (1/mui(a,b))*g2(x,y)-g1(x,y);
    Up=@(x,y,a,b) v(x,pi,a,b);
    Down=@(x,y,a,b) v(x,0,a,b);
    Zuo=@(x,y,a,b) v(0,y,a,b);
    You=@(x,y,a,b) v(pi,y,a,b);
end
if(func_num==2)
    x0=0;
    x1=pi;
    y0=0;
    y1=pi;
    mui=@(a,b) 0.5*(a+sqrtm(a^2-4*b));
    f=@(x,y,a,b) (1+a+b)*(x*sin(y)+y*sin(x));
    
    u_xy=@(x,y) x*sin(y)+y*sin(x);
    g1=@(x,y) x*sin(y)+y*sin(x);
    g2=@(x,y) -x*sin(y)-y*sin(x);
    
    v=@(x,y,a,b) (1/mui(a,b))*g2(x,y)-g1(x,y);
    Up=@(x,y,a,b) v(x,pi,a,b);
    Down=@(x,y,a,b) v(x,0,a,b);
    Zuo=@(x,y,a,b) v(0,y,a,b);
    You=@(x,y,a,b) v(pi,y,a,b);
end

%% initiate the parameter
a=ia;
b=ib;
N=input('������ȷַ���:');
h=(x1-x0)/N;
alpha1=-(1/(h^2))*mui(a,b);
alpha2=-(1/(h^2))*mui(a,b);
alpha3=-(1/(h^2))*mui(a,b);
alpha4=-(1/(h^2))*mui(a,b);

alpha0=alpha1+alpha2+alpha3+alpha4-b;

F=zeros((N-1)^2,1); %AU=F
A=zeros(N-1,N-1);   %AΪ(n-1)*(n-1)�׷���
%% initiate the matrix-A
%�����A��ֵ
for i=1:(N-1)^2  %�Խ���ȫΪ4
    A(i,i)=alpha0;
end
for i=1:N-1 %�Ծ���B��ֵ-1
    for j=2+(i-1)*(N-1):(N-1)+(i-1)*(N-1)
        A(j,j-1)=-alpha1;
    end
end
for i=1:N-1 %�Ծ���B��ֵ-1
    for j=2+(i-1)*(N-1):(N-1)+(i-1)*(N-1)
        A(j-1,j)=-alpha3;
    end
end
for i=N:(N-1)^2 %�Ծ���-I��ֵ-1
    A(i,i-N+1)=-alpha2;
end
for i=N:(N-1)^2 %�Ծ���-I��ֵ-1
    A(i-N+1,i)=-alpha4;
end

%% initiate the matrix-F
for i=1:N-3 %����߽��޹ص�fij��ֵ ,i means y axis, j means x axis
    for j=2:N-2
        F((N-1)*i+j)=f(j*h,(i+1)*h,a,b);
    end
end

%��ʼ��4���ǵ�
F(1)=f(x0+h,y0+h,a,b)+alpha1*Zuo(x0,y0+h,a,b)+alpha1*Down(x0+h,y0,a,b); %���½� (x0+h,y0+h)
F(N-1)=f(x1-h,y0+h,a,b)+alpha1*You(x1,y0+h,a,b)+alpha1*Down(x1-h,y0,a,b); %���½� (x1-h,y0+h)
F((N-1)*(N-2)+1)=f(x0+h,y1-h,a,b)+alpha1*Zuo(x0,y1-h,a,b)+alpha1*Up(x0+h,y1,a,b) ; %���Ͻ� (x0+h,y1-h)    
F((N-1)*(N-1))=f(x1-h,y1-h,a,b)+alpha1*Up(x1-h,y1,a,b)+alpha1*You(x1,y1-h,a,b); %���Ͻ� (x1-h,y1-h)


for i=1:N-3 %��ʼ�������µ�fij (x0+(i+1)*h,y0+h)
    F(i+1)=f(x0+(i+1)*h,y0+h,a,b)+alpha1*Down(x0+(i+1)*h,y0,a,b);
end
for i=1:N-3 %��ʼ�������ϵ�fij (x0+(i+1)*h,y1-h)
    F((N-1)*(N-2)+1+i)=f(x0+(i+1)*h,y1-h,a,b)+alpha1*Up(x0+(i+1)*h,y1,a,b);
end
for i=1:N-3 %��ʼ���������fij (x0+h,y0+(i+1)*h)
    F((N-1)*i+1)=f(x0+h,y0+(i+1)*h,a,b)+alpha1*Zuo(x0,y0+(i+1)*h,a,b);
end
for i=1:N-3 %��ʼ�������ҵ�fij (x1-h,y0+(i+1)*h)
    F((N-1)*(i+1))=f(x1-h,y0+(i+1)*h,a,b)+alpha1*You(x1,y0+(i+1)*h,a,b);
end

%% compute
u=A\F;

U1=zeros(N-1,N-1);%U1���ڴ������u��Ԫ�ض�ά���������

for i=1:N-1%u��һά���ݷŵ���ά����U1��
    for j=1:N-1
        U1(i,j)=u((N-1-i)*(N-1)+j);
    end
end

U1

%% Display the v
X=linspace(x0,x1,N-1);
Y=linspace(y0,y1,N-1);


%% Begin section 2, refresh the parameter
alpha1=-(1/(h^2))*(1/mui(a,b));
alpha2=-(1/(h^2))*(1/mui(a,b));
alpha3=-(1/(h^2))*(1/mui(a,b));
alpha4=-(1/(h^2))*(1/mui(a,b));

alpha0=alpha1+alpha2+alpha3+alpha4-1;

Up=@(x,y,a,b) g1(x,pi);
Down=@(x,y,a,b) g1(x,0);
Zuo=@(x,y,a,b) g1(0,y);
You=@(x,y,a,b) g1(pi,y);
if(func_num == 2)
    Up=@(x,y,a,b) pi*sin(x);
    Down=@(x,y,a,b) 0;
    Zuo=@(x,y,a,b) 0;
    You=@(x,y,a,b) pi*sin(y);
end
F=zeros((N-1)^2,1); %AU=F
A=zeros(N-1,N-1);   %AΪ(n-1)*(n-1)�׷���
%% Section 2 -- initiate the matrix-A
%�����A��ֵ
for i=1:(N-1)^2  
    A(i,i)=alpha0;
end
for i=1:N-1 %�Ծ���B��ֵ-1
    for j=2+(i-1)*(N-1):(N-1)+(i-1)*(N-1)
        A(j,j-1)=-alpha1;
    end
end
for i=1:N-1 %�Ծ���B��ֵ-1
    for j=2+(i-1)*(N-1):(N-1)+(i-1)*(N-1)
        A(j-1,j)=-alpha3;
    end
end
for i=N:(N-1)^2 %�Ծ���-I��ֵ-1
    A(i,i-N+1)=-alpha2;
end
for i=N:(N-1)^2 %�Ծ���-I��ֵ-1
    A(i-N+1,i)=-alpha4;
end


%% Section2-- initiate the matrix-F
for i=1:N-3 %����߽��޹ص�fij��ֵ ,i means y axis, j means x axis
    for j=2:N-2
        F((N-1)*i+j)=U1(N-1-i,j);
    end
end

%��ʼ��4���ǵ�
F(1)=U1(N-1,1)+alpha1*Zuo(x0,y0+h,a,b)+alpha1*Down(x0+h,y0,a,b); %���½� (x0+h,y0+h)
F(N-1)=U1(N-1,N-1)+alpha1*You(x1,y0+h,a,b)+alpha1*Down(x1-h,y0,a,b); %���½� (x1-h,y0+h)
F((N-1)*(N-2)+1)=U1(1,1)+alpha1*Zuo(x0,y1-h,a,b)+alpha1*Up(x0+h,y1,a,b) ; %���Ͻ� (x0+h,y1-h)    
F((N-1)*(N-1))=U1(1,N-1)+alpha1*Up(x1-h,y1,a,b)+alpha1*You(x1,y1-h,a,b); %���Ͻ� (x1-h,y1-h)

for i=1:N-3 %��ʼ�������µ�fij (x0+(i+1)*h,y0+h)
    F(i+1)=U1(N-1,i+1)+alpha1*Down(x0+(i+1)*h,y0,a,b);
end
for i=1:N-3 %��ʼ�������ϵ�fij (x0+(i+1)*h,y1-h)
    F((N-1)*(N-2)+1+i)=U1(1,i+1)+alpha1*Up(x0+(i+1)*h,y1,a,b);
end
for i=1:N-3 %��ʼ���������fij (x0+h,y0+(i+1)*h)
    F((N-1)*i+1)=U1(i+1,1)+alpha1*Zuo(x0,y0+(i+1)*h,a,b);
end
for i=1:N-3 %��ʼ�������ҵ�fij (x1-h,y0+(i+1)*h)
    F((N-1)*(i+1))=U1(i+1,N-1)+alpha1*You(x1,y0+(i+1)*h,a,b);
end
% F
%% compute
u=A\F;

U3=zeros(N-1,N-1);%U1���ڴ������u��Ԫ�ض�ά���������
U2=zeros(N-1,N-1);%U2���ڴ�����

for i=1:N-1%u��һά���ݷŵ���ά����U1��
    for j=1:N-1
        U3(i,j)=u((N-1-i)*(N-1)+j);
    end
end
for i=1:N-1%�������Ľڵ㺯��ֵ�ŵ�U2��
    for j=1:N-1
        U2(i,j)=u_xy(j*h,(N-i)*h);
    end
end

U3
U2
%% Display the v
set(gcf,'unit','centimeters','position',[3,3,25,25]);
X=linspace(x0,x1,N-1);
Y=linspace(y0,y1,N-1);
subplot(2,1,1)
mesh(X,Y,U3);
xlabel('x');
ylabel('y');
title(['��ֵ��ͼ��  a=',num2str(a),'  b=',num2str(b),'  h=',num2str(h),'  �ȷַ���',num2str(N)]);
subplot(2,1,2)
mesh(X,Y,U2);
xlabel('x');
ylabel('y');
title(['��ȷ��ͼ��  a=',num2str(a),'  b=',num2str(b),'  h=',num2str(h),'  �ȷַ���',num2str(N)]);


r=sum(sum(abs(U3-U2)))/((N-1)*(N-1));










