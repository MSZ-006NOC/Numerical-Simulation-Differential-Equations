%% 问题1.1
clc
close all
clear
a = 3;%a(可调节参�?
b = 1;%b(可调节参�?
XArea = [0, pi];
TArea = [0, pi];
h = 0.05*pi;
t = 0.05*pi;
TargetFun = @(x, y) sin(x)*sin(y);
DefaultParam = 0;
mu = 0.5*(a+ (a^2 - 4*b)^0.5);
MarginFun = @(x, y, mu) -(2/mu+1)*sin(x)*sin(y);
Fun = @(i, j, h, t, a, b, XArea, TArea, DefaultParam)(4+ 2*a+ b)*sin((i-1)*h+XArea(1))*sin((j-1)*t+TArea(1));
v = CentralDifference(TargetFun, MarginFun, XArea, TArea, h, t, mu, Fun,a, b, 0, 0);
v=v';
mu = 1/mu;
MarginFun = @(x, y, mu) sin(x)*sin(y);
Fun = @(i, j, h, t, a, b, XArea, TArea, DefaultParam) DefaultParam(i, j);
u1 = CentralDifference(TargetFun, MarginFun, XArea, TArea, h, t, mu, Fun,a, 1, 1, v);
%% 问题1.2
clc
close all
clear
a = 3;%a(可调节参�?
b = 1;%b(可调节参�?
XArea = [0, pi];
TArea = [0, pi];
h = 0.05*pi;%步长 (可调节参�?
t = 0.05*pi;%时间间隔(可调节参�?
TargetFun = @(x, y) sin(x)*y+sin(y)*x;

mu = 0.5*(a+ (a^2 - 4*b)^0.5);
MarginFun = @(x, y, mu) -(1/mu+1)*(y*sin(x)+x*sin(y));
Fun = @(i, j, h, t, a, b, XArea, TArea, DefaultParam)(1+ a+ b)*(((j-1)*t+TArea(1))*sin((i-1)*h+XArea(1))+((i-1)*h+XArea(1))*sin((j-1)*t+TArea(1)));
v = CentralDifference(TargetFun, MarginFun, XArea, TArea, h, t, mu, Fun,a, b, 0, 0);

mu = 1/mu;
MarginFun = @(x, y, mu) sin(x)*y+sin(y)*x;
Fun = @(i, j, h, t, a, b, XArea, TArea, DefaultParam) DefaultParam(j, i);
u = CentralDifference(TargetFun, MarginFun, XArea, TArea, h, t, mu, Fun,a, 1, 1, v);
%% 问题2.1
clc
clear
D = 16;;%扩散系数(可调节参�?
u = 1;%对流速度(可调节参�?
h = 0.01;%步长 (可调节参�?
t = 0.01;%时间间隔(可调节参�?
XArea = [0,1];
TArea = [0,1];
TLeftFun = @(x,t) 0;
TRightFun = @(x,t) 0;
XDownFun = @(x,t) sin(pi*x);
NumericalResult=DiffusionAndConvection1(TLeftFun, TRightFun, XDownFun, XArea, TArea, D, u, h, t);
%% 问题2.2
clc
clear
D = 16;;%扩散系数(可调节参�?
u = 1;%对流速度(可调节参�?
h = 0.01;%步长 (可调节参�?
t = 0.01;%时间间隔(可调节参�?
XArea = [0,1];
TArea = [0,1];
XDownFun = @(x,t) sin(pi*x);
NumericalResult = DiffusionAndConvection2(XDownFun, XArea, TArea, D, u, h, t);
%% 问题2.3
clc
close all
clear
t = 0.001;%时间间隔(可调节参�?
u = 1; %对流速度(可调节参�?
D = 1;%扩散系数(可调节参�?
v = 10;%收油船船�?可调节参�?
XArea = [0, 1]; 
XDonwFun = @(x, t) sin(pi*x);
NumericalResult = DiffusionAndConvection3(t, u, D, v, XArea, XDonwFun);
%% 问题2.4
clc
close all
clear
t = 0.01;%时间间隔(可调节参�?
u = 10;%对流速度(可调节参�?
D = 1;%扩散系数(可调节参�?
w = 1;%围油拦移动�?�?可调节参�?
XArea = [0, 1]; 
XDonwFun = @(x, t) sin(pi*x);
NumericalResult = DiffusionAndConvection4(t, u, D, 1/w, XArea, XDonwFun);
%% Q2.5
clc
close all
clear
n=50;%时间层等分数(可调节参�?
u = 1;%对流速度(可调节参�?
D = 1;%扩散系数(可调节参�?
w = 1;%围油拦移动�?�?可调节参�?
v = 5;%收油船船�?可调节参�?
XArea = [0, 1]; 
XDownFun = @(x, t) sin(pi*x);
NumericalResult = DiffusionAndConvection5(n, u, D, w, v, XArea, XDownFun);