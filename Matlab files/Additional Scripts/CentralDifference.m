function NumericalResult = CentralDifference(TargetFun, MarginFun, XArea, TArea, h, t, mu, Fun, a, b, isPlot, DefaultParam)
%% 初始�?
X = XArea(1): h: XArea(2);
T = TArea(1): t: TArea(2);
XSize = size(X, 2);
TSize = size(T, 2);
HSize = (TSize-2)*(XSize-2);

TargetResult = zeros(TSize, XSize);
NumericalResult = zeros(TSize, XSize);
H = zeros(HSize);
g = zeros(HSize,1);
ResultVector = zeros(HSize,1);

alpha = h^2*mu;
beta = t^2*mu;
theta = -t^2*h^2*b-2*h^2*mu-2*t^2*mu;
%% 生成求解矩阵、向�?
for i = 1 : HSize
    H(i, i) = theta;
    if mod(i, XSize-2) ~= 0
        H(i+1, i) = beta;
        H(i, i+1) = beta;
    end
end
for i = 1:HSize - XSize + 2
    H(XSize+i-2, i) = alpha;
    H(i, XSize+i-2) = alpha;
end

for i = 1:TSize-2
    for j = 1:XSize-2
        g((XSize-2)*(i-1) + j) = t^2*h^2*Fun(j+1, i+1, h, t, a, b, XArea, TArea, DefaultParam);
        if i == 1       % 下半部分
            g((XSize-2)*(i-1) + j) = g((XSize-2)*(i-1) + j) - alpha*MarginFun(X(j+1), T(1), mu);
        end
        if i == TSize-2 % 上半部分
            g((XSize-2)*(i-1) + j) = g((XSize-2)*(i-1) + j) - alpha*MarginFun(X(j+1), T(TSize), mu);
        end
        if j == XSize-2 % 右半部分
            g((XSize-2)*(i-1) + j) = g((XSize-2)*(i-1) + j) - beta*MarginFun(X(XSize), T(i+1), mu);
        end
        if j == 1       % 左半部分
            g((XSize-2)*(i-1) + j) = g((XSize-2)*(i-1) + j) - beta*MarginFun(X(1), T(i+1), mu);
        end
    end
end

%% 求解方程 && 将解矩阵�?
ResultVector= inv(H)*g;
for i = 2:TSize-1
    for j = 2:XSize-1
        NumericalResult(i, j) = ResultVector((i-2)*(XSize-2) + j-1);
    end
end
for  i = 1:XSize
    NumericalResult(1, i) = TargetFun(X(i), T(1));
    NumericalResult(TSize, i) = TargetFun(X(i), T(TSize));
end
for  i = 1:TSize
    NumericalResult(i, 1) = TargetFun(X(1),T(i));
    NumericalResult(i, XSize) = TargetFun(X(XSize),T(i));
end

%% 求精确解
if isPlot == 1
    for i = 1:TSize
        for j = 1:XSize
            TargetResult(i, j) = TargetFun(X(j),T(i));
        end
    end
%% 求解误差 && 作图
    ErrorMatrix = abs(TargetResult-NumericalResult);
    AvarageError = sum(sum(ErrorMatrix))/((TSize-1)*(XSize-1));
    disp(AvarageError)
    subplot(3,1,1);
    mesh(X, T, NumericalResult);
    subtitle("Numerical Result");
    xlabel("x");ylabel("t");
    subplot(3,1,2);
    mesh(X, T, TargetResult);
    subtitle("Exact Result");
    xlabel("x");ylabel("t");
    subplot(3,1,3);
    mesh(X, T, ErrorMatrix);
    title("Absolute Error");
    xlabel("x");ylabel("t");zlabel("e");
end
end