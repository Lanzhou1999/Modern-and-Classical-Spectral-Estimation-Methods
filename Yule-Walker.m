clear all;close all 
%--------------------------------------------------------------
% 现代功率谱估计：计算Yule-Walker方程的MSE曲线
% 邢兴润
%--------------------------------------------------------------

%参数设定：N为采样点数,p为方程的阶数,且p一般取[N/3, N/2].
n=1:128; N=length(n);
p=80; 
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);
xn = awgn(xn_,1,'measured');     %信号按照不同信噪比加入噪声
xn = [xn, zeros(1,N)];               %在信号以外补零

%求xn的自相关函数：Rx(0),Rx(2),...,Rx(N-1).
rx = zeros(1,N);
for i=1:N
    rx(i) = xn(1:N)*xn(i:N+i-1)'./N;
end

%求xn的自相关矩阵：pxp方阵.
for i=1:p
     for j=i:p
         R(i,j) = rx(j-i+1);
         R(j,i) = rx(j-i+1);
     end
end
%求Yule-Walker方程右边的向量 RX(1),RX(2),...,Rx(p).
b = -rx(2:p+1)';

%求Yule-Walker方程的系数a0,a1,a2,...,ap.
a(1) = 1;
a(2:p+1) = (R\b)';

%由自相关函数和参数a1,a2,...,ap求系统增益的平方
G2 = rx(1)+rx(2:p+1)*a(2:p+1)';

%计算系统频率响应和系统输出
[H,w] = freqz(1,a',1000);
out = G2.*abs(H).^2;

figure(1);
subplot(2,1,1); plot(w/1000./2,out);
title('Yule-Walker方程估算AR模型功率谱(N = 128; M = 50)');xlabel('f/2pi');
subplot(2,1,2); plot(w./1000./2,out1);
title('处理过的Yule-Walker方程估算AR模型功率谱(N = 128; M = 50)');xlabel('f/2pi');ylabel('dB');











