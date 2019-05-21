clear all;close all 
%--------------------------------------------------------------
% 使用Levinson-Durbin估计AR模型功率谱 
% 邢兴润 
%--------------------------------------------------------------
n=0:128; N=length(n); M=50;
xn = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);
xn = awgn(xn,10,'measured');

%求xn的自相关函数
xn_ = [xn, zeros(1,N)];
rx = zeros(1,N);
for i=1:N
    rx(i) = xn*xn_(i:N+i-1)'./N;
end

%循环初始化,m=0,1阶时的值
a = zeros(1,M+1); a(1) = 1; % a00 
a(2) = -rx(2)./rx(1); % a11
G2 = rx(1).*(1-a(2).^2);

for m=2:M
    a(m+1) = -(a(2:m)*fliplr(rx(2:m))' + rx(m+1))./G2;
    a(2:m) = a(2:m) + a(m+1).*fliplr(a(2:m));
    G2 = G2.*(1-a(m+1).^2);
end

%计算系统频率响应和系统输出
[H,w] = freqz(1,a',1000);
out = G2.*abs(H).^2;

%归一化并转化为分贝
out=out./abs(max(out));
out1=(10*log10(abs(out)));

figure(1);
subplot(2,1,1); plot(w/1000./2,out);
title('Levinson-Durbin算法估算AR模型功率谱(N = 128; M = 50)');xlabel('f/2pi');
subplot(2,1,2); plot(w./1000./2,out1);
title('处理过的Levinson-Durbin算法估算AR模型功率谱(N = 128; M = 50)');xlabel('f/2pi');ylabel('dB');
