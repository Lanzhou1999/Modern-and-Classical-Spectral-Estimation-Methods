clear all;close all 
%--------------------------------------------------------------
% 使用burg算法估算AR模型功率谱 
% 邢兴润 通信二班 20162410233 
%--------------------------------------------------------------
n=0:128; N=length(n);
xn = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);
xn = awgn(xn,10,'measured');

%0阶时循环初始化条件
ef = xn;
eb = xn;
a = 1;
G2 = xn*xn'./N;
for m=1:50
    efm = ef(2:end);               %m-1阶前向预测误差的有用部分,n属于[m,N-1]
    ebm = eb(1:end - 1);           %m-1阶后向预测误差的有用部分,n-1属于[m-1,N-2]

    km = (-2.*sum(ebm.*efm))./sum(efm.*efm + ebm.*ebm);      % 第m阶反射系数
    
    %更新前后项误差
    ef = efm + km.*ebm;
    eb = ebm + km.*efm;
    %更新系数矩阵和预估误差功率
    a = [a; 0] + km*[0; flipud(a)];
    G2 = (1 - km*km)*G2;
end

%计算系统频率响应和系统输出
[H,w] = freqz(1,a');               %默认值为512点，如果需要可以设定为1024等，增加图像分辨率
out = G2.*abs(H).^2;               %系统函数和系统增益相乘，系统增益包含了噪声功率

%后处理，归一化并转化为分贝
out=out./abs(max(out));
out1=(10*log10(abs(out)));

figure(1);
subplot(2,1,1); plot(w./512.*0.5,out);
title('burg算法估算AR模型功率谱(N = 128; M = 50)');xlabel('f/2pi');
subplot(2,1,2); plot(w./512.*0.5,out1);
title('处理过的burg算法估算AR模型功率谱(N = 128; M = 50)');xlabel('f/2pi');ylabel('dB');
