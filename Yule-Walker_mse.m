clear all;close all 
%--------------------------------------------------------------
% 作业二,现代功率谱估计：计算Yule-Walker方程的MSE曲线
% 邢兴润 通信二班 20162410233
%--------------------------------------------------------------

%参数设定：N为采样点数,p为方程的阶数,且p一般取[N/3, N/2].
n=1:128; N=length(n);
p=80; 
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);

for t0=1:41                                  %t0决定信噪比取值
    for t1=1:1000                            %t1决定每个信噪比下的计算次数
        xn = awgn(xn_,t0-31,'measured');     %信号按照不同信噪比加入噪声
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

        %判断峰值点位置
        [max1,local1(t1)]=max(out(1:425));        %寻找最大值估计0.2*2pi对应峰值，并得到该局部最大值下标local1
        [max2,local2(t1)]=max(out(425:end));      %寻找最大值估计0.213*2pi对应峰值
        local2(t1)=local2(t1)+424;                %得到f2峰值对应局部最大值的对应下标

    end
    
    %把[1,1000]的位置坐标转化为对2pi的归一化频率
    f1s= local1./1000.*0.5;
    f2s= local2./1000.*0.5;
    f1(t0)=sum(f1s)./1000;
    f2(t0)=sum(f2s)./1000;
    
    %计算均方误差，并对重复计算的次数取平均
    mse_f1(t0)= sum((f1s-0.2).^2)./1000
    mse_f2(t0)= sum((f2s-0.213).^2)./1000
end

figure(1);
subplot(2,1,1); plot(-30:10,mse_f1);title('Yule-Walker f1=0.2:MSE-SNR(N=128,p=80)');xlabel('SNR/dB');ylabel('MSE');
subplot(2,1,2); plot(-30:10,mse_f2);title('Yule-Walker f2=0.213:MSE-SNR(N=128,p=80)');xlabel('SNR/dB');ylabel('MSE'); 
