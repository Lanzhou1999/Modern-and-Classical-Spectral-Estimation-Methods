clear all;close all 
%--------------------------------------------------------------
% 计算Levinson-Durbin算法的MSE曲线 
% 邢兴润 通信二班 20162410233 
%--------------------------------------------------------------

n=0:128; N=length(n); M=50;
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);

for t0=1:41                                   %t0决定信噪比取值
    for t1=1:1000                             %t1决定每个信噪比下的计算次数
        xn = awgn(xn_,t0-21,'measured');      %信号按照不同信噪比加入噪声
        %求xn的自相关函数
        xn_1 = [xn, zeros(1,N)];
        rx = zeros(1,N);
        for i=1:N
            rx(i) = xn*xn_1(i:N+i-1)'./N;
        end

        %循环初始化,m=0,1阶时的值
        a = zeros(1,M+1); a(1) = 1;    % a00 
        a(2) = -rx(2)./rx(1);          % a11
        G2 = rx(1).*(1-a(2).^2);       %一阶时的功率密度

        for m=2:M
            a(m+1) = -(a(2:m)*fliplr(rx(2:m))' + rx(m+1))./G2; %更新am
            a(2:m) = a(2:m) + a(m+1).*fliplr(a(2:m));          %更新a1,a2,...,am-1
            G2 = G2.*(1-a(m+1).^2);                            %更新功率密度
        end

        %计算系统频率响应和系统输出
        [H,w] = freqz(1,a',1000);
        out = G2.*abs(H).^2;

        %后处理，归一化并转化为分贝
        out=out./abs(max(out));
        out1=(10*log10(abs(out)));

        %判断峰值点位置
        [max1,local1(t1)]=max(out1(1:425));        %寻找最大值估计0.2*2pi对应峰值，并得到该局部最大值下标local1
        [max2,local2(t1)]=max(out1(425:end));      %寻找最大值估计0.213*2pi对应峰值
        local2(t1)=local2(t1)+424;                 %得到f2峰值对应局部最大值的对应下标
        
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
subplot(2,1,1); plot(-20:20,mse_f1);title('Levinson-Durbin f1=0.2:MSE-SNR;M=50,N=128');xlabel('SNR/dB');ylabel('MSE');
subplot(2,1,2); plot(-20:20,mse_f2);title('Levinson-Durbin f2=0.213:MSE-SNR;M=50,N=128');xlabel('SNR/dB');ylabel('MSE'); 
   
    
