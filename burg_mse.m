clear all;close all 
%--------------------------------------------------------------
% 计算burg算法的MSE曲线 
% 邢兴润 通信二班 20162410233 
%--------------------------------------------------------------
n=0:128; N=length(n);
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);

for t0=1:61                                  %t0决定信噪比取值
    for t1=1:1000                             %t1决定每个信噪比下的计算次数
        xn = awgn(xn_,t0-31,'measured');     %信号按照不同信噪比加入噪声
        
        %0阶时循环初始化条件
        ef = xn;
        eb = xn;
        a = 1;
        G2 = xn*xn'./N;
        for m=1:100                          %m为burg算法取的阶次
            efm = ef(2:end);                 %m-1阶前向预测误差的有用部分,n属于[m,N-1]
            ebm = eb(1:end - 1);             %m-1阶后向预测误差的有用部分,n-1属于[m-1,N-2]

            km = (-2.*sum(ebm.*efm))./sum(efm.*efm + ebm.*ebm); % 第m阶反射系数
    
            %更新前后项误差
            ef = efm + km.*ebm;
            eb = ebm + km.*efm;
            %更新系数矩阵和预估误差功率
            a = [a; 0] + km*[0; flipud(a)];
            G2 = (1 - km*km)*G2;
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
subplot(2,1,1); plot(-30:30,mse_f1);title('burg f1=0.2:MSE-SNR');xlabel('SNR/dB');ylabel('MSE');
subplot(2,1,2); plot(-30:30,mse_f2);title('burg f2=0.213:MSE-SNR');xlabel('SNR/dB');ylabel('MSE'); 
