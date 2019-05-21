clear all;close all 
%--------------------------------------------------------------
%作业一,古典功率谱估计：画出在不同信噪比下，估计频率的均方误差曲线 
%邢兴润 通信二班 20162410233 
%--------------------------------------------------------------

n=0:255; 
xn=exp(i*(2*pi*0.1*n+pi/3))+10*exp(i*(2*pi*0.01*n+pi/4));

%定义全零矩阵占位，以提高计算速度 
mse_f1=zeros(1,50); 
mse_f2=zeros(1,50); 
local1=zeros(1,1000); 
local2=zeros(1,1000); 
f1_=zeros(1,20); 
f2_=zeros(1,20);

for j=1:50
    for i=1:1000        
        xn_=awgn(xn,j-31,'measured');  %给信号xn加上不同信噪比的噪声：[-30,20]dB
        xk=abs(fft(xn_,256));          %对信号做256点傅里叶变换
                
        [max1,local1(i)]=max(xk(1:20));        %在[1/128,20/128]pi频率内寻找最大值可以找到0.01*2pi对应峰值，并得到该局部最大值下标local1
        [max2,local2(i)]=max(xk(21:100));      %在[21/128,100/128]pi频率内寻找最大值可以找到0.1*2pi对应峰值
        local2(i)=local2(i)+20;                %得到峰值对应局部最大值对应的下标
    end
    
    f1=local1/256;                             %f1数字域频率转化为频率
    f2=local2/256;                             %f2数字域频率转化为频率
    f1_(j)=sum(f1)/1000;                       %f1的均值
    f2_(j)=sum(f2)/1000;                       %f2的均值
    mse_f1(j)=sum((f1-0.01).^2)/1000;          %f1的均方误差
    mse_f2(j)=sum((f2-0.1).^2)/1000;           %f2的均方误差
end

figure(1);
subplot(2,1,1); plot(-30:19,mse_f1); title('f1=0.01:MSE-SNR');xlabel('SNR/dB');ylabel('MSE'); 
subplot(2,1,2); plot(-30:19,mse_f2); title('f2=0.1:MSE-SNR');xlabel('SNR/dB');ylabel('MSE'); 
