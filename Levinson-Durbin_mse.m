clear all;close all 
%--------------------------------------------------------------
% ����Levinson-Durbin�㷨��MSE���� 
% ������ ͨ�Ŷ��� 20162410233 
%--------------------------------------------------------------

n=0:128; N=length(n); M=50;
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);

for t0=1:41                                   %t0���������ȡֵ
    for t1=1:1000                             %t1����ÿ��������µļ������
        xn = awgn(xn_,t0-21,'measured');      %�źŰ��ղ�ͬ����ȼ�������
        %��xn������غ���
        xn_1 = [xn, zeros(1,N)];
        rx = zeros(1,N);
        for i=1:N
            rx(i) = xn*xn_1(i:N+i-1)'./N;
        end

        %ѭ����ʼ��,m=0,1��ʱ��ֵ
        a = zeros(1,M+1); a(1) = 1;    % a00 
        a(2) = -rx(2)./rx(1);          % a11
        G2 = rx(1).*(1-a(2).^2);       %һ��ʱ�Ĺ����ܶ�

        for m=2:M
            a(m+1) = -(a(2:m)*fliplr(rx(2:m))' + rx(m+1))./G2; %����am
            a(2:m) = a(2:m) + a(m+1).*fliplr(a(2:m));          %����a1,a2,...,am-1
            G2 = G2.*(1-a(m+1).^2);                            %���¹����ܶ�
        end

        %����ϵͳƵ����Ӧ��ϵͳ���
        [H,w] = freqz(1,a',1000);
        out = G2.*abs(H).^2;

        %������һ����ת��Ϊ�ֱ�
        out=out./abs(max(out));
        out1=(10*log10(abs(out)));

        %�жϷ�ֵ��λ��
        [max1,local1(t1)]=max(out1(1:425));        %Ѱ�����ֵ����0.2*2pi��Ӧ��ֵ�����õ��þֲ����ֵ�±�local1
        [max2,local2(t1)]=max(out1(425:end));      %Ѱ�����ֵ����0.213*2pi��Ӧ��ֵ
        local2(t1)=local2(t1)+424;                 %�õ�f2��ֵ��Ӧ�ֲ����ֵ�Ķ�Ӧ�±�
        
    end

    %��[1,1000]��λ������ת��Ϊ��2pi�Ĺ�һ��Ƶ��
    f1s= local1./1000.*0.5;
    f2s= local2./1000.*0.5;
    f1(t0)=sum(f1s)./1000;
    f2(t0)=sum(f2s)./1000;
    
    %��������������ظ�����Ĵ���ȡƽ��
    mse_f1(t0)= sum((f1s-0.2).^2)./1000
    mse_f2(t0)= sum((f2s-0.213).^2)./1000
end

figure(1);
subplot(2,1,1); plot(-20:20,mse_f1);title('Levinson-Durbin f1=0.2:MSE-SNR;M=50,N=128');xlabel('SNR/dB');ylabel('MSE');
subplot(2,1,2); plot(-20:20,mse_f2);title('Levinson-Durbin f2=0.213:MSE-SNR;M=50,N=128');xlabel('SNR/dB');ylabel('MSE'); 
   
    
