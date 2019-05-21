clear all;close all 
%--------------------------------------------------------------
% ��ҵ��,�ִ������׹��ƣ�����Yule-Walker���̵�MSE����
% ������ ͨ�Ŷ��� 20162410233
%--------------------------------------------------------------

%�����趨��NΪ��������,pΪ���̵Ľ���,��pһ��ȡ[N/3, N/2].
n=1:128; N=length(n);
p=80; 
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);

for t0=1:41                                  %t0���������ȡֵ
    for t1=1:1000                            %t1����ÿ��������µļ������
        xn = awgn(xn_,t0-31,'measured');     %�źŰ��ղ�ͬ����ȼ�������
        xn = [xn, zeros(1,N)];               %���ź����ⲹ��

        %��xn������غ�����Rx(0),Rx(2),...,Rx(N-1).
        rx = zeros(1,N);
        for i=1:N
            rx(i) = xn(1:N)*xn(i:N+i-1)'./N;
        end

        %��xn������ؾ���pxp����.
        for i=1:p
            for j=i:p
                R(i,j) = rx(j-i+1);
                R(j,i) = rx(j-i+1);
            end
        end
        %��Yule-Walker�����ұߵ����� RX(1),RX(2),...,Rx(p).
        b = -rx(2:p+1)';

        %��Yule-Walker���̵�ϵ��a0,a1,a2,...,ap.
        a(1) = 1;
        a(2:p+1) = (R\b)';

        %������غ����Ͳ���a1,a2,...,ap��ϵͳ�����ƽ��
        G2 = rx(1)+rx(2:p+1)*a(2:p+1)';

        %����ϵͳƵ����Ӧ��ϵͳ���
        [H,w] = freqz(1,a',1000);
        out = G2.*abs(H).^2;

        %�жϷ�ֵ��λ��
        [max1,local1(t1)]=max(out(1:425));        %Ѱ�����ֵ����0.2*2pi��Ӧ��ֵ�����õ��þֲ����ֵ�±�local1
        [max2,local2(t1)]=max(out(425:end));      %Ѱ�����ֵ����0.213*2pi��Ӧ��ֵ
        local2(t1)=local2(t1)+424;                %�õ�f2��ֵ��Ӧ�ֲ����ֵ�Ķ�Ӧ�±�

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
subplot(2,1,1); plot(-30:10,mse_f1);title('Yule-Walker f1=0.2:MSE-SNR(N=128,p=80)');xlabel('SNR/dB');ylabel('MSE');
subplot(2,1,2); plot(-30:10,mse_f2);title('Yule-Walker f2=0.213:MSE-SNR(N=128,p=80)');xlabel('SNR/dB');ylabel('MSE'); 
