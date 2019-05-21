clear all;close all 
%--------------------------------------------------------------
% ����burg�㷨��MSE���� 
% ������ ͨ�Ŷ��� 20162410233 
%--------------------------------------------------------------
n=0:128; N=length(n);
xn_ = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);

for t0=1:61                                  %t0���������ȡֵ
    for t1=1:1000                             %t1����ÿ��������µļ������
        xn = awgn(xn_,t0-31,'measured');     %�źŰ��ղ�ͬ����ȼ�������
        
        %0��ʱѭ����ʼ������
        ef = xn;
        eb = xn;
        a = 1;
        G2 = xn*xn'./N;
        for m=1:100                          %mΪburg�㷨ȡ�Ľ״�
            efm = ef(2:end);                 %m-1��ǰ��Ԥ���������ò���,n����[m,N-1]
            ebm = eb(1:end - 1);             %m-1�׺���Ԥ���������ò���,n-1����[m-1,N-2]

            km = (-2.*sum(ebm.*efm))./sum(efm.*efm + ebm.*ebm); % ��m�׷���ϵ��
    
            %����ǰ�������
            ef = efm + km.*ebm;
            eb = ebm + km.*efm;
            %����ϵ�������Ԥ������
            a = [a; 0] + km*[0; flipud(a)];
            G2 = (1 - km*km)*G2;
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
subplot(2,1,1); plot(-30:30,mse_f1);title('burg f1=0.2:MSE-SNR');xlabel('SNR/dB');ylabel('MSE');
subplot(2,1,2); plot(-30:30,mse_f2);title('burg f2=0.213:MSE-SNR');xlabel('SNR/dB');ylabel('MSE'); 
