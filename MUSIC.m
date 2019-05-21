clear all;close all 
%--------------------------------------------------------------
% ʹ������ֵ�ֽ��MUSIC�㷨����ARģ�͹����� 
% ������
%--------------------------------------------------------------

%�����趨��NΪ��������,pΪ����ؾ������.
n=0:100; N=length(n);

xn = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);
xn = awgn(xn,1,'measured');
xn = [xn, zeros(1,N)];

%��xn������غ�����Rx(0),Rx(2),...,Rx(N-1).
rx = zeros(1,N);
for i=1:N
    rx(i) = xn(1:N)*xn(i:N+i-1)'./N;
end

%��xn������ؾ���pxp����.
for i=1:N-1
    for j=i:N-1
        R(i,j) = rx(j-i+1);
        R(j,i) = rx(j-i+1);
    end
end

%������ؾ��������ֵ
[V,D_]=eig(R);                                    %V���������� D_Ϊ�ԽǾ�����ʽ������ֵ
D=eig(R);                                         %����������ʽ��ʾ����ֵ

%��������ֵ�Ӵ�С��˳���������������
[D_sort,D_index]= sort(D,'descend');              %����D_sort������������ֵ��D_index�������ԭ���
V_sort=V(:, D_index);                             %V_sort��Ӧ��������������

%ȥ������ֵ�ϴ��������������
V = V_sort(:,3:end);

ew = zeros(1,N-1);
for w=1:500                                       %w�Ĳ�������Ϊ500
    %������e(w)
    for i=1:N-1
        ew(i) = exp(1i*(i-1).*w.*2.*pi./1000);
    end
    
    %����Px(w)�ķ�ĸ
    for i=1:N-3
        t = abs(ew*V(:,i)).^2;
    end
    
    Px(w) = 1./t;                                  %����Px(w)
end

figure(1);
subplot(1,1,1); plot((1:500)/1000,Px);title('MUSIC�㷨����ARģ�͹�����(N = 128; p=80)');xlabel('f/2pi');




