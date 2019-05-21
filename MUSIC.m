clear all;close all 
%--------------------------------------------------------------
% 使用特征值分解的MUSIC算法估算AR模型功率谱 
% 邢兴润
%--------------------------------------------------------------

%参数设定：N为采样点数,p为自相关矩阵阶数.
n=0:100; N=length(n);

xn = sqrt(20)*sin(2*pi*0.2*n) + sqrt(2)*sin(2*pi*0.213*n);
xn = awgn(xn,1,'measured');
xn = [xn, zeros(1,N)];

%求xn的自相关函数：Rx(0),Rx(2),...,Rx(N-1).
rx = zeros(1,N);
for i=1:N
    rx(i) = xn(1:N)*xn(i:N+i-1)'./N;
end

%求xn的自相关矩阵：pxp方阵.
for i=1:N-1
    for j=i:N-1
        R(i,j) = rx(j-i+1);
        R(j,i) = rx(j-i+1);
    end
end

%求自相关矩阵的特征值
[V,D_]=eig(R);                                    %V是特征向量 D_为对角矩阵形式的特征值
D=eig(R);                                         %用向量的形式表示特征值

%根据特征值从大到小的顺序对特征向量排序
[D_sort,D_index]= sort(D,'descend');              %排序，D_sort是排序后的特征值，D_index是排序的原序号
V_sort=V(:, D_index);                             %V_sort对应排序后的特征向量

%去掉特征值较大的两个特征向量
V = V_sort(:,3:end);

ew = zeros(1,N-1);
for w=1:500                                       %w的采样点数为500
    %求序列e(w)
    for i=1:N-1
        ew(i) = exp(1i*(i-1).*w.*2.*pi./1000);
    end
    
    %求函数Px(w)的分母
    for i=1:N-3
        t = abs(ew*V(:,i)).^2;
    end
    
    Px(w) = 1./t;                                  %计算Px(w)
end

figure(1);
subplot(1,1,1); plot((1:500)/1000,Px);title('MUSIC算法估算AR模型功率谱(N = 128; p=80)');xlabel('f/2pi');




