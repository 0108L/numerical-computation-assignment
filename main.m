clear
clc
% MIMO System
% 发射天线数NT,接收天线数NR,发射矩阵长度L
% x = H*c+v
NT=4;NR=4;L=1000;
SNR = 0:1:20; % 信噪比,单位为(dB)
c_real=randi([0,1],NT,L); % 生成内部元素为0和1的NT*L发射信号
c=zeros(NT, L); % 经过V-BLAST计算后的发射信号

% 实际发射信号的0转化为-1, 1仍保持为1
X = (-1).^(c_real + 1);

% 快速衰弱的NR*NT*L维瑞利通道
H=sqrt(1/2)*(randn(NR,NT,L)+1i*randn(NR,NT,L));
% 服从正态分布的高斯白噪声v,均值为0,方差为1,NR*1维
v=sqrt(1/2)*(randn(NR,L)+1i*randn(NR,L));

% 未叠加噪声的噪声信号
x=zeros(NR,L);
for i=1:L
    x(:,i)=sqrt(1/2)*H(:,:,i)*X(:,i);
end

% 先绘制第一个图床
figure;
disp('文献1中的算法: ');

%%%%%%%%% V-blast算法 %%%%%%%%%%
disp('V-blast算法');
% 不同信噪比下的误码率
vector=[];
% 叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    % 解码得到的信号
    c=V_blast(H,x_noised);                     % 在V_blast.m函数文件中
    % 计算USQR算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    vector=[vector,err_ratio];
end
semilogy(SNR,vector,'x-r'); % 红色叉
hold on;

%%%%%%%%%% 无排序QR算法 %%%%%%%%%%%
disp('未排序的QR算法');
%不同信噪比下的误码率
vector_usqr=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=USQR(H,x_noised);                         %在USQR.m函数文件中
    %计算USQR算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    vector_usqr=[vector_usqr,err_ratio];
end
semilogy(SNR,vector_usqr,'--k'); %黑色星号
hold on;    

%%%%%%%%%%% 排序QR算法 %%%%%%%%%%
disp('排序QR算法');
%不同信噪比下的误码率
vector_sqrd=[];
%叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    %经解码得到的信号
    c=SQRD(H,x_noised);                     % 在SQRD.m函数文件中
    %计算SQRD算法的误码率
    [errbit,err_ratio]=biterr(c_real,c);
    vector_sqrd=[vector_sqrd,err_ratio];
end
semilogy(SNR,vector_sqrd,'o-g'); %绿色圆圈
hold on; 

xlabel('信噪比SNR');
ylabel('比特误码率BER');
title('NT=4，NR=4时,V-blast算法和有无排序的QR算法的比特误码率和信噪比关系曲线');
legend('V-blast','未排序 QR','排序 QR');

% 绘制第二个图床
figure;
disp('文献2中的算法: ');

%%%%%%%%%% MMSE算法 %%%%%%%%%%
disp('MMSE算法');
% 不同信噪比下的误码率
vector=[];
% 叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    % 解码得到的信号
    c=MMSE(H,x_noised,snr);                     % 在MMSE.m函数文件夹中
    % 计算V-blast算法的误码率
    frame_error = 0;
    for j = 1:L
        if any(biterr(c_real(:,j), c(:,j)))
            frame_error = frame_error + 1;
        end
    end
    err_ratio = frame_error / L;
    [errbit,err_ratio]=biterr(c_real,c);
    vector=[vector,err_ratio];
end
semilogy(SNR,vector,'d-r'); % 红色菱形
hold on;

%%%%%%%%%% ZF算法 %%%%%%%%%%%
disp('ZF算法');
% 不同信噪比下的误码率
vector=[];
% 叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    % 解码后得到信号
    c=ZF(H,x_noised);                       % 在ZF.m函数文件夹中
    % 计算ZF算法的误码率
    frame_error = 0;
    for j = 1:L
        if any(biterr(c_real(:,j), c(:,j)))
            frame_error = frame_error + 1;
        end
    end
    err_ratio = frame_error / L;
    [errbit,err_ratio]=biterr(c_real,c);
    vector=[vector,err_ratio];
end
semilogy(SNR,vector,'o-g'); % 绿色虚线
hold on;

%%%%%%%%%% MMSE_QR算法 %%%%%%%%%%
disp('MMSE-QR算法');
% 不同信噪比下的误码率
vector=[];
% 叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    % 经解码得到的信号
    c=MMSE_QR(H,x_noised,snr);                      %在MMSE_QR.m
    % 计算MMSE_QR算法的误码率
    frame_error = 0;
    for j = 1:L
        if any(biterr(c_real(:,j), c(:,j)))
            frame_error = frame_error + 1;
        end
    end
    err_ratio = frame_error / L;
    [errbit,err_ratio]=biterr(c_real,c);
    vector=[vector,err_ratio];
end
semilogy(SNR,vector,'--m'); % 紫红色实体
hold on;

%%%%%%%%%%% MMSE-SQRD算法 %%%%%%%%%%
disp('排序MMSE-QR算法');
% 不同信噪比下的误码率
vector=[];
% 叠加噪声
for m=SNR
    snr=10^(m/10);
    x_noised=x+sqrt(1/snr)*v;
    % 经解码得到的信号
    c=MMSE_SQRD(H,x_noised,snr);                        % 在MMSE_SQRD.m函数文件中
    % 计算MMSE_QR算法的帧误码率
    frame_error = 0;
    for j = 1:L
        if any(biterr(c_real(:,j), c(:,j)))
            frame_error = frame_error + 1;
        end
    end
    err_ratio = frame_error / L;
    [errbit,err_ratio]=biterr(c_real,c);
    vector=[vector,err_ratio];
end
semilogy(SNR,vector,'*-k'); % 黑色实线
hold on;

xlabel('信噪比SNR');
ylabel('帧误码率FER');
title('NT=4，NR=4时,MMSE各种算法与ZF算法的帧误码率和信噪比关系曲线');
legend('MMSE', 'ZF', 'MMSE-QR','排序MMSE-QR');
