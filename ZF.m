function c = ZF(H, x)
    % 使用零强制（ZF）准则对空时信号进行解码
    % H -- NR*NT维的瑞利信道矩阵
    % x -- 接收到的信号向量
    % c -- 解码后的信号向量

    [NR, NT, L] = size(H);
    c = zeros(NT, L);

    for j = 1:L
        HH = H(:,:,j);
        % G -- 滤波矩阵
        G = inv(HH' * HH) * HH';
        y = G * x(:, j);
        % 根据阈值决策生成解码信号
        c(:, j) = (y >= 0) - (y < 0);
    end

    % 将解码信号转换为二进制形式
    c = (c + 1) / 2;
end
