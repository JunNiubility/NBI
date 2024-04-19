% OFDM消除NBI仿真代码
clear; clc; close all

%% OFDM系统参数
M = 32; % 子载波数
L_cp = 8; % 循环前缀长度
N_sym = 100; % OFDM符号数量
Mod_Order = 16; % 16QAM为4位符号

%% NBI参数
B_I = 4; % NBI带宽
% SIR_db = 3; % 信号与干扰比（分贝）
% SIR = 10 ^ (SIR_db / 10); % 信号与干扰比（线性）

%% 生成和调制OFDM数据
data = randi([0, Mod_Order - 1], M, N_sym);
modulated_data = qammod(data, Mod_Order, 'UnitAveragePower', true);
ofdm_signal = ifft(modulated_data, M);

%% 添加循环前缀
ofdm_signal_with_cp = [ofdm_signal(M - L_cp + 1:end, :); ofdm_signal];

%% 生成窄带干扰信号(NBI)
nbi_data = randi([0, 3], round(B_I * N_sym / 4), 1); % QPSK调制，每个符号2bit信息
nbi_modulated = pskmod(nbi_data, 4, pi / 4); % 使用默认的Gray编码
P = round(M + L_cp); % 采样倍增因子P，这里是一个示例值，需要根据NBI的初始采样点数来调整
nbi_upsampled = upfirdn(nbi_modulated, rcosdesign(0.35, 4, round(P)), P, 1);

%% 调整NBI的采样点以使其与OFDM信号匹配
nbi_signal = zeros(M + L_cp, N_sym);
nbi_signal(:) = nbi_upsampled(1:(M + L_cp) * N_sym);

%% 调整干扰信号功率以匹配SIR
for SIR_db = 15:-3:15 % 信号与干扰比（分贝）
    fprintf("--------当前SIR=%ddB--------\n", SIR_db);
    SIR = 10 ^ (SIR_db / 10); % 信号与干扰比（线性）
    nbi_signal = nbi_signal / norm(nbi_signal) * norm(ofdm_signal_with_cp) / sqrt(SIR);

    %% 假设接收信号为OFDM信号与窄带干扰的叠加
    rx_signal = ofdm_signal_with_cp + nbi_signal;

    %% 接收端处理
    y = rx_signal(L_cp + 1:end, :);

    %% 动态计算预测器系数alpha
    % 初始化接收信号的估计值和干扰的估计值
    z = zeros(M, N_sym);
    d_hat = zeros(M, N_sym);
    delta_d = zeros(1, M);
    SINR_est = zeros(M, N_sym);

    %% FFT转换到频域
    y = fft(y, M);
    % 算法迭代
    for j = 1:N_sym
        % rx_fft = fft(z(:, j), M);
        % P_signal = abs(rx_fft) .^ 2;
        % P_interference = abs(d_hat(:, j)) .^ 2; % 若已知噪声功率，则需要加上噪声功率
        % SINR_est(:, j) = P_signal ./ (P_interference + 1); % SINR估计，加1为噪声功率的简化处理
        % %% 对每个符号的子载波根据SINR排序
        % [~, sorted_indices] = sort(SINR_est, 'descend');

        for i = 1:M
            % idx = sorted_indices(i, j);
            % 移除干扰估计
            z(i, j) = y(i, j) - d_hat(i, j);
            s_hat = qammod(qamdemod(z(i, j), Mod_Order, 'UnitAveragePower', true), Mod_Order, 'UnitAveragePower', true); % 解调为原始调制符号
            % 计算干扰残差
            delta_d(i) = y(i, j) - s_hat;
            % 更新干扰估计
            if i < M
                % 计算R_d_k和q_k参数
                R = xcorr(delta_d(1:i), conj(delta_d(1:i)'), i - 1);
                q = xcorr(delta_d(1:i), y(i + 1, j), i - 1);
                R_d_k = toeplitz(R(i:end));
                q_k = q(i:end)';
                % 根据alpha_k = R_d_k^-1 * q_k公式来更新滤波器系数alpha
                if det(R_d_k) ~= 0
                    alpha = R_d_k \ q_k;
                else
                    alpha = 0;
                end

                % 线性预测干扰
                d_hat(i + 1, j) = alpha' * delta_d(1:i)';
            else

            end

        end

    end

    %% 解调
    %使用对消NBI算法
    demodulated_data_remove = qamdemod(z, Mod_Order, 'UnitAveragePower', true);
    %不消除干扰
    demodulated_data = qamdemod(y, Mod_Order, 'UnitAveragePower', true);

    %% 评估性能
    %使用对消NBI算法
    SER = sum(sum(data ~= demodulated_data_remove)) / (M * N_sym); % 计算误符号率
    disp(['Symbol Error Rate Of remove NBI is: ', num2str(SER)]);
    %不消除干扰
    SER = sum(sum(data ~= demodulated_data)) / (M * N_sym); % 计算误符号率
    disp(['Symbol Error Rate is: ', num2str(SER)]);
end
