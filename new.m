clear; clc; close all

%% OFDM系统参数
M = 32; % 子载波数
L_cp = 8; % 循环前缀长度
N_sym = 100; % OFDM符号数量
Mod_Order = 16; % 16QAM为4位符号

%% NBI参数
B_I = 4; % NBI带宽
SIR_db = 0; % 信号与干扰比（分贝）
SIR = 10 ^ (SIR_db / 10); % 信号与干扰比（线性）

%% 迭代次数和算法精度门限（根据具体算法实现调整）
max_iter = 10;

%% 生成和调制OFDM数据
data = randi([0, Mod_Order - 1], M, N_sym);
modulated_data = qammod(data, Mod_Order, 'UnitAveragePower', true);
ofdm_signal = ifft(modulated_data, M);

%% 添加循环前缀
ofdm_signal_with_cp = [ofdm_signal(M - L_cp + 1:end, :); ofdm_signal];

%% 生成窄带干扰信号(NBI)
nbi_data = randi([0, 3], round(B_I * N_sym / 4), 1); % QPSK调制，每个符号2bit信息
nbi_modulated = pskmod(nbi_data, 4, pi / 4); % 使用默认的Gray编码
P = round(M + L_cp); % 采样倍增因子P，P=子载波数量+循环前缀长度
nbi_upsampled = upfirdn(nbi_modulated, rcosdesign(0.35, 4, round(P)), P, 1);

%% 调整NBI的采样点以使其与OFDM信号匹配
nbi_signal = zeros(M + L_cp, N_sym);
nbi_signal(:) = nbi_upsampled(1:(M + L_cp) * N_sym);

%% 调整干扰信号功率以匹配SIR
nbi_signal = nbi_signal / norm(nbi_signal) * norm(ofdm_signal_with_cp) / sqrt(SIR);

%% 假设接收信号为OFDM信号与窄带干扰的叠加
rx_signal = ofdm_signal_with_cp + nbi_signal;

%% 接收端处理
y = rx_signal(L_cp + 1:end, :);

%% 初始化估计值
z = zeros(M, N_sym);
d_hat = zeros(M, N_sym);

%% 估计SINR并排序子载波

%% 迭代消除干扰
for k = 1:max_iter
    %% 线性预测与消除干扰
    SINR_est = zeros(M, N_sym);
    for n = 1:N_sym
        z(:, n) = y(:, n) - d_hat(:, n);
        rx_fft = fft(z(:, n), M);
        P_signal = abs(rx_fft).^2;
        P_interference = abs(d_hat(:, n)).^2; % 若已知噪声功率，则需要加上噪声功率
        SINR_est(:, n) = P_signal ./ (P_interference + 1); % SINR估计，加1为噪声功率的简化处理
    end
    % 对每个符号的子载波根据SINR排序
    [~, sorted_indices] = sort(SINR_est, 'descend');

    for n = 1:N_sym
        % 使用前一个子载波的估计值来预测当前子载波的干扰
        for i = 1:M
            idx = sorted_indices(i, n);
            if i == 1
                d_hat(idx, n) = 0; % 对于第一个子载波，没有前一个子载波的信息，可以设置为0或上一个OFDM符号的对应子载波值
            else
                pre_idx = sorted_indices(i-1, n);
                d_hat(idx, n) = d_hat(pre_idx, n); % 使用上一个排好序的子载波的干扰估计值来预测当前子载波的干扰
            end
        end
        % 用估计出的干扰值来消除干扰
        z(:, n) = y(:, n) - d_hat(:, n);
    end
end

%% FFT转换到频域并解调
z_fft = fft(z, M);
demodulated_data = qamdemod(z_fft, Mod_Order, 'UnitAveragePower', true);

%% 评估性能
SER = sum(sum(data ~= demodulated_data)) / (M * N_sym); % 计算误符号率
disp(['Symbol Error Rate is: ', num2str(SER)]);