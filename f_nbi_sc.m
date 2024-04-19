clear; clc; close all

% OFDM系统参数
M = 32; % 子载波数
L_cp = 8; % 循环前缀长度
N_sym = 100; % OFDM符号数量
Mod_Order = 16; % 16QAM为4位符号

% NBI参数
B_I = 4; % NBI带宽
SIR_db = 0; % 信号与干扰比（分贝）
SIR = 10 ^ (SIR_db / 10); % 信号与干扰比（线性）

% 迭代次数和算法精度门限（根据具体算法实现调整）
max_iter = 10;
mu = 0.01; % LMS算法步长

% 生成和调制OFDM数据
data = randi([0, Mod_Order - 1], M, N_sym);
modulated_data = qammod(data, Mod_Order, 'UnitAveragePower', true);
ofdm_signal = ifft(modulated_data, M);

% 添加循环前缀
ofdm_signal_with_cp = [ofdm_signal(M - L_cp + 1:end, :); ofdm_signal];

% 生成窄带干扰信号(NBI)
nbi_data = randi([0, 3], round(B_I * N_sym / 4), 1); % QPSK调制，每个符号2bit信息
nbi_modulated = pskmod(nbi_data, 4, pi / 4); % 使用默认的Gray编码
P = round(M + L_cp); % 采样倍增因子P，这里是一个示例值，需要根据NBI的初始采样点数来调整
nbi_upsampled = upfirdn(nbi_modulated, rcosdesign(0.35, 4, round(P)), P, 1);

% 调整NBI的采样点以使其与OFDM信号匹配
nbi_signal = zeros(M + L_cp, N_sym);
nbi_signal(:) = nbi_upsampled(1:(M + L_cp) * N_sym);

% 调整干扰信号功率以匹配SIR
nbi_signal = nbi_signal / norm(nbi_signal) * norm(ofdm_signal_with_cp) / sqrt(SIR);

% 假设接收信号为OFDM信号与窄带干扰的叠加
rx_signal = ofdm_signal_with_cp + nbi_signal;

%% 接收端处理
% 移除循环前缀
y = rx_signal(L_cp + 1:end, :);

% F-NBI-SC算法初始化
z = zeros(1, M); % 初始化接收信号的估计值
d_hat = zeros(1, M); % 初始化干扰的估计值
alpha = zeros(1, M); % 初始化每个子载波的alpha值

% 算法迭代
for sym_idx = 1:N_sym
    for k = 1:max_iter
        for i = 1:M
            z(i) = y(i, sym_idx) - d_hat(i); % 移除干扰估计
            s_hat = qamdemod(z(i), Mod_Order, 'UnitAveragePower', true); % 解调为原始调制符号
            e = s_hat - z(i); % 误差计算
            alpha(i) = alpha(i) + mu * e * conj(y(i, sym_idx)); % 自适应更新alpha值
            d_hat(i) = alpha(i) * y(i, sym_idx); % 更新干扰估计值
        end
    end
end

% z= zeros(M,N_sym);
% for i=1:M
%     z(i,:)=y(i,:) - d_hat(i);
% end

% FFT转换到频域
rx_fft = fft(y, M);

% 对接收的OFDM信号进行16-QAM解调
demodulated_data = qamdemod(rx_fft, Mod_Order, 'UnitAveragePower', true);

% 评估性能
SER = sum(sum(data ~= demodulated_data)) / (M * N_sym); % 计算误符号率
disp(['Symbol Error Rate is: ', num2str(SER)]);

%% 1. 显示原始OFDM信号的频谱
figure;
subplot(4, 1, 1);
ofdm_signal_freq = fftshift(fft(ofdm_signal, M));
plot(1:M, abs(ofdm_signal_freq(:,1)));
title('Original OFDM Signal Spectrum');
xlabel('Subcarrier Index');
ylabel('Amplitude');

%% 2. 显示带有NBI的接收信号频谱
subplot(4, 1, 2);
rx_signal_freq = fftshift(fft(rx_signal, M));
plot(1:M, abs(rx_signal_freq(:,1))); % 只展示一个OFDM符号的频谱，以保持图示简洁
title('Received Signal Spectrum with NBI');
xlabel('Subcarrier Index');
ylabel('Amplitude');

%% 3. 展示干扰消除后的接收信号频谱
subplot(4, 1, 3);
y_no_cp_fft = fftshift(fft(y, M));
plot(1:M, abs(y_no_cp_fft(:,1))); % 只展示一个OFDM符号的频谱，以保持图示简洁
title('Received Signal Spectrum After NBI Cancellation');
xlabel('Subcarrier Index');
ylabel('Amplitude');

%% 4. 比较原始数据和解调后数据的符号错误率（SER）
subplot(4, 1, 4);
scatter(real(modulated_data(:)), imag(modulated_data(:)), 'b');
hold on;
scatter(real(rx_fft(:)), imag(rx_fft(:)), 'r*');
legend('Original 16QAM Constellation', 'Demodulated Constellation after NBI Cancellation');
title('Constellation Diagram Comparison');
xlabel('In-phase Amplitude');
ylabel('Quadrature Amplitude');