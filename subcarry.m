% 假设
% OFDM子载波间隔为 DeltaF
% 中心频率为 Fc（相对于OFDM的基带位置）
% NBI带宽为 BW_NBI

N = 256; % 例子中的总子载波数
DeltaF = 15e3; % 子载波间隔
Fc = 3.75e5; % NBI信号的中心频率
BW_NBI = 30e3; % NBI带宽

% 计算子载波的频率点
subcarrier_freqs = ((-(N/2):((N/2)-1)) * DeltaF).';

% 确定NBI的频率范围
nbi_freq_range = [(Fc - BW_NBI/2) (Fc + BW_NBI/2)];

% 查找受影响的子载波
affected_subcarriers = subcarrier_freqs(subcarrier_freqs >= nbi_freq_range(1) & subcarrier_freqs <= nbi_freq_range(2));