%% Basic Settings
addpath("tools/");
addpath("work in progress/")
%% Environment variable
sim_num = 3;
M = 4;

% EbN0_dB = 2:2:26;
EsN0_dB = -10:5:25;
EsN0 = db2pow(EsN0_dB);
EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);
iteration = 10^4;  

%% Simulate
BER = zeros(2, length(EbN0_dB));

[BER_tmp, ~] = simulate_golden_code(M, iteration, EbN0_dB);
BER(1,:) = BER_tmp';

[BER_tmp, ~] = simulate_mld(M, iteration, EbN0_dB);
BER(2, :) = BER_tmp';

eta = 2;
n = 2;
BER_tmp = simulate_modulation_diversity(eta, n, iteration, EbN0_dB)
BER(3,:) = BER_tmp';

%% Plot
BER_Title = sprintf("BER for %d-QAM", M);
% SER_Title = sprintf("SER for %d-QAM", M);
x_axis = "Es/N0 (dB)";

legend_order = ["Golden Code (MLD)", "Uncoded (MLD)", "Mod Diversity"];

myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1]);
xlim([2 26]);
% myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
% ylim([10^(-6) 1])