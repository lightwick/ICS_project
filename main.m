close all
clc
clear all

%% Basic Settings
addpath("tools/");
addpath("work in progress/")
%% Environment variable
sim_num = 4;
M = 4;

% EbN0_dB = 2:2:26;
EsN0_dB = 0:4:24;
EsN0 = db2pow(EsN0_dB);
EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);
iteration = 10^4;

%% Simulate
BER = zeros(sim_num, length(EbN0_dB));

[BER_tmp, ~] = simulate_golden_code(M, iteration, EbN0_dB);
BER(1,:) = BER_tmp';

[BER_tmp, ~] = simulate_mld(M, iteration, EbN0_dB);
BER(2, :) = BER_tmp';

eta = 2;
n = 2;
BER_tmp = simulate_modulation_diversity(eta, n, iteration, EbN0_dB);
BER(3,:) = BER_tmp';

eta = 2;
n = 2;
BER_tmp = simulate_modulation_diversity_reduced(eta, n, iteration, EbN0_dB);
BER(4,:) = BER_tmp';

%% Plot
BER_Title = sprintf("BER for %d-QAM", M);
% SER_Title = sprintf("SER for %d-QAM", M);
x_axis = "Es/N0 (dB)";

legend_order = ["Golden Code (MLD)", "Uncoded (MLD)", "Mod Diversity", "ZF"];

myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1]);
xlim([0 24]);
% myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
% ylim([10^(-6) 1])