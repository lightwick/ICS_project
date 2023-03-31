close all
clc
clear all

%% Basic Settings
addpath("../../tools/");

%% Environment variable
sim_num = 2;
M = 4;

EsN0_dB = 0:4:24;
EsN0 = db2pow(EsN0_dB);
EbN0 = EsN0 / log2(M);
EbN0_dB = pow2db(EbN0);

iteration = 10^2;

%% Simulate
BER = zeros(sim_num, length(EbN0_dB));

% [BER_tmp, ~] = simulate_golden_code(M, iteration, EbN0_dB);
% BER(1,:) = BER_tmp';

eta = 2;
n = 4;
TimeFrame = 3;
BER_tmp = simulate_modulation_diversity_reduced(eta, n, TimeFrame, iteration, EbN0_dB);
BER(1,:) = BER_tmp';

% n = 4;
% [BER_tmp, ~] = simulate_mld(M, n, iteration, EbN0_dB);
% BER(2, :) = BER_tmp';
% 
% n = 8;
% [BER_tmp, ~] = simulate_mld(M, n, iteration, EbN0_dB);
% BER(3, :) = BER_tmp';
% eta = 2;
% n = 2;
% BER_tmp = simulate_modulation_diversity(eta, n, iteration, EbN0_dB);
% BER(3,:) = BER_tmp';

%% Plot
% BER_Title = sprintf("BER for %d-QAM", M);
% SER_Title = sprintf("SER for %d-QAM", M);
title = "MIMO \zeta=4";
x_axis = "Es/N0 (dB)";

legend_order = ["Reduced Mod Div","Uncoded 4x4 (ML)", "Uncoded 8x8 (ML)"];

myplot(EsN0_dB, BER, title, x_axis, "BER", legend_order);
ylim([10^(-6) 1]);
xlim([0 24]);
% myplot(EsN0_dB, SER, SER_Title, x_axis, "SER", legend_order);
% ylim([10^(-6) 1])