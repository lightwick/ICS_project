close all
clear all
dbstop if error
% dbstop if warning

addpath('../tools/')

%% Environment Varible
M = 4

NormalizationFactor = sqrt(2/3*(M-1));
colors = ["#0000FF", "#FF0000", "#4DBEEE", "#D95319", "#77AC30", "#EDB120", "#7E2F8E"];

%% (A) Capacity CDF
Nr = 2;
Nt = 2;

EsN0_dB = 10;
EsN0 = db2pow(EsN0_dB);

plotCapacityCDF(Nr, Nt, EsN0);
xlim([1 10]);

%% (B)
EsN0_dB = 0:5:20;
MarkerSize = 15;
Marker = '.--';
Nr = [1 2 1 2 4];
Nt = [1 1 2 2 4];

%% Ergodic Capacity
figure();
hold on
for ii=1:length(Nr)
    plotErgodicCapcity(Nr(ii), Nt(ii), EsN0_dB, Marker, colors(ii), MarkerSize);
end
grid on

legend({"N_T=1, N_R=1", "N_T=1, N_R=2", "N_T=2, N_R=1", "N_T=2, N_R=2", "N_T=4, N_R=4"}, 'Location','northwest')
title("Ergodic Capacity");
xlabel("Es/N0 (dB)");
ylabel("Ergodic Capacity (bps/Hz)");

%% Outage Capacity
figure();
hold on
for ii=1:length(Nr)
    plotOutageCapcity(Nr(ii), Nt(ii), EsN0_dB, Marker, colors(ii), MarkerSize);
end
grid on
legend({"N_T=1, N_R=1", "N_T=1, N_R=2", "N_T=2, N_R=1", "N_T=2, N_R=2", "N_T=4, N_R=4"}, 'Location','northwest')
title("10% Outage Capacity");
xlabel("Es/N0 (dB)");
ylabel("10% Outage Capacity (bps/Hz)");

%% (C)
EsN0_dB = 0:20;
Nt = 4;
Nr = 4;
figure();
Marker = '-';
hold on;
plotErgodicCapcity(4, 4, EsN0_dB, Marker, colors(1), MarkerSize);
grid on;