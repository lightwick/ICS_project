close all
clear all
dbstop if error
% dbstop if warning

addpath('../tools/')

%% Environment Varible
M = 4

NormalizationFactor = sqrt(2/3*(M-1));
colors = ["#0000FF", "#FF0000", "#4DBEEE", "#D95319", "#77AC30", "#EDB120", "#7E2F8E"];

iTotal = 10^4;

%% (A) Capacity CDF
Nr = 2;
Nt = 2;

EsN0_dB = 10;
EsN0 = db2pow(EsN0_dB);

[EntireCapacity_CU, ErgodicCapacity_CU, OutageCapacity_CU] = getCapcity_CU(Nr, Nt, EsN0_dB, iTotal)
[EntireCapacity_CK, ErgodicCapacity_CK, OutageCapacity_CK] = getCapcity_CK(Nr, Nt, EsN0_dB, iTotal)

figure();
hold on
cdfplot(EntireCapacity_CU);
cdfplot(EntireCapacity_CK);

plot([0 OutageCapacity_CU], [0.1 0.1], "k");
plot([0 OutageCapacity_CK], [0.1 0.1], "k");

plot([OutageCapacity_CU OutageCapacity_CU], [0 0.1], "k");
plot([OutageCapacity_CK OutageCapacity_CK], [0 0.1], "k");

text(OutageCapacity_CK, 0.1, '\leftarrow10% outage capacity');

text(OutageCapacity_CU, 0.1+0.025, '\downarrow', 'HorizontalAlignment', 'center');
text(OutageCapacity_CU, 0.1+0.07, '10% outage capacity', 'HorizontalAlignment', 'center');

ErgodicLine_CU = xline(ErgodicCapacity_CU, '-', {'Ergodic Capacity'});
ErgodicLine_CK = xline(ErgodicCapacity_CK, '-', {'Ergodic Capacity'});

xlim([1 10]);
grid on
xlabel("Rate (bps/Hz)");
ylabel("CDF");
title("CDF of Channel Rate");
legend({"Channel Unknown","Channel Known"}, 'Location','northwest')
% plotCapacityCDF(Nr, Nt, EsN0);

%% (B)
EsN0_dB = 0:5:20;
MarkerSize = 15;
Marker = '.--';
Nr = [1 2 1 2 4];
Nt = [1 1 2 2 4];

%% Ergodic Capacity
Ergodic_CU = figure();
xlabel("Es/N0 (dB)"); ylabel("Ergodic Capacity (bps/Hz)");
title("Ergodic Capacity (Channel Unknown)");
grid on
ylim([0 25]);

Outage_CU = figure();
xlabel("Es/N0 (dB)"); ylabel("Ergodic Capacity (bps/Hz)");
title("10% Outage Capacity (Channel Unknown)");
grid on
ylim([0 25]);

Ergodic_CK = figure();
xlabel("Es/N0 (dB)"); ylabel("Ergodic Capacity (bps/Hz)");
title("Ergodic Capacity (Channel Known)");
grid on
ylim([0 25]);

Outage_CK = figure();
xlabel("Es/N0 (dB)"); ylabel("Ergodic Capacity (bps/Hz)");
title("10% Outage Capacity (Channel Known)");
grid on
ylim([0 25]);

for ii=1:length(Nr)
    [EntireCapacity_CU, ErgodicCapacity_CU, OutageCapacity_CU] = getCapcity_CU(Nr(ii), Nt(ii), EsN0_dB, iTotal);
    figure(Ergodic_CU); hold on;
    plot(EsN0_dB, ErgodicCapacity_CU, Marker, 'Color', colors(ii), 'MarkerSize', MarkerSize);
    figure(Outage_CU); hold on;
    plot(EsN0_dB, OutageCapacity_CU, Marker, 'Color', colors(ii), 'MarkerSize', MarkerSize);
    
    [EntireCapacity_CK, ErgodicCapacity_CK, OutageCapacity_CK] = getCapcity_CK(Nr(ii), Nt(ii), EsN0_dB, iTotal);
    figure(Ergodic_CK); hold on;
    plot(EsN0_dB, ErgodicCapacity_CK, Marker, 'Color', colors(ii), 'MarkerSize', MarkerSize);
    figure(Outage_CK); hold on;
    plot(EsN0_dB, OutageCapacity_CK, Marker, 'Color', colors(ii), 'MarkerSize', MarkerSize);
end

grid on
figure(Ergodic_CU);
legend({"N_T=1, N_R=1", "N_T=1, N_R=2", "N_T=2, N_R=1", "N_T=2, N_R=2", "N_T=4, N_R=4"}, 'Location','northwest');
figure(Outage_CU); 
legend({"N_T=1, N_R=1", "N_T=1, N_R=2", "N_T=2, N_R=1", "N_T=2, N_R=2", "N_T=4, N_R=4"}, 'Location','northwest');
figure(Ergodic_CK);
legend({"N_T=1, N_R=1", "N_T=1, N_R=2", "N_T=2, N_R=1", "N_T=2, N_R=2", "N_T=4, N_R=4"}, 'Location','northwest');
figure(Outage_CK);
legend({"N_T=1, N_R=1", "N_T=1, N_R=2", "N_T=2, N_R=1", "N_T=2, N_R=2", "N_T=4, N_R=4"}, 'Location','northwest');

%% (C)
EsN0_dB = 0:20;
Nt = 4;
Nr = 4;

Marker = '.-';
MarkerSize = 10;
[EntireCapacity_CU, ErgodicCapacity_CU, OutageCapacity_CU] = getCapcity_CU(Nr, Nt, EsN0_dB, iTotal);
[EntireCapacity_CK, ErgodicCapacity_CK, OutageCapacity_CK] = getCapcity_CK(Nr, Nt, EsN0_dB, iTotal);

figure();
hold on;
grid on;
plot(EsN0_dB, ErgodicCapacity_CU, Marker, 'Color', colors(1), 'MarkerSize', MarkerSize);
plot(EsN0_dB, ErgodicCapacity_CK, Marker, 'Color', colors(2), 'MarkerSize', MarkerSize);
legend({"Channel Unknown to Tx", "Channel Known to Tx"}, 'Location','northwest');
xlabel("Es/N0(dB)");
ylabel("Ergodic Capacity (bps/Hz)");
title("Ergodic Capacity by SNR");
ylim([0 25]);

figure();
hold on;
grid on;
plot(EsN0_dB, OutageCapacity_CU, Marker, 'Color', colors(1), 'MarkerSize', MarkerSize);
plot(EsN0_dB, OutageCapacity_CK, Marker, 'Color', colors(2), 'MarkerSize', MarkerSize);
legend({"Channel Unknown to Tx", "Channel Known to Tx"}, 'Location','northwest');
xlabel("Es/N0(dB)");
ylabel("Outage Capacity (bps/Hz)");
title("Outage Capacity by SNR")
ylim([0 25]);