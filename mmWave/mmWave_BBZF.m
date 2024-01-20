clc
clear all
%% Description %%
% 수업자료 17페이지: mmWave BB-ZF

%% Base Config
Nt = 64;

lambda = 1; % Wavelength
d = lambda / 2; % Inter-antenna spacing
k = 2 * pi / lambda; % Array response coefficient

%% Configuration
EsN0_dB = -40:4:0;
EsN0 = db2pow(EsN0_dB);

Nu = 4; % User
Nus = 1;% User Stream; Antenna per user
Nts = Nu * Nus; % total stream

iTotal = 10^3;
RayNumber = 10;

Tx_W = sqrt(Nt);
Tx_H = sqrt(Nt);
Rx_W = sqrt(Nus);
Rx_H = sqrt(Nus);

Rx_RF = 1;
Tx_RF = 1;


DataRate = zeros(length(EsN0), 3); % ZFBF, R-ZF BF, Beamsteering

%% Timer
FivePercent = ceil(iTotal/20);

for iteration = 1:iTotal
    if mod(iteration-100, FivePercent)==0
        tic
    end
    
    Nr=1;
    for i1 = 1 : Nu % Nu or Nts NOT SURE!!
        [H(Nr * (i1 - 1) + 1 :  Nr * i1, 1 : Nt), UPA_Tx(1 : Nt, RayNumber * (i1 - 1) + 1 :  RayNumber * i1), UPA_Rx(1 : Nr, RayNumber * (i1 - 1) + 1 :  RayNumber * i1), alpha( : , i1)]...
            = mmWave_channel_realization(RayNumber, d, k, Nt, Nr, Tx_W, Tx_H, Rx_W, Rx_H);
    end
    % alpha = ray x user
    [~, desc_idx] = sort(abs(alpha), 1, 'descend');
    desc_idx = desc_idx(1:Tx_RF, :); % keeping only as much values as RF-chain numbers
    for i1=1:Nts % NOT SURE IF N_u or N_ts
            F_RF(:, Tx_RF * (i1-1) + 1 : Tx_RF * i1) = UPA_Tx(:, RayNumber*(i1-1)+desc_idx(:, i1));
    end
    H_eff = H * F_RF;
    % F_bb_zf = H' * inv(H*H');
    F_bb_zf = inv(H_eff' * H_eff) * H_eff';

    %% Normalization
    % for i2 = 1 : Nu
    %     F_zf( : , i2) = F_zf( : , i2) / norm(F_zf( : , i2), 'fro');
    % end
    
    H_eff_zf = H * F_zf;
    tmp_gain_zf = H_eff_zf * H_eff_zf';
    ZF_Gain = diag(tmp_gain_zf);

    %% Channel Realization
    for SNR_idx = 1:length(EsN0)
        F_mmse = H' * inv(H * H' + Nts / EsN0(SNR_idx) * eye(Nts));
        for i2=1:Nu
            F_mmse( : , i2) = F_mmse( : , i2) / norm(F_mmse( : , i2), 'fro');
        end
        H_eff_mmse = H * F_mmse;
        tmp_gain_mmse = H_eff_mmse * H_eff_mmse';
        MMSE_Gain = diag(tmp_gain_mmse);

        for user=1:Nu
            tmp = eye(Nus) + ZF_Gain(user) / (sum(tmp_gain_zf(user, :), 'all') - ZF_Gain(user)+Nts/EsN0(SNR_idx)*eye(Nus));
            rate = log2(det(tmp));
            DataRate(SNR_idx, 1) = DataRate(SNR_idx, 1) + rate;

            tmp = eye(Nus) + MMSE_Gain(user) / (sum(tmp_gain_mmse(user, :), 'all') - MMSE_Gain(user)+Nts/EsN0(SNR_idx)*eye(Nus));
            rate = log2(det(tmp));
            DataRate(SNR_idx, 2) = DataRate(SNR_idx, 2) + rate;
        end
    end

    if mod(iteration-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iTotal-iteration)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iteration/iTotal*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end

DataRate = DataRate / iTotal;
figure;
plot(EsN0_dB, DataRate(:,1), '.--', 'MarkerSize', 15);
hold on
plot(EsN0_dB, DataRate(:,2), '.--', 'MarkerSize', 15);


function [H, Tx_UPA, Rx_UPA, alpha] = mmWave_channel_realization_MU(Ray_number, d, k, Nt, Nus, Tx_W, Tx_H, Rx_W, Rx_H, Nu)
    %% mmWave Channel Realization
    %% Initialization
    Tx_azimuth_angle = zeros(1, Ray_number);
    Tx_elevation_angle = zeros(1, Ray_number);
    Rx_azimuth_angle = zeros(1, Ray_number);
    Rx_elevation_angle = zeros(1, Ray_number);
    Tx_UPA = zeros(Tx_W * Tx_H, Ray_number);
    Rx_UPA = zeros(Rx_W * Rx_H, Ray_number, Nu);
    %% L-ray channel realization
    for ray = 1 : Ray_number;
        %% AoA & AoD generation
        Tx_azimuth_angle(ray) = 2 * pi * rand;
        Tx_elevation_angle(ray) = 2 * pi * rand;
        Rx_azimuth_angle(ray) = 2 * pi * rand;
        Rx_elevation_angle(ray) = 2 * pi * rand;
    %     %% Test
    %     for i1 = 0 : Tx_W - 1
    %         temp(1, i1 * Tx_W + 1 : (i1 + 1) * Tx_W) = (i1 * 1i) + (0 : 1 : Tx_H - 1);
    %     end
        %% Transmitter UPA generation
        for i1 = 0 : Tx_W - 1
            Tx_UPA(i1 * Tx_W + 1 : (i1 + 1) * Tx_W, ray) = exp(1i * k * d * (i1 * sin(Tx_azimuth_angle(ray)) * sin(Tx_elevation_angle(ray))...
                 + (0 : 1 : Tx_H - 1) * cos(Tx_elevation_angle(ray))));
        end
    %     %% Test (using Kronecker product)
    %     temp_1 = exp(1i * k * d * cos(Tx_elevation_angle(ray)) * (0 : 1 : Tx_H - 1));
    %     temp_2 = exp(1i * k * d * sin(Tx_azimuth_angle(ray)) * sin(Tx_elevation_angle(ray)) * (0 : 1 : Tx_W - 1));
    %     Tx_UPA_temp( : , ray) = kron(temp_1, temp_2);
        %% Receiver UPA generation
        for user=1:Nu
            for i1 = 0 : Rx_W - 1
                Rx_UPA(i1 * Rx_W + 1 : (i1 + 1) * Rx_W, ray, user) = exp(1i * k * d * (i1 * sin(Rx_azimuth_angle(ray)) * sin(Rx_elevation_angle(ray))...
                     + (0 : 1 : Rx_H - 1) * cos(Rx_elevation_angle(ray))));
            end
        end
    end % ray iteration
    % Tx_UPA_temp = 1 / sqrt(Nt) * Tx_UPA;
    Tx_UPA = 1 / sqrt(Nt) * Tx_UPA; % At matrix
    Rx_UPA = 1 / sqrt(Nus) * Rx_UPA; % Ar matrix
    
    %% Summing each propagation paths
    H = zeros(Nus * Nu, Nt); % Initialization of channel
    alpha = zeros(Ray_number, Nu); % initialization of complex gains    
    for user=1:Nu
        for ray = 1 : Ray_number;
            alpha(ray, user) = (randn + randn * 1i) / sqrt(2);
            temp_H = sqrt(Nt * Nus / Ray_number) * alpha(ray, user) * Rx_UPA( : , ray, user) * Tx_UPA( : , ray)';
            H([Nus*(user-1)+1:Nus*user], :) = H([Nus*(user-1)+1:Nus*user], :) + temp_H;
        end
    end
    % alpha = alpha.';
    %% Test (using matrix product)
    % D = sqrt(Nt * Nr / Ray_number) * diag(alpha);
    % H2 = Rx_UPA * D * Tx_UPA';
end