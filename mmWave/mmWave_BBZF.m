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
EsN0_dB = -20:5:20;
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
Tx_RF_USER = 4;


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
    desc_idx = desc_idx(1:Tx_RF_USER, :); % keeping only as much values as RF-chain numbers
    for i1=1:Nts % NOT SURE IF N_u or N_ts
            F_RF(:, Tx_RF_USER * (i1-1) + 1 : Tx_RF_USER * i1) = UPA_Tx(:, RayNumber*(i1-1)+desc_idx(:, i1));
    end
    H_eff = H * F_RF;
    % F_bb_zf = H' * inv(H*H');
    % F_BB_ZF = inv(H_eff' * H_eff) * H_eff';
    F_BB_ZF = H_eff' * inv(H_eff * H_eff');
    
    %% Normalization
    F = F_RF * F_BB_ZF;
    F = F/norm(F, 'fro')*sqrt(Nts);
    
    % tmp_BB_ZF_Gain = H*F*F'*H';
    % BB_ZF_Gain = diag(tmp_BB_ZF_Gain);

    tmp_gain_BB_ZF = H * F * F' * H';
    BB_ZF_HF_Gain = diag(tmp_gain_BB_ZF);
    BB_ZF_Gain = abs(diag(H*F)).^2;
    
    %% Full dimension ZF
    F_zf = H' * inv(H*H');
    
    for i2 = 1 : Nu;
        F_zf( : , i2) = F_zf( : , i2) / norm(F_zf( : , i2), 'fro');
    end
    
    H_eff_zf = H * F_zf;
    tmp_gain_zf = H_eff_zf * H_eff_zf';
    ZF_Gain = diag(tmp_gain_zf);

    %% Beam Steering
    F_BS = F_RF;
    for i1=1:Nts % NOT SURE IF N_u or N_ts
            % F_RF(:, Tx_RF_USER * (i1-1) + 1 : Tx_RF_USER * i1) = UPA_Tx(:, RayNumber*(i1-1)+desc_idx(:, i1));
            % HONESTLY NOT SURE IF IT'S VALID
            W_BS(:, Tx_RF_USER * (i1-1) + 1 : Tx_RF_USER * i1) = UPA_Rx(:, RayNumber*(i1-1)+desc_idx(:, i1));
    end
    % NOTE: NOT INTUITIVE BUT IT WORKS
    H_eff_BS = H * F_BS;
    tmp_gain_BS = H_eff_BS * H_eff_BS';
    BS_HF_Gain = diag(tmp_gain_BS);
    BS_Gain = abs(diag(H_eff_BS)).^2;

    %% simulation
    for SNR_idx = 1:length(EsN0)
        for user=1:Nu
            %% BB-ZF
            % tmp = eye(Nus) + ZF_Gain(user) / (sum(tmp_gain_zf(user, :), 'all') - ZF_Gain(user)+Nts/EsN0(SNR_idx)*eye(Nus));
            % tmp = eye(Nus) + BB_ZF_Gain(user) / (sum(tmp_BB_ZF_Gain(user, :), 'all') - BB_ZF_Gain(user) + Nts/EsN0(SNR_idx)*eye(Nus));
            % tmp = eye(Nus) + BB_ZF_Gain(user) / (BB_ZF_HF_Gain(user) - BB_ZF_Gain(user)+Nts/EsN0(SNR_idx)*eye(Nus));
            gain = H(user, :) * F(:, user) * F(:, user)' * H(user, :)';
            interference = 0;
            for i1=1:Nu
                if i1==user
                    continue
                end
                interference = interference + H(user, :) * F(:, i1) * F(:, i1)' * H(user, :)';
            end
            tmp = eye(Nus) + gain / (interference+Nts/EsN0(SNR_idx)*eye(Nus));
            rate = log2(det(tmp));
            DataRate(SNR_idx, 1) = DataRate(SNR_idx, 1) + rate;

            %% ZF
            tmp = eye(Nus) + ZF_Gain(user) / (sum(tmp_gain_zf(user, :), 'all') - ZF_Gain(user)+Nts/EsN0(SNR_idx)*eye(Nus));
            rate = log2(det(tmp));
            DataRate(SNR_idx, 2) = DataRate(SNR_idx, 2) + rate;

            %% BEAM STEERING
            % 현재 UPA_Rx가 모두 1로 나오는 문제 있음. 그거 수정하면 아래 코드도 수정 필요
            tmp = eye(Nus) + BS_Gain(user) / (BS_HF_Gain(user) - BS_Gain(user)+Nts/EsN0(SNR_idx)*eye(Nus));
            rate = log2(det(tmp));
            DataRate(SNR_idx, 3) = DataRate(SNR_idx, 3) + rate;
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
grid on
plot(EsN0_dB, DataRate(:,2), '.--', 'MarkerSize', 15);
plot(EsN0_dB, DataRate(:,3), '.--', 'MarkerSize', 15);