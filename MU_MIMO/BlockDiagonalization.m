clear
clc

addpath('tools/');

%% BEGIN
iTotal = 10^4;

EsN0_dB = 0:5:20;
% EsN0_dB = 0;
EsN0 = db2pow(EsN0_dB);

Nt = 4;
K = 2; % Number of users
Nr = 2; % Number of Rx per user

ChannelSize = [Nr*K, Nt];
CapacitySum = zeros(length(EsN0), 1);

for iteration = 1:iTotal
    for SNR_idx = 1:length(EsN0)
        H_s = (randn(ChannelSize)+1j*randn(ChannelSize))/sqrt(2);
%         H_s = mod(randn(ChannelSize),5);

        % When debugging, this is should be the first thing to be
        % considered
        H_V_tilde = zeros(Nr, Nt - (Nr*K-Nr), K); % Unsure if this is right; This is assuming L_j_tilde = rank(H_j_tilde) = Nr*K-Nr
        for User=1:K
            StartIndex = (User-1)*Nr+1;
            EndIndex = User*Nr;

            H_j = H_s(StartIndex:EndIndex, :);
            H_j_tilde = H_s([1:StartIndex-1, EndIndex+1:size(H_s, 2)], :);
%             [~,~,V] = svd(H_s(((User-1)*Nr+1:User*Nr), :));
            L_j_tilde = rank(H_j_tilde);
            [~,~,V] = svd(H_j_tilde);
            V_tilde = V(:, (L_j_tilde+1:Nt));

            H_V_tilde(:,:,User) = H_j * V_tilde;
            % H_j_tilde * V_tilde should equal 0 matrix
            
            %% Power Distribution
            [~,S,~] = svd(H_V_tilde(:,:,User));

            EigenValues = diag(S).^2;

            p = 1;
            while true
                r = rank(H_V_tilde(:,:,User));
                PowerDistribution = Nr/(r-p+1)*(1+1/EsN0(SNR_idx)*sum(1./EigenValues, 'all')) - Nr/EsN0(SNR_idx)*(1./EigenValues);
                if PowerDistribution(r-p+1)>=0
                    break;
                else
                    p = p+1;
                    EigenValues(length(EigenValues))=[];
                end
            end

            Capacity = sum(log2(1+EsN0(SNR_idx)/Nr * PowerDistribution.*EigenValues), 'all');

            CapacitySum(SNR_idx) = CapacitySum(SNR_idx) + Capacity;
        end
%         H_s_prime = diag(H_V_tilde);
    end
end
ErgodicCapacity = CapacitySum / iTotal;

legend_ = sprintf('N_T=%d, K=%d, N_R=%d', Nt, K, Nr);

plot(EsN0_dB, ErgodicCapacity.');
title('Sum rate performance of BD');
xlabel('Es/N0 (dB)');
ylabel('Sum Rate [Bits/Transmission]');
legend(legend_);

%% Things I don't quite understand fully
% rank 개념

%% 불확실
% When doing water-filling, should Nt be changed with Nr?
% Meaning, does each user get the same power constraint?? (In the logical antenna set; meaning after the modulation)