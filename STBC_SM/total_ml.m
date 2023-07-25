%close all
clear all
clc

%% TODO
% check if input configuration is viable with the code
% condition 1: in M=16, n_t<=8

%% Current STATUS
% When Nt=2, no BER==0; meaning there is a problem when SM is added

%% DEBUGGING FEATURES
DEBUG = false;
NOISE = true;

%% User Defined Configuration
addpath('../tools');
M = 4;
EsN0_dB =  0:2:18;
if DEBUG==true
    EsN0_dB=0;
end

Nt = 4;
Nr = Nt;

iteration = 10^4;


%% Basic Configuration
EsN0_config % configuration script in tools

num = 1;
BEC = zeros(num, length(EsN0_dB));
TBEC = zeros(num, length(EsN0_dB));
SBEC = zeros(num, length(EsN0_dB));

NormalizationFactor = sqrt(2/3*(M-1) * 2); % Normally, this would be 'NormalizationFactor = sqrt(2/3*(M-1) * Nt);', but Spatial Modulation is added here

%% STBC-SM
map = codebook_gen();
codebook=map(Nt);

c = Nt*(Nt-1)/2;
while bitand(c, c-1) % a way to check if a number is the power of 2
    c=c-1;
end
n = size(codebook,1); % number of codebooks
a = size(codebook{1},2); % number of codewords per codebook

last_a = c-a*(n-1); % Has no purpose in this version of implementation(n_t=3,4,8)

if Nt<=4
    switch M
        case 2
            theta(2) = 1.57
        case 4
            theta(2) = 0.61
        case 16
            theta(2) = 9.05
        case 64
            theta(2) = 8.23
        otherwise
            error('M value needs to be among the following: 2, 4, 16, 64')
    end
else
    for k=2:n
        switch M
            case 2
                theta(k)=(k-1)*pi/(n);
            case {4, 16}
                theta(k)=(k-1)*pi/(2*n);
            otherwise
                error('M value needs to be 2,4, or 16 in n_t>4 case')
        end
    end
end

%% Modulated Candidate Generation
% Candidates = zeros(M^2, 2);
% 
%  for ii = 0:M^2-1
%      for jj = 1:2
%          Candidates(ii+1,jj) = mod(floor(ii/M^(2-jj)),M);
%      end
%  end
%  Candidates = qammod(Candidates', M) / NormalizationFactor; % Each column is a Candidate
Candidates = qammod(0:M-1, M) / NormalizationFactor;

%% Main simulation process
% Timer
FivePercent = ceil(iteration/20);

for iTotal = 1:iteration
    if mod(iTotal-100, FivePercent)==0
        tic
    end
    %% Bit Generation
    % Transmitter Defined Bits
    TransmitterDecimal = randi([0 c-1]);
    if Nt~=2
        TransmitterBinary = de2bi(TransmitterDecimal, log2(c), 'left-msb');
    end

    CodebookIndex = floor(TransmitterDecimal/a)+1;
    CodewordIndex = mod(TransmitterDecimal, a)+1;

    SignalSequence = randi([0 M-1], 2, 1);
    SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
    SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
    
    %% Coding
    STBC = [SymbolSequence.'; -conj(SymbolSequence(2,1)) conj(SymbolSequence(1,1))];
    STBC_SM = zeros(2,Nt);
    
    STBC_SM(:, codebook{CodebookIndex}{CodewordIndex}) = STBC * exp(1j*theta(CodebookIndex));
    
    H = (randn(Nt, Nr) + 1j * randn(Nt, Nr)) ./ sqrt(2); % Receiver x Transmitter
    if NOISE==false
        Noise = zeros(2,Nr);
    else
        Noise = (randn(2, Nr) + 1j * randn(2, Nr)) / sqrt(2);
    end

    H_l = zeros(Nr*2, 2, c);

    %% Create H_l23
    for book=1:n
        assert(a==last_a, 'The current implementation needs a==last_a. Check configuration or improve code')
        for word=1:a % This should be fixed if a~=last_a
            H_t = H.'; % transpose so Nr x Nt
            tmp_1 = H_t(:, codebook{book}{word}) * exp(1j*theta(book));
            tmp_2 = conj(tmp_1(:, [2, 1]));
            tmp_2(:, 2) = -tmp_2(:, 2);

            H_l(1:2:Nr*2, :, word+(book-1)*a) = tmp_1;
            H_l(2:2:Nr*2, :, word+(book-1)*a) = tmp_2;
        end
    end


    for SNR_idx = 1 : length(EsN0)
        ReceivedSignal = STBC_SM * H + Noise / sqrt(EsN0(SNR_idx)); % shape of 2xNr
        y = reshape(ReceivedSignal, [], 1);
        y(2:2:size(y,1), 1) = conj(y(2:2:size(y,1), 1));

        % SM_Candidates = zeros(2, Nt, c*M^2);
        % FroNorm = zeros(c*M^2,1);
        index = 1;

        Candidate_hi = zeros(4, M^2*c);
        for book=1:n
            for word=1:a
                for x1_candidate=1:M
                    for x2_candidate=1:M
                        x1 = Candidates(x1_candidate);
                        x2 = Candidates(x2_candidate);
                        STBC_Candidate = [x1 x2; -conj(x2) conj(x1)];
                        % SM_Candidates(:, codebook{book}{word}, book*a+word) = STBC_Candidate * exp(1j*theta(book));
                        SM_Candidates = zeros(2, Nt);
                        SM_Candidates(:, codebook{book}{word}) = STBC_Candidate * exp(1j*theta(book));
                        Receive_Candidate = SM_Candidates*H;
                        RC = reshape(Receive_Candidate, [], 1);
                        RC(2:2:size(RC,1), 1) = conj(RC(2:2:size(RC,1), 1));
                        
                        FroNorm(index) = norm(y-RC,'fro')^2;
                        Candidate_hi(:, index) = [book; word; x1_candidate; x2_candidate];
                        index = index+1;
                    end
                end
            end
        end

        [~, detected_index] = min(FroNorm, [], 'all');
        detected_x1 = Candidate_hi(3, detected_index)-1;
        detected_x2 = Candidate_hi(4, detected_index)-1;
        detected_l = (Candidate_hi(1, detected_index)-1)*a + Candidate_hi(2, detected_index) - 1;

        DetectedSequence = [detected_x1; detected_x2];
        DetectedBinary = de2bi(DetectedSequence, log2(M), 'left-msb');
        DetectedTransmitterBinary = de2bi(detected_l, log2(c), 'left-msb');
        
        TransmitError = sum((TransmitterBinary~=DetectedTransmitterBinary), 'all');
        SignalError = sum((SignalBinary~=DetectedBinary), 'all');
        ErrorCount = TransmitError + SignalError;

        % if DEBUG==true
        %     c_tmp = CodewordIndex+(CodebookIndex-1)*a;
        %     assert(isequal(H_l(:, 1, c_tmp)*Candidates(SignalSequence(1,1)+1) + H_l(:, 2, c_tmp)*Candidates(SignalSequence(2,1)+1),y), 'non match');
        % end
        % if ErrorCount~=0
        %     disp("err")
        % end
        
        BEC(1, SNR_idx) = BEC(1, SNR_idx) + ErrorCount; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        TBEC(1, SNR_idx) = TBEC(1, SNR_idx) + TransmitError; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        SBEC(1, SNR_idx) = SBEC(1, SNR_idx) + SignalError; %bitcount(bitxor([detected_x1; detected_x2], SignalSequence)));
        
        % if ErrorCount~=0
        %     disp(sum((SignalBinary~=DetectedBinary), 'all'))
        %     disp(sum((TransmitterBinary~=DetectedTransmitterBinary), 'all'))
        % end
    end

    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end
TotalTransmitBits = (2 * log2(M) + log2(c))* iteration;

BER(1,:) = BEC(1,:)/TotalTransmitBits;
TBER(1,:) = TBEC(1,:)/(log2(c)*iteration);
SBER(1,:) = SBEC(1,:)/(2*log2(M)*iteration);

%% Plotting
BER_Title = sprintf("BER for %d-QAM", M);
x_axis = "SNR (dB)";
% legend_order = ["STBC-SM, n_T=8, 16-QAM, 6 bits/s/Hz"];
legend_order = ["foo"];
myplot(EsN0_dB, BER, BER_Title, x_axis, "BER", legend_order);
ylim([10^(-6) 1])