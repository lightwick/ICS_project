%% TODO
% check if input configuration is viable with the code
% condition 1: in M=16, n_t<=8

%% Basic Configuration
addpath('../tools');
M = 16;
EsN0_dB = 0:10;

EsN0_config % configuration script in tools

Nt = 8;
Nr = Nt;

iteration = 1;

NormalizationFactor = sqrt(2/3*(M-1) * Nt);

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
    CodebookIndex = floor(TransmitterDecimal/a)+1;
    CodewordIndex = mod(TransmitterDecimal, a)+1;

    SignalSequence = randi([0 M-1], 2, 1);
    SignalBinary = de2bi(SignalSequence, log2(M), 'left-msb');
    SymbolSequence = qammod(SignalSequence, M) / NormalizationFactor;
    
    
    %% Coding
    STBC = [SymbolSequence'; -conj(SymbolSequence(2,1)) conj(SymbolSequence(1,1))].';
    STBC_SM = zeros(2,Nt);
    
    STBC_SM(:,codebook{CodebookIndex}{CodewordIndex}) = STBC * theta(CodebookIndex);
    
    H = (randn(Nt, Nr) + 1j * randn(Nt, Nr)) ./ sqrt(2); % Receiver x Transmitter
    Noise = (randn(2, Nr) + 1j * randn(2, Nr)) / sqrt(2);
    ReceivedSignal = STBC_SM * H + Noise;
    %% Decoding
    


    if mod(iTotal-100, FivePercent)==0
        ElapsedTime = toc;
        EstimatedTime = (iteration-iTotal)*ElapsedTime;
        disp(sprintf("%d%%, estimated wait time %d minutes %d seconds", round(iTotal/iteration*100), floor(EstimatedTime/60), floor(mod(EstimatedTime, 60))))
    end
end