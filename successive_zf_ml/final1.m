close all
clear
clc

Antenna_T = 4;
Antenna_R = 4;
M=4;
NumberIteration = 10^5;

Normalization_power = sqrt(2/3 * (M-1));

Eb = 1;
EbN0_dB = 0 : 5 : 25;
EbN0 = db2pow(EbN0_dB);
EsN0 = EbN0 .* log2(M);

LengthBitSequence = Antenna_T * log2(M) ;

ErrorCount_ZF_1 = zeros(1, length(EbN0_dB));
ErrorCount_ML = zeros(1, length(EbN0_dB));

ErrorCount_ZF_S_1 = zeros(1, length(EbN0_dB));
ErrorCount_ZF_S_2 = zeros(1, length(EbN0_dB));
ErrorCount_ZF_S_3 = zeros(1, length(EbN0_dB));
ErrorCount_ZF_S_4 = zeros(1, length(EbN0_dB));
ErrorCount_Symbol_S = zeros(1, length(EbN0_dB));

ErrorCount_Bit_S = zeros(1, length(EbN0_dB));
ErrorCount_ZF1 = zeros(1, length(EbN0_dB));
ErrorCount_ML1 = zeros(1, length(EbN0_dB));

MinQammod_ML_all = zeros(Antenna_T,1);

QammodArray = zeros(1,1);

Detection_Symbol = zeros(Antenna_T,1);

% ML-대조군
Candidate1 = zeros(M^Antenna_T,Antenna_T);
CandidateSample1 = 0 : 1 :  M^Antenna_T - 1;
CandidateStr1 = dec2base(CandidateSample1.',M,Antenna_T);

for indx_row = 1 : M^Antenna_T
    for indx_col = 1 : Antenna_T
        if (M == 4)
            Candidate1(indx_row, indx_col) = str2double(CandidateStr1(indx_row, indx_col));
        end
        if (M == 16)
            Candidate1(indx_row, indx_col) = hex2dec(CandidateStr1(indx_row, indx_col));
        end
    end
end
QammodArray1 = qammod(Candidate1,M,'InputType','integer').* (1 ./ sqrt(2/3 .* (M-1))) ;

% ML-Simulation
Candidate = zeros(M^(Antenna_T-1),Antenna_T-1);
CandidateSample = 0 : 1 :  M^(Antenna_T-1) - 1;
CandidateStr = dec2base(CandidateSample.',M,Antenna_T-1);

for indx_row = 1 : M^(Antenna_T-1)
    for indx_col = 1 : Antenna_T-1
        if (M == 4)
            Candidate(indx_row, indx_col) = str2double(CandidateStr(indx_row, indx_col));
        end
        if (M == 16)
            Candidate(indx_row, indx_col) = hex2dec(CandidateStr(indx_row, indx_col));
        end
    end
end
QammodArray_Pre = qammod(Candidate,M,'InputType','integer').* (1 ./ sqrt(2/3 .* (M-1))) ;

for iTotal = 1 : NumberIteration
    tic
     
    Signal = randi([0 M-1],Antenna_T,1);
    Bit= de2bi(Signal,log2(M),'left-msb');
    S= qammod(Signal, M) / Normalization_power ;
    N = 1/sqrt(2)*(randn(Antenna_R,1) + 1j*randn(Antenna_R,1)); % N (n) Generation
    H = 1/sqrt(2)*randn(Antenna_R, Antenna_T*1) + 1j*randn(Antenna_R, Antenna_T*1); % Channel (h) Generation
    
    Signal1 = randi([0 M-1],Antenna_T,1);
    S1= qammod(Signal1, M) / Normalization_power ;
    Bit1 = de2bi(Signal1,log2(M),'left-msb');
    N1 = 1/sqrt(2)*(randn(Antenna_R,1) + 1j*randn(Antenna_R,1)); 
    H1 = 1/sqrt(2)*randn(Antenna_R, Antenna_T) + 1j*randn(Antenna_R, Antenna_T);
    
    %대조군----------------------------------------------------------------
    for indx_EbN0 = 1 : length(EbN0)
        
             y1 = sqrt(EsN0(indx_EbN0)/ Antenna_T) * H * S + N;
       
               % ZF------------
                w_ZF1 = sqrt(Antenna_T / EsN0(indx_EbN0)) * pinv(H);
                z_ZF1 = w_ZF1 * y1;
       
               % ML--------------
                Min_1_1 = norm((y1 - sqrt(EsN0(indx_EbN0) / Antenna_T) * H * QammodArray1(1,:).'),'fro');
                MinQammod_ML1 = QammodArray1(1,:);
 
                for indx_ML = 1 : M^Antenna_T
            
                  Comparison_ML1 = norm((y1 - sqrt(EsN0(indx_EbN0) / Antenna_T) * H * QammodArray1(indx_ML,:).'),'fro');
                    if (Comparison_ML1 <= Min_1_1)
                      Min_1_1 = Comparison_ML1;
                     MinQammod_ML1 = QammodArray1(indx_ML,:).';
                    end
                end
                DetectionSignal_ZF1 = qamdemod((z_ZF1 * Normalization_power) , M);
                DetectionSignal_ML1 = qamdemod(MinQammod_ML1 * Normalization_power, M);
                
                DetectionBitSequence_ZF1 = de2bi(DetectionSignal_ZF1, log2(M),'left-msb');
                DetectionBitSequence_ML1 = de2bi(DetectionSignal_ML1, log2(M),'left-msb');
        
                ErrorCount_Tmp_ZF1 = sum(DetectionBitSequence_ZF1 ~= Bit,'all'); 
                ErrorCount_Tmp_ML1 = sum(DetectionBitSequence_ML1 ~= Bit,'all'); 
                
                ErrorCount_ZF1(1, indx_EbN0) = ErrorCount_ZF1(1, indx_EbN0) + ErrorCount_Tmp_ZF1;
                ErrorCount_ML1(1, indx_EbN0) = ErrorCount_ML1(1, indx_EbN0) + ErrorCount_Tmp_ML1;
                
                
    end
    
    
    %Simulation------------------------------------------------------------
    %Simulation Step 1
    for indx_EbN0 = 1 : length(EbN0)
        
            y(:,1) = sqrt(EsN0(indx_EbN0)/ Antenna_T) * H( : ,(1 : Antenna_T) ) * S( : ,1) + N( : ,1);
            
            %ZF--------------
            w_ZF( : ,( 1 : Antenna_T)) = sqrt(Antenna_T / EsN0(indx_EbN0)) * pinv(H( : ,( 1 : Antenna_T) ));
            z_ZF(:,1) = w_ZF( : ,( 1 : Antenna_T)) * y(:,1);  
            
            %Detection
            DetectionSignal_ZF = qamdemod((z_ZF * Normalization_power) , M); 
            DetectionBitSequence_ZF = de2bi(DetectionSignal_ZF, log2(M),'left-msb');
                      
            %SER
            ErrorCount_Tmp_ZF_S_1 = sum(DetectionSignal_ZF(1,1) ~= Signal(1,1),'all'); % Error Count
            ErrorCount_ZF_S_1(1, indx_EbN0) = ErrorCount_ZF_S_1(1, indx_EbN0) + ErrorCount_Tmp_ZF_S_1;
            ErrorCount_Tmp_ZF_S_2 = sum(DetectionSignal_ZF(2,1) ~= Signal(2,1),'all'); 
            ErrorCount_ZF_S_2(1, indx_EbN0) = ErrorCount_ZF_S_2(1, indx_EbN0) + ErrorCount_Tmp_ZF_S_2;
            ErrorCount_Tmp_ZF_S_3 = sum(DetectionSignal_ZF(3,1) ~= Signal(3,1),'all'); 
            ErrorCount_ZF_S_3(1, indx_EbN0) = ErrorCount_ZF_S_3(1, indx_EbN0) + ErrorCount_Tmp_ZF_S_3;
            ErrorCount_Tmp_ZF_S_4 = sum(DetectionSignal_ZF(4,1) ~= Signal(4,1),'all'); 
            ErrorCount_ZF_S_4(1, indx_EbN0) = ErrorCount_ZF_S_4(1, indx_EbN0) + ErrorCount_Tmp_ZF_S_4;

            %Error min (ZF)
            Error_1 = sum(ErrorCount_ZF_S_1);
            Error_2 = sum(ErrorCount_ZF_S_2);
            Error_3 = sum(ErrorCount_ZF_S_3);
            Error_4 = sum(ErrorCount_ZF_S_4);
            

    end

   %BEST Antenna
            Error = [Error_1 Error_2 Error_3 Error_4];
            
            Best_Antenna = find(Error == min(Error));
            
            if numel(Best_Antenna) > 1
                Best_Antenna = Best_Antenna(1,1);         
            end


   %Next Iteration 을 위한 초기화
   ErrorCount_ZF_S_1 = zeros(1,indx_EbN0);
   ErrorCount_ZF_S_2 = zeros(1,indx_EbN0);
   ErrorCount_ZF_S_3 = zeros(1,indx_EbN0);
   ErrorCount_ZF_S_4 = zeros(1,indx_EbN0);

   %ML QammodArray 배열에 BestAntenna Symbol값 열 삽입
   Best_Antenna_S=repmat(S(Best_Antenna,1),[M^(Antenna_T-1) 1]);
   insert_row = Best_Antenna;
   if insert_row == 1
        QammodArray = horzcat(Best_Antenna_S, QammodArray_Pre);
   elseif insert_row == size(QammodArray_Pre, 2)+1
        QammodArray = horzcat(QammodArray_Pre, Best_Antenna_S);
   elseif insert_row > 1 && insert_row <= size(QammodArray_Pre, 2)
        QammodArray = horzcat(QammodArray_Pre(:,1:insert_row-1), Best_Antenna_S, QammodArray_Pre(:,insert_row:end));
   end

   %Simulation Step 2
   for indx_EbN0 = 1 : length(EbN0)
        
            y2(:,1) = sqrt(EsN0(indx_EbN0)/ (Antenna_T)) * H( : ,(1 : Antenna_T) ) * S( : ,1) + N( : ,1);
            
            % ML--------------
            Min_1 = norm((y2(:,1) - sqrt(EsN0(indx_EbN0) / (Antenna_T)) * H( : ,(1 : Antenna_T) ) * QammodArray(1,:).'),'fro');
            MinQammod_ML = QammodArray(1,:);
 
            for indx_ML = 1 : M^(Antenna_T-1)
            
                Comparison_ML = norm((y2(:,1) - sqrt(EsN0(indx_EbN0) / (Antenna_T)) * H( : ,(1 : Antenna_T) ) * QammodArray(indx_ML,:).'),'fro');
                if (Comparison_ML <= Min_1)
                   Min_1 = Comparison_ML;
                   MinQammod_ML = QammodArray(indx_ML,:).';
                end
                
                MinQammod_ML_all(:,1) = MinQammod_ML;
                
            end  

            ML_qamdemod = qamdemod(MinQammod_ML_all * Normalization_power, M);
            DetectionBitSequence_ML = de2bi(ML_qamdemod, log2(M),'left-msb');
            
            ErrorCount_Tmp_Bit_S = sum(DetectionBitSequence_ML ~= Bit,'all');
            ErrorCount_Bit_S(1, indx_EbN0) = ErrorCount_Bit_S(1, indx_EbN0) + ErrorCount_Tmp_Bit_S;

   end
   
   
  
    t = toc;
    if mod(iTotal, NumberIteration * 0.01) == 0
        fprintf('%5.1f초 남았습니다.(%4.1f%% 수행)\n', t * (NumberIteration - iTotal), iTotal * 100 / NumberIteration);
    end
end
    


SER_Simulation_ZF_1 = ErrorCount_ZF_S_1 / (NumberIteration);
SER_Simulation_ZF_2 = ErrorCount_ZF_S_2 / (NumberIteration);
SER_Simulation_ZF_3 = ErrorCount_ZF_S_3 / (NumberIteration);
SER_Simulation_ZF_4 = ErrorCount_ZF_S_4 / (NumberIteration);

BER_Simulation     = ErrorCount_Bit_S / (LengthBitSequence * NumberIteration);
BER_Simulation_ZF1 = ErrorCount_ZF1 / (LengthBitSequence * NumberIteration);
BER_Simulation_ML1 = ErrorCount_ML1 / (LengthBitSequence * NumberIteration);

ber1 = berfading(EbN0_dB,'qam',M,4);
% ber2 = berfading(EbN0_dB,'qam',M,1);

ML_4X3 = [0.102776500000000 0.0169305000000000 0.000711833333333333 1.21666666666667e-05 1.66666666666667e-07 0];
ML_4X4 = [0.155303750000000	0.0418362500000000	0.00254875000000000	5.00000000000000e-05	2.50000000000000e-06	0];
BER_Simulation_ML1 = ML_4X4;

% Plot_BER
figure()
semilogy(EbN0_dB, BER_Simulation, '-_b', 'LineWidth',2);
hold on
semilogy(EbN0_dB,BER_Simulation_ZF1, '--squarer', 'LineWidth',1);
semilogy(EbN0_dB, BER_Simulation_ML1, '--pentagramr', 'LineWidth',1);
semilogy(EbN0_dB, ber1, '--vk', 'LineWidth',1);
% semilogy(EbN0_dB, ber2, '--^k', 'LineWidth',1);

semilogy(EbN0_dB, ML_4X3, '--.', 'LineWidth',1, 'MarkerSize', 15);
% semilogy(EbN0_dB, 4x4_ML, '--.', 'LineWidth',1, 'MarkerSize', 15);
axis([0 25 10^-5 1]);
grid on

legend('Proposed Algorithm','ZF Receiver ','ML Receiver', 'Analysis (Diversity Order 4)', '4x3 ML');

%legend('Proposed Algorithm','ZF Receiver ','ML Receiver', 'Analysis (Diversity Order 4)');
xlabel('Eb/No [dB]');
ylabel('BER');
title(string(Antenna_T)+' x '+string(Antenna_R) + ' MIMO System with '+string(M)+' - QAM ');