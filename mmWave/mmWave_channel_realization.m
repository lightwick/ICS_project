function [H, Tx_UPA, Rx_UPA, alpha] = mmWave_channel_realization(Ray_number, d, k, Nt, Nr, Tx_W, Tx_H, Rx_W, Rx_H)
%% mmWave Channel Realization
%% Initialization
Tx_azimuth_angle = zeros(1, Ray_number);
Tx_elevation_angle = zeros(1, Ray_number);
Rx_azimuth_angle = zeros(1, Ray_number);
Rx_elevation_angle = zeros(1, Ray_number);
Tx_UPA = zeros(Tx_W * Tx_H, Ray_number);
Rx_UPA = zeros(Rx_W * Rx_H, Ray_number);
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
    for i1 = 0 : Rx_W - 1
        Rx_UPA(i1 * Rx_W + 1 : (i1 + 1) * Rx_W, ray) = exp(1i * k * d * (i1 * sin(Rx_azimuth_angle(ray)) * sin(Rx_elevation_angle(ray))...
             + (0 : 1 : Rx_H - 1) * cos(Rx_elevation_angle(ray))));
    end
end % ray iteration
% Tx_UPA_temp = 1 / sqrt(Nt) * Tx_UPA;
Tx_UPA = 1 / sqrt(Nt) * Tx_UPA; % At matrix
Rx_UPA = 1 / sqrt(Nr) * Rx_UPA; % Ar matrix

%% Summing each propagation paths
H = zeros(Nr, Nt); % Initialization of channel
alpha = zeros(1, Ray_number); % initialization of complex gains
for ray = 1 : Ray_number;
    alpha(ray) = (randn + randn * 1i) / sqrt(2);
    temp_H = sqrt(Nt * Nr / Ray_number) * alpha(ray) * Rx_UPA( : , ray) * Tx_UPA( : , ray)';
    H = H + temp_H;
end
alpha = alpha.';
%% Test (using matrix product)
% D = sqrt(Nt * Nr / Ray_number) * diag(alpha);
% H2 = Rx_UPA * D * Tx_UPA';
return