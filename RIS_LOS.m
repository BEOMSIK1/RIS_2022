clc, clear
%% Parameters
M = 64;    % BS ant.
N = 128;    % RIS ele.
N1 = 16;    % RIS x-axis
N2 = N/N1;  % RIS y-axis
U = 4;      % User
RA = 1;     % RX ant.
S = U*RA;   % Data stream
L1 = 4;     % BS-RIS paths
L2 = 4;     % RIS-UE paths
L3 = 4;     % BS-UE paths
scatter = 10;
snr = 0:3:30;
iter = 300;
%% SCM parameters
model_BR = SCM();
model_BR.n_path = L1;
model_BR.n_mray = scatter;
model_BR.tx_ant = [M, 1, 0.5, 0.5];
model_BR.rx_ant = [N1, N2, 0.5, 0.5];
model_BR.asd = 10;
model_BR.zsd = 10;
model_BR.asa = 10;
model_BR.zsa = 10;
model_BR.los = 0;
model_BR.fc = 28e9;
model_BR.fs = 20e6;

model_RU = SCM();
model_RU.n_path = L2;
model_RU.n_mray = scatter;
model_RU.tx_ant = [N1, N2, 0.5, 0.5];
model_RU.rx_ant = [U, 1, 10, 10];
model_RU.asd = 10;
model_RU.zsd = 10;
model_RU.asa = 10;
model_RU.zsa = 10;
model_RU.los = 0;
model_RU.fc = 28e9;
model_RU.fs = 20e6;

model_BU = SCM();
model_BU.n_path = L3;
model_BU.n_mray = scatter;
model_BU.tx_ant = [M, 1, 0.5, 0.5];
model_BU.rx_ant = [U, 1, 0.5, 0.5];
model_BU.asd = 10;
model_BU.zsd = 10;
model_BU.asa = 10;
model_BU.zsa = 10;
model_BU.los = 0;
model_BU.fc = 28e9;
model_BU.fs = 20e6;
%% init.
SR = zeros(1,length(snr));
BER = zeros(1,length(snr));
%% RIS sim.
for i = 1:iter
    %% channel gene.
    BR_temp = model_BR.FD_channel(1);
    H_BR(:,:) = BR_temp(1,1,:,:);
    for snr_i = 1:length(snr)
    end
end