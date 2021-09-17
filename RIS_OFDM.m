clc, clear
%% Parameters
fft_size = 64;
mod_type = 4;                     %1 - BPSK, 2 - QPSK, 4 - 16QAM, 6 - 64QAM, 8 - 256QAM
cp_size = fft_size / 4;
data_size = fft_size*mod_type;

M = 16;     % BS ant.
N = 64;     % RIS ele.
N1 = 8;     % RIS x-axis
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
model_BU.rx_ant = [U, 1, 10, 10];
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
    tic
    %% channel gene. (BS-RIS)
    BR_temp = model_BR.FD_channel(fft_size + cp_size);
    h_BR(:,:,:) = BR_temp(:,1,:,:);
    H_BR = fft(h_BR, fft_size, 1);
    %% channel gene. (RIS-UE)
    RU_temp = model_RU.FD_channel(fft_size + cp_size);
    h_RU(:,:,:) = RU_temp(:,1,:,:);
    H_RU = fft(h_RU, fft_size, 1);
    %% channel gene. (BS-UE)
    BU_temp = model_BU.FD_channel(fft_size + cp_size);
    h_BU(:,:,:) = BU_temp(:,1,:,:);
    H_BU = fft(h_BU, fft_size, 1);
    %% data gene.
    data = randi([0 1], S, data_size);
    sym = base_mod(data, mod_type);
    %% RIS matrix gene. (random)
    ang_ran = diag(exp(1j*2*pi*rand(1,N)));
    %% RIS matrix gene. (opt.)
    K_ = zeros(N,N);
    for k = 1:fft_size
        H_BR_(:,:) = H_BR(k,:,:);
        H_RU_(:,:) = H_RU(k,:,:);
        for i_ = 1:N
            for j_ = 1:N
                K_(i_,j_) = H_RU_(:,i_)'*H_RU_(:,j_)*H_BR_(j_,:)*H_BR_(j_,:)';
            end
        end
        K_ = K_ + K_;
    end
    K_ = K_/fft_size;
    [eig_v,~] = eig(K_);
    ang_opt = diag(exp(1j*2*pi*angle(eig_v(:,1))));
    %% iter(SNR)
    for snr_i = 1:length(snr)
        %% init
        n = 1;
        gamma_ran = zeros(1,fft_size);
        gamma_opt = zeros(1,fft_size);
        gamma_d = zeros(1,fft_size);
        
        for k = 1:fft_size
            H_BR_(:,:) = H_BR(k,:,:);
            H_RU_(:,:) = H_RU(k,:,:);
            H_BU_(:,:) = H_BU(k,:,:);
            %% eff. channel
            He_ran = H_RU_ * ang_ran * H_BR_;
            He_opt = H_RU_ * ang_opt * H_BR_;
            %% precoding (ZF)
            G_ran = He_ran' * inv(He_ran * He_ran');
            G_opt = He_opt' * inv(He_opt * He_opt');
            G_d = H_BU_' * inv(H_BU_ * H_BU_');
            gamma_ran(n) = trace(G_ran * G_ran');
            gamma_opt(n) = trace(G_opt * G_opt');
            gamma_d(n) = trace(G_d * G_d');
            pc_ran(:,k) = G_ran * (sym(:,k));
            pc_opt(:,k) = G_opt * (sym(:,k));
            pc_d(:,k) = G_d * (sym(:,k));
            n = n+1;
        end
        pc_ran = pc_ran./sqrt(gamma_ran);
        pc_opt = pc_opt./sqrt(gamma_opt);
        pc_d = pc_d./sqrt(gamma_d);
        ofdm_ran = ifft(pc_ran,fft_size,2)*sqrt(fft_size);
        ofdm_opt = ifft(pc_opt,fft_size,2)*sqrt(fft_size);
        ofdm_d = ifft(pc_d,fft_size,2)*sqrt(fft_size);
        cp_ran = [ofdm_ran(:,fft_size-cp_size+1:end) ofdm_ran];
        cp_opt = [ofdm_opt(:,fft_size-cp_size+1:end) ofdm_opt];
        cp_d = [ofdm_d(:,fft_size-cp_size+1:end) ofdm_d];
        %% pass channel
        for r = 1 : U
            for t = 1 : M
                receive_d(t,:) = conv(cp_d(t,:),h_BU(:,r,t).');
            end
            hx_d(r,:) = sum(receive_d,1);
        end
        
    end
    toc
end