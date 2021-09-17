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
L1 = 3;     % BS-RIS paths
L2 = 3;     % RIS-UE paths
L3 = 3;     % BS-UE paths
scatter = 10;
snr = -10:5:30;
iter = 50;
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
BER_ran = zeros(1,length(snr));
BER_opt = zeros(1,length(snr));
%% RIS sim.
tic
for i = 1:iter
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
        
        for k = 1:fft_size
            H_BR_(:,:) = H_BR(k,:,:);
            H_RU_(:,:) = H_RU(k,:,:);
            H_BU_(:,:) = H_BU(k,:,:);
            %% eff. channel
            He_ran(k,:,:) = H_BU_ + H_RU_ * ang_ran * H_BR_;
            He_opt(k,:,:) = H_BU_ + H_RU_ * ang_opt * H_BR_;
            He_ran_(:,:) = He_ran(k,:,:);
            He_opt_(:,:) = He_opt(k,:,:);
            %% precoding (ZF)
            G_ran = He_ran_' * inv(He_ran_ * He_ran_');
            gamma_ran(n) = trace(G_ran * G_ran');
            pc_ran(:,k) = G_ran * (sym(:,k));
            G_opt = He_opt_' * inv(He_opt_ * He_opt_');
            gamma_opt(n) = trace(G_opt * G_opt');
            pc_opt(:,k) = G_opt * (sym(:,k));
            n = n+1;
        end
        pc_ran = pc_ran./sqrt(gamma_ran);
        ofdm_ran = ifft(pc_ran,fft_size,2)*sqrt(fft_size);
        cp_ran = [ofdm_ran(:,fft_size-cp_size+1:end) ofdm_ran];
        pc_opt = pc_opt./sqrt(gamma_opt);
        ofdm_opt = ifft(pc_opt,fft_size,2)*sqrt(fft_size);
        cp_opt = [ofdm_opt(:,fft_size-cp_size+1:end) ofdm_opt];
        %% pass channel
        he_ran = ifft(He_ran,fft_size,1);
        he_ran_ = he_ran(1:L1+L2-1,:,:);
        he_opt = ifft(He_opt,fft_size,1);
        he_opt_ = he_opt(1:L1+L2-1,:,:);
        for r = 1 : U
            for t = 1 : M
                receive_ran(t,:) = conv(cp_ran(t,:),he_ran_(:,r,t).');
                receive_opt(t,:) = conv(cp_opt(t,:),he_opt_(:,r,t).');
            end
            hx_ran(r,:) = sum(receive_ran,1);
            hx_opt(r,:) = sum(receive_opt,1);
        end
        %% Rx
        [y_ran, No_ran] = awgn_noise( hx_ran, snr(snr_i) );
        cp_remove_ran = y_ran(:,cp_size+1:fft_size+cp_size);
        y_hat_ran = cp_remove_ran.* sqrt(gamma_ran);
        ofdm_sym_rx_ran = fft(y_hat_ran,fft_size,2)/sqrt(fft_size);
        rx_data_ran = base_demod(ofdm_sym_rx_ran, mod_type);
        %BER_ran(snr_i) = BER_ran(snr_i) + sum( sum( data ~= rx_data_ran ) )/ (data_size * S);
        num_error_ran(i,snr_i) = biterr(data,rx_data_ran);
        
        [y_opt, No_opt] = awgn_noise( hx_opt, snr(snr_i) );
        cp_remove_opt = y_opt(:,cp_size+1:fft_size+cp_size);
        y_hat_opt = cp_remove_opt.* sqrt(gamma_opt);
        ofdm_sym_rx_opt = fft(y_hat_opt,fft_size,2)/sqrt(fft_size);
        rx_data_opt = base_demod(ofdm_sym_rx_opt, mod_type);
        %BER_opt(snr_i) = BER_opt(snr_i) + sum( sum( data ~= rx_data_opt ) )/ (data_size * S);
        num_error_opt(i,snr_i) = biterr(data,rx_data_opt);
    end
    
end
toc
%BER_ran = BER_ran / iter;
%BER_opt = BER_opt / iter;
BER_ran = (sum(num_error_ran,1)/(data_size*S))/iter;
BER_opt = (sum(num_error_opt,1)/(data_size*S))/iter;
%% figure
semilogy(snr, BER_ran, '-d');
hold on
semilogy(snr, BER_opt, '-d');
title('BER Performance')
legend('Ran','Opt')
ylabel('BER')
xlabel('SNR (dB)')
grid on