close all
clear 

rng(2);

% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 20; % number of resource blocks
    prm.K = 12*prm.NRB;
    
    prm.U = 8; % U per RB
    prm.N_T = 1; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.N_s = prm.N_T*prm.Nofdm;
    prm.MCS = 16; %modulation order

% % % Receive Array Definitions

    prm.RxPos = [0; 0; 0]; % Polar
    prm.N = 16;
    prm.DeltaRX = .5; % Half wavelength for now, ideally M*.5 but might introduce ambiguity 
    prm.NumRxElements = prod(prm.N);
    prm.RxAZlim = [-60 60];
    prm.RxELlim = [-90 0];

% % % Scene Construction
    
    prm.N_theta = 128;
    thetaMin = prm.RxAZlim(1); thetaMax = prm.RxAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;

    prm.L = 150;
    prm.rMin = 20; prm.rMax = 70;

    limit_rangeRes = false;
    if (limit_rangeRes)
        prm.delta_R = prm.PropagationSpeed./(2*prm.Delta_f*prm.K); % nominal range resolution
    else
        % Block for changing rangeRes to a non-nominal value
        prm.delta_R = 2;
    end
    prm.WholeRange = 0:prm.delta_R:(prm.K-1)*prm.delta_R;
    minIndex = find(prm.WholeRange < prm.rMin, 1, 'last')+1;
    maxIndex = find(prm.WholeRange > prm.rMax, 1, 'first')-1;
    
    prm.RangeBins = prm.WholeRange(minIndex:maxIndex);
    prm.N_R = numel(prm.RangeBins);

    [H_tens, Z, ScatPosPol, UserPosPol, Psi] = genMultistaticChannel(prm);
    
    prm.Pt_dBm = 20;
    prm.Pt_W = 10^(prm.Pt_dBm/10)*1e-3; % Watts

    [txGrid] = genUplinkFreqTxGrid(prm.U, prm.MCS, prm.N_T, prm.Nofdm, prm.K, prm.Pt_W); % (Nofdm * N_T) x K
    
    
    Y_tens = zeros(prm.N, prm.Nofdm * prm.N_T, prm.K);
    
    A_Theta = zeros(prm.N * prm.N_s * prm.K, prm.N_R*prm.N_theta);

    for k = 1:prm.K
        X_k = txGrid(:, :, k);
        Y_k = zeros(prm.N, prm.N_s);
        for n_s = 1:prm.N_s
            Y_k(:, n_s) = H_tens(:, :, k) * X_k(:, n_s); % this has multiuser multiplexed into N RX antennas
        end
        Y_tens(:, :, k) = Y_k;

        Y_kron_k = reshape(Y_k, [prm.N * prm.N_s, 1]);
        A_Theta_k = zeros(prm.N*prm.N_s, prm.N_R*prm.N_theta);
        for u = 1:prm.U
            A_Theta_k = A_Theta_k + kron(X_k(u, :).', eye(prm.N)) * squeeze(Psi(:, u, k, :));
            % A_Theta_uk = kron(X_k(u, :).', eye(prm.N)) * squeeze(Psi(:, u, k, :)); % CHange for multiuser
        end
        % A_Theta_k = kron(X_k.', eye(prm.N)) * squeeze(Psi(:, 1, k, :)); % CHange for multiuser

        start_k_index = (k-1) * prm.N * prm.N_s; % no idea if this works
        A_Theta(start_k_index + 1:start_k_index + prm.N*prm.N_s, :) = A_Theta_k;
    end
    Z_hat_vec = omp(A_Theta, Y_tens(:), prm.L, 1e-20);
    Z_hat = reshape(Z_hat_vec, [prm.N_R prm.N_theta]);

    peaksnr = psnr(abs(Z_hat).^2, abs(Z).^2, max(abs(Z).^2, [], "all"));

    figure;
    imagesc(prm.AzBins, prm.RangeBins, abs(Z));
    xlabel('Azimuth [$^\circ$]', 'Interpreter','latex', 'FontSize',14);
    ylabel('Range $[m]$', 'Interpreter','latex', 'FontSize',14)
    title('True','FontSize',14, 'Interpreter','latex');
    % c=colorbar;
    % c.Label.String = '$|\textbf Z|$';
    % c.Label.Interpreter = 'Latex';
    % c.Label.Rotation = 360;
    % c.Label.FontSize = 18;

    figure;
    imagesc(prm.AzBins, prm.RangeBins, abs(Z_hat));
    xlabel('Azimuth [$^\circ$]', 'Interpreter','latex', 'FontSize',14);
    ylabel('Range $[m]$', 'Interpreter','latex', 'FontSize',14)
    title({ ...
        ['Estimate'], ...
        ['$PSNR$ = ', num2str(peaksnr, 4), ' dB '] ...
        }, 'Interpreter','latex','FontSize',14)


