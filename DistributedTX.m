close all
clear 

rng(52);

% % % OFDM Signal Params

    prm.CenterFreq = 28e9;
    prm.PropagationSpeed = physconst('LightSpeed');
    prm.lam = prm.PropagationSpeed/prm.CenterFreq;

    prm.Delta_f = 120*1e3; % SCS in KHz
    
    prm.NRB = 60; % number of resource blocks
    prm.K = 12*prm.NRB;
    
    prm.U = 1; % U per RB
    prm.N_T = 1; % number of time slots
    prm.Nofdm = 14; %number of OFDM symbols per slot
    prm.MCS = 16; %modulation order

    
% % % Receive Array Definitions

    prm.RxPos = [0; 0; 0]; % Polar
    prm.N = 16;
    prm.DeltaRX = .5; % Half wavelength for now, ideally M*.5 but might introduce ambiguity 
    prm.NumRxElements = prod(prm.RxArraySize);
    prm.RxAZlim = prm.BsAZlim;
    prm.RxELlim = [-90 0];

% % % Scene Construction
    
    prm.N_theta = 32;
    thetaMin = prm.RxAZlim(1); thetaMax = prm.RxAZlim(2); %in Azimuth
    prm.AzBins = thetaMin:(thetaMax-thetaMin)/(prm.N_theta-1):thetaMax;

    prm.L = 10;
    prm.rMin = 20; prm.rMax = 70;

    limit_rangeRes = true;
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

    [H_tens, Z, ScatPosPol] = genMultistaticChannel(prm);

    [txGrid] = genUplinkFreqTxGrid(prm.NumUsers, prm.MCS, prm.N_T, prm.Nofdm, prm.K); % (Nofdm * N_T) x M x K
    



