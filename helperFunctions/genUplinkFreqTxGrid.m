function [txGrid] = genUplinkFreqTxGrid(U, MCS, N_T, Nofdm, K, Pt_W)
    %Generates baseband equivalent frequency-domain signaling
    % txGrid output is (Nofdm * N_T) x M x K
    % txGrid output is (M x Nofdm * N_T x K)
    txGrid = zeros(U, Nofdm * N_T, K);
    
    sIndices = randi([0 MCS-1], [U, Nofdm * N_T, K]);     % per user symbols given as  Nofdm * N_T x U x K
    s = qammod(sIndices, MCS, 'UnitAveragePower', true); % this should take care of E[ss^H] = I_U / U ?
      
    % precoded transmit symbols given as Nofdm * N_T x M x K

    % Loop over subcarriers and slots because idk how to tensorize this
    % txCodebook = txCodebook(:, 1:200);
    % txCodebook = txCodebook(:, 128);
        for n_T = 1:N_T
            for nofdm = 1:Nofdm
                for k = 1:K
                    startTimeIndex = (Nofdm * (n_T - 1));
                    s_slice = squeeze(s(:, startTimeIndex + nofdm, k)); % U x 1
                    unnormed_symbol_vector = s_slice;
                    a = Pt_W * norm(unnormed_symbol_vector*unnormed_symbol_vector', 'fro');
                    txGrid(:, startTimeIndex + nofdm, k) = a*unnormed_symbol_vector; 
                end
            end
        end
end