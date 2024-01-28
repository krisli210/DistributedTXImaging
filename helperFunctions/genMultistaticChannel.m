function [H_tens, Z, ScatPosPol, UserPosPol, Psi, tau, alpha] = genMultistaticChannel(prm)
    % H_tens output as N x U x K
    % Z as N_R x N_Theta
    % ScatPosPol as 2 x L
    % UserPosPol as 2 x U x K
    % Psi as N x K x U x (N_R * N_Theta)
    % tau as U x K
    % alpha as U x K
    
    azInd = randperm(prm.N_theta, prm.L);
    rangeInd = randperm(prm.N_R, prm.L);

    azValues = prm.AzBins(azInd);
    rangeValues = prm.RangeBins(rangeInd);

    ScatPosPol = [rangeValues; azValues];
    ScatCoeff = ones(1, prm.L) .* complex(1, 1) ./ sqrt(2);

    Z = zeros(prm.N_R, prm.N_theta);

    H_tens = zeros(prm.N, prm.U, prm.K);
    
    % Channel Construction
    for l = 1:prm.L
        
        for k = 1:prm.K

            for u = 1:prm.U
                % For now assume that the user is randomly placed in grid

                theta_uk = prm.AzBins(randperm(prm.N_theta, 1));
                r_uk = prm.RangeBins(randperm(prm.N_R, 1));
                UserPosPol(:, u, k) = [r_uk, theta_uk];
                d_uk = sqrt(r_uk^2 + rangeValues(l)^2 - 2*r_uk*rangeValues(l)*cosd(theta_uk-azValues(l))) + rangeValues(l);
                tau_uk = d_uk / prm.PropagationSpeed;
                
                alpha_uk = (4*pi*d_uk / prm.lam)^-2;
                
                H_tens(:, u, k) = H_tens(:, u, k) + alpha_uk * ScatCoeff(l) ... 
                                  * exp(-1j * 2*pi * (k * prm.Delta_f*tau_uk)) ...
                                  .* exp(1j*2*pi*prm.DeltaRX * (0:prm.N-1).' * sind(azValues(l))); 
                
            end
        end

    end
