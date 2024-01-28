function [H_tens, Z, ScatPosPol] = genMultistaticChannel(prm)
    % H_tens output as N x U x K
    % Z as N_R x N_Theta
    % ScatPosPol as 3 x L

    azInd = randperm(prm.N_theta, prm.L);
    rangeInd = randperm(prm.N_R, prm.L):

    azValues = prm.AzBins(azInd);
    rangeValues = prm.RangeBins(rangeInd);

    ScatPosPol = [rangeValues; azValues; zeros(1, prm.L)];
    ScatCoeff = ones(1, prm.L) .* complex(1, 1) ./ sqrt(2);

    Z = zeros(prm.N_R, prm.N_theta);

    H_tens = zeros(prm.N, prm.U, prm.K);
    
    for l = 1:prm.L
        
        for u = 1:prm.U

        end

    end

end