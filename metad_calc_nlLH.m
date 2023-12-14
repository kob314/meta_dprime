function lLH = metad_calc_nlLH(inp,nR_S0,nR_S1)
    dprime      = inp(1);
    meta_dprime = inp(2);
    cT1         = inp(3);
    a           = inp(4);
    b           = inp(5);

    cT2_0 = [cT1-a-b, cT1-a];
    cT2_1 = [cT1+a, cT1+a+b];
    N_samp = 10^5;

    % [p_S0, p_S1] = metad_calc_sampler_verMetanoise(dprime, meta_dprime, cT1, cT2, N_samp);
    [p_S0, p_S1] = metad_calc_sampler_verSOM(dprime, meta_dprime, cT1, cT2_0, cT2_1, N_samp);
    % p_S0
    % nR_S0

    size(p_S0)
    size(nR_S0)

    lLH = sum(log(p_S0) .* nR_S0 + log(p_S1) .* nR_S1);

    % rtx
    if abs(lLH) > 10^10
        lLH = 10^10;
    end

    lLH = -lLH;
end