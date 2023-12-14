function [P_z0, P_z1] = metad_calc_sampler_verSOM(dprime, meta_dprime, cT1, cT2_0, cT2_1, N_samp)

    N_r2 = length(cT2_1)+1
    cT2_1
    for iz = 1:2
        z = (iz-1) .* ones(N_samp,1);

        [r1, r2] = second_ord_model_sampler(z,dprime,meta_dprime,cT1, cT2_0, cT2_1);

        N_r0 = flipud( accumarray(r2(r1==0),1,[N_r2 1]) );
        N_r1 = accumarray(r2(r1==1),1,[N_r2 1]);
        N_z = [N_r0' N_r1'];

        epsilon = 0.1;
        if iz==1
            P_z0 = (N_z + epsilon/N_r2) /(N_samp+epsilon);
        else
            P_z1 = (N_z + epsilon/N_r2) /(N_samp+epsilon);
        end

    end

