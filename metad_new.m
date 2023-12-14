%%
clear all
close all

set(0,'DefaultFigureWindowStyle','docked')

dprime = 1.16;
means_t1 = [-dprime/2, dprime/2];
sigma_t1 = 1;
p = 0.5;
metadprime = [dprime*1.1, 1.3, 1.5];
Number_trials = 10^6;

% type 1 decision boundary
t1c = 0;

% t2r1 = [.05:.05:1];
% t2r1 = [0.1 0.3 0.6];
% type 2 decision boundaries when type 1 response was 1
t2r1 = [0.2 0.4];
% t2r1 = [0.3];
% type 2 decision boundaries when type 1 response was 0
t2r0 = fliplr(-t2r1);
% t2r0 = [-0.3];


%% matrix initialization
stimID   = zeros(length(metadprime), Number_trials);
response = zeros(length(metadprime), Number_trials);
rating   = zeros(length(metadprime), Number_trials);

for iM = 1:length(metadprime)
    metadprime_now = metadprime(iM);
    sde = sqrt((dprime/metadprime_now)^2-1);

    % simulate stimuli
    S = binornd(1, p, 1, Number_trials);
    stimID(iM, :) = S;

    % simulate responses
    [r_t1,r_t2]=second_ord_model_sampler(S,dprime,metadprime_now,t1c,t2r0,t2r1);

    response(iM, :) = r_t1;
    rating(iM, :) = r_t2;
end

nRating = length(t2r0)+1;
md_idx = 1;
[nR_S1_0, nR_S2_0] = trials2counts(stimID(md_idx,:), response(md_idx,:), rating(md_idx,:), nRating, 1);


% sim = metad_sim_Fleming(1.5, 1.3, 0, -.3, .3, Number_trials);
% nR_S1_0 = sim.nR_S1;
% nR_S2_0 = sim.nR_S2;

%% fitting the model

options.UncertaintyHandling = 1;    % Tell BADS that the objective is noisy
options.MaxIter = 200;%150;
options.MaxFunEvals = 200;

l_v = [0.1  0.1  -2   0    0];
u_v = [3    3     2   1.5  1.5];
for rep = 1:1 

    init_v = l_v;
    for i = 1:length(l_v)
        init_v(i) = rand*(u_v(i)-l_v(i))+l_v(i);
    end

    [par_,fval,exitflag,output] = bads(@(x) metad_calc_nlLH(x,nR_S1_0,nR_S2_0),init_v,l_v,u_v,[],[],[],options);

    par_m(rep,:) = par_
    fval_v(rep)  = fval
end
%%
disp('ours: '+string((par_m(2)/metadprime(1)-1)*100))


res = fit_meta_d_MLE(nR_S1_0, nR_S2_0);
% 
% (res.meta_da/metadprime(1)-1)*100
disp('theirs: '+string((res.meta_da/metadprime(1)-1)*100))
return
% 
% grid_x = linspace(-10,10,100);
% 
% [P_rT2e1_G_rT1e0_z, P_rT2e1_G_rT1e1_z] = calc_probs(grid_x',1.5,0.99,0,-.3,.3);
% P_R = [P_rT2e1_G_rT1e0_z; 1 - P_rT2e1_G_rT1e0_z; 1 - P_rT2e1_G_rT1e1_z; P_rT2e1_G_rT1e1_z]
% 
% [p_S0, p_S1, xT1_G_z, xT2_G_z] = new_metad_calc_sampler(1.5,0.99,0,.3, 10^5);
% 
% % histogram(xT1_G_z)
% % hold on
% close
% hv1(1,:) =histcounts(xT2_G_z(xT1_G_z(:,1)<0,1),grid_x);
% hv2(1,:) =histcounts(xT2_G_z(xT1_G_z(:,2)<0,2),grid_x);
% 
% hv1(2,:) = histcounts(xT2_G_z(xT1_G_z(:,1)>0,1),grid_x);
% hv2(2,:) = histcounts(xT2_G_z(xT1_G_z(:,2)>0,2),grid_x);
% 
% hv1 = hv1./sum(hv1,'all');
% hv2 = hv2./sum(hv2,'all');
% 
% plot(hv1')
% hold on
% plot(hv2')

