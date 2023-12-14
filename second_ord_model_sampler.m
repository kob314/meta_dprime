% [a b]=second_ord_model_sampler1(1+[1 1 1 0 0]',1.5,2,0,.3)

function [r1, r2] = second_ord_model_sampler(z,dprime,meta_dprime,cT1, cT2_0, cT2_1)

    sigma_x1  = 1;                  % type 1 std
    sigma_x2  = dprime/meta_dprime; % type 2 std

    rho = min(1/sigma_x2,sigma_x2);

    cov_m = [1, rho*sigma_x2; rho*sigma_x2, sigma_x2.^2]; % bivariate normal covariance matrix

    out = mvnrnd([1,1]*dprime.*(z(:)-.5),cov_m);
    x1  = out(:,1); % type 1 observation
    x2  = out(:,2); % type 2 observation

    r1  = (x1>cT1);  % type 1 decision
    r2  = nan(size(r1));

    r2(r1==0)  = discretize(-x2(r1==0),[-inf flipud(-cT2_0(:))' inf]);
    r2(r1==1)  = discretize(x2(r1==1),[-inf cT2_1(:)' inf]);


end


