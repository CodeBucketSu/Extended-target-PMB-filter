function ggiw_mb_hat = extendedVMB_onlyExist(ggiw_mbm,Wmbm,n_exist,model)
% variational approxation
% ggiw_mb_hat: merged MB
% ggiw_mbm: MBM to be merged
% n: number of pre-existing tracks

num_mb = length(ggiw_mbm);

% remove bernoulli with zero existence probability
for i = 1:num_mb
    idx = ggiw_mbm{i}.r ~= 0;
    ggiw_mbm{i}.r = ggiw_mbm{i}.r(idx);
    ggiw_mbm{i}.alpha = ggiw_mbm{i}.alpha(idx);
    ggiw_mbm{i}.beta = ggiw_mbm{i}.beta(idx);
    ggiw_mbm{i}.x = ggiw_mbm{i}.x(:,idx);
    ggiw_mbm{i}.P = ggiw_mbm{i}.P(:,:,idx);
    ggiw_mbm{i}.v = ggiw_mbm{i}.v(idx);
    ggiw_mbm{i}.V = ggiw_mbm{i}.V(:,:,idx);
end

% Total num of Bernoullis for pre-existing tracks
len = num_mb*n_exist;

% Extract all the Bernoullis from MBM
h_r = zeros(len,1);
h_alpha = zeros(len,1);
h_beta = zeros(len,1);
h_x = zeros(4,len);
h_P = zeros(4,4,len);
h_v = zeros(len,1);
h_V = zeros(2,2,len);

for i = 1:num_mb
    indices = 1+(i-1)*n_exist:n_exist*i;
    h_r(indices) = ggiw_mbm{i}.r(1:n_exist);
    h_alpha(indices) = ggiw_mbm{i}.alpha(1:n_exist);
    h_beta(indices) = ggiw_mbm{i}.beta(1:n_exist);
    h_x(:,indices) = ggiw_mbm{i}.x(:,1:n_exist);
    h_P(:,:,indices) = ggiw_mbm{i}.P(:,:,1:n_exist);
    h_v(indices) = ggiw_mbm{i}.v(1:n_exist);
    h_V(:,:,indices) = ggiw_mbm{i}.V(:,:,1:n_exist);
end

[~,uniqueIdx] = unique(h_x(1,:)'+h_alpha+h_v);
h_r = h_r(uniqueIdx);
h_alpha = h_alpha(uniqueIdx);
h_beta = h_beta(uniqueIdx);
h_x = h_x(:,uniqueIdx);
h_P = h_P(:,:,uniqueIdx);
h_v = h_v(uniqueIdx);
h_V = h_V(:,:,uniqueIdx);

% Compute marginal data association probability
H = length(h_r);
phi_hat = zeros(H,n_exist);
for h = 1:H
    for j = 1:n_exist
        for i = 1:num_mb
            if isequal(h_x(1,h)+h_alpha(h)+h_v(h),...
                    ggiw_mbm{i}.x(1,j)+ggiw_mbm{i}.alpha(j)+ggiw_mbm{i}.v(j))
                phi_hat(h,j) = phi_hat(h,j) + Wmbm(i);
            end
        end
    end
end

ph = sum(phi_hat,2);

minCost = inf;
numIterations = 0;
maxIterations = 1e1;

% Variational approximation using Linear Programming
while(numIterations < maxIterations)
    numIterations = numIterations + 1;
    % M-step
    [r_temp,alpha_temp,beta_temp,x_temp,P_temp,v_temp,V_temp] = Mstep2...
        (h_r,h_alpha,h_beta,h_x,h_P,h_v,h_V,phi_hat,n_exist,model);
    % E-step
    [Cmin,phi_hat] = Estep2(r_temp,alpha_temp,beta_temp,x_temp,P_temp,v_temp,V_temp,...
        h_r,h_alpha,h_beta,h_x,h_P,h_v,h_V,ph,H,n_exist);
    
    if minCost - Cmin < 1e-3
        break;
    end
    if Cmin < minCost
        minCost = Cmin;
        ggiw_mb_hat.r = r_temp;ggiw_mb_hat.alpha = alpha_temp;ggiw_mb_hat.beta = beta_temp;
        ggiw_mb_hat.x = x_temp;ggiw_mb_hat.P = P_temp;ggiw_mb_hat.v = v_temp;ggiw_mb_hat.V = V_temp;
    end
end

% Merge Bernoullis that correspond to the same new target, i.e.,
% KLD-based greedy MBM merging
ggiw_mb_new_hat = newTrackMerging(ggiw_mbm,Wmbm,n_exist,model);

ggiw_mb_hat = catenate(ggiw_mb_hat,ggiw_mb_new_hat);

end


function [r_hat,alpha_hat,beta_hat,x_hat,P_hat,v_hat,V_hat] = Mstep2...
    (h_r,h_alpha,h_beta,h_x,h_P,h_v,h_V,phi,n,model)

r_hat = zeros(n,1);
alpha_hat = zeros(n,1);
beta_hat = zeros(n,1);
x_hat = zeros(4,n);
P_hat = zeros(4,4,n);
v_hat = zeros(n,1);
V_hat = zeros(2,2,n);

for i = 1:n
    w = phi(:,i);
    nonzero_idx = w~=0;
    if length(find(nonzero_idx==true)) == 1
        r_hat(i) = w(nonzero_idx).*h_r(nonzero_idx);
        alpha_hat(i) = h_alpha(nonzero_idx); beta_hat(i) = h_beta(nonzero_idx);
        x_hat(:,i) = h_x(:,nonzero_idx); P_hat(:,:,i) = h_P(:,:,nonzero_idx);
        v_hat(i) = h_v(nonzero_idx); V_hat(:,:,i) = h_V(:,:,nonzero_idx);
    else
        w_mixture = w(nonzero_idx).*h_r(nonzero_idx);
        alpha_mixture = h_alpha(nonzero_idx);
        beta_mixture = h_beta(nonzero_idx);
        x_mixture = h_x(:,nonzero_idx);
        P_mixture = h_P(:,:,nonzero_idx);
        v_mixture = h_v(nonzero_idx);
        V_mixture = h_V(:,:,nonzero_idx);
        % select the Bernoulli with highest (valid) existence probability
        [~,maxIdx] = max(w_mixture);
        I = maxIdx;
        % find other Bernoullis that has small enough KL divergence to the
        % selected Bernoulli
        for l = 1:length(w_mixture)
            if l~= maxIdx
                KLD = GGIW_KLdiff2(alpha_mixture(maxIdx),...
                    beta_mixture(maxIdx),x_mixture(:,maxIdx),P_mixture(:,:,maxIdx),...
                    v_mixture(maxIdx),V_mixture(maxIdx),alpha_mixture(l),...
                    beta_mixture(l),x_mixture(:,l),P_mixture(:,:,l),...
                    v_mixture(l),V_mixture(l));
                if KLD < model.GGIWMergingThreshold
                    I = [I l];
                end
            end
        end
        I = 1:length(w_mixture);
        r_hat(i) = sum(w_mixture);
        
        [~,alpha_hat(i),beta_hat(i),x_hat(:,i),P_hat(:,:,i),v_hat(i),V_hat(:,:,i)] = ...
            GGIW_merge(w_mixture(I),alpha_mixture(I),beta_mixture(I),...
            x_mixture(:,I),P_mixture(:,:,I),v_mixture(I),V_mixture(:,:,I));
    end
end

end

function [Cmin,phi_hat] = Estep2(r_hat,alpha_hat,beta_hat,x_hat,P_hat,v_hat,V_hat,...
    h_r,h_alpha,h_beta,h_x,h_P,h_v,h_V,ph,H,n)

% avoid numerical issues brought by adding floating numbers
r_hat(r_hat > 1) = 1;

C = zeros(H,n);
C1 = zeros(H,n);
C2 = zeros(H,n);
C3 = zeros(H,n);
C4 = zeros(H,n);

% Compute cross entropy bewteen two Bernoullis with GGIW density
for h = 1:H
    for i = 1:n
        if r_hat(i)==1 && h_r(h) == 1
            C1(h,i) = 0;
        else
            C1(h,i) = -(1-h_r(h))*log(1-r_hat(i)) - h_r(h)*log(r_hat(i));
        end

        C2(h,i) = crossEntropyGaussian(h_x(:,h),h_P(:,:,h),x_hat(:,i),P_hat(:,:,i));
        C3(h,i) = crossEntropyGamma(h_alpha(h),h_beta(h),alpha_hat(i),beta_hat(i));
        C4(h,i) = crossEntropyIW(h_v(h),h_V(:,:,h),v_hat(i),V_hat(:,:,i));

        C(h,i) = C1(h,i) + h_r(h)*(C2(h,i)+C3(h,i)+C4(h,i));
    end
end

% Solve the Linear Programming
[Cmin,phi_hat] = LP_transport(C,ph);

end
