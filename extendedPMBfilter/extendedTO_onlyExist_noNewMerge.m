function ggiw_mb_hat = extendedTO_onlyExist_noNewMerge(ggiw_mbm,Wmbm,n_exist,model)
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

% reconstruct MBM parameters
r_exist = zeros(num_mb,n_exist);
alpha_exist = zeros(num_mb,n_exist);
beta_exist = zeros(num_mb,n_exist);
x_exist = zeros(4,num_mb,n_exist);
P_exist = zeros(4,4,num_mb,n_exist);
v_exist = zeros(num_mb,n_exist);
V_exist = zeros(2,2,num_mb,n_exist);

for i = 1:num_mb
    r_exist(i,:) = ggiw_mbm{i}.r(1:n_exist);
    alpha_exist(i,:) = ggiw_mbm{i}.alpha(1:n_exist);
    beta_exist(i,:) = ggiw_mbm{i}.beta(1:n_exist);
    x_exist(:,i,:) = ggiw_mbm{i}.x(:,1:n_exist);
    P_exist(:,:,i,:) = ggiw_mbm{i}.P(:,:,1:n_exist);
    v_exist(i,:) = ggiw_mbm{i}.v(1:n_exist);
    V_exist(:,:,i,:) = ggiw_mbm{i}.V(:,:,1:n_exist);
end

% Merge Bernoullis that correspond to the same pre-target, i.e.,
% track-oriented MBM merging
ggiw_mb_hat = newMultiBernoulli(n_exist,Wmbm,r_exist,alpha_exist,...
    beta_exist,x_exist,P_exist,v_exist,V_exist,model);

% Merge Bernoullis that correspond to the same new target, i.e.,
% track-oriented MBM merging
ggiw_mb_new_hat = newTrackMerging_to(ggiw_mbm,Wmbm,n_exist);

ggiw_mb_hat = catenate(ggiw_mb_hat,ggiw_mb_new_hat);

end