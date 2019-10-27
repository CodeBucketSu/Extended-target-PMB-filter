function ggiw_mb_hat = extendedOA_onlyExist(ggiw_mbm,Wmbm,n_exist,model)
% variational approx
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

% Compute new optimal assignment
[newAssignment,costSum] = optimalAssignment(ggiw_mbm,ggiw_mb_hat,n_exist,num_mb,Wmbm);

maxIterations = 1e1;
numIterations = 1;
minCost = costSum;
while(numIterations < maxIterations)
    for i = 1:num_mb
        r_exist(i,:) = ggiw_mbm{i}.r(newAssignment(i,:));
        alpha_exist(i,:) = ggiw_mbm{i}.alpha(newAssignment(i,:));
        beta_exist(i,:) = ggiw_mbm{i}.beta(newAssignment(i,:));
        x_exist(:,i,:) = ggiw_mbm{i}.x(:,newAssignment(i,:));
        P_exist(:,:,i,:) = ggiw_mbm{i}.P(:,:,newAssignment(i,:));
        v_exist(i,:) = ggiw_mbm{i}.v(newAssignment(i,:));
        V_exist(:,:,i,:) = ggiw_mbm{i}.V(:,:,newAssignment(i,:));
    end
    % M-step
    ggiw_mb_hat = newMultiBernoulli(n_exist,Wmbm,r_exist,alpha_exist,...
        beta_exist,x_exist,P_exist,v_exist,V_exist,model);
    % E-step
    [newAssignment,costSumTemp] = optimalAssignment(ggiw_mbm,ggiw_mb_hat,n_exist,num_mb,Wmbm);
    if minCost - costSumTemp < 1e-3
        break;
    end
    if costSumTemp < minCost
        minCost = costSumTemp;
    end
    numIterations = numIterations + 1;
end

% Merge Bernoullis that correspond to the same new target, i.e.,
% KLD-based greedy MBM merging
ggiw_mb_new_hat = newTrackMerging(ggiw_mbm,Wmbm,n_exist,model);

ggiw_mb_hat = catenate(ggiw_mb_hat,ggiw_mb_new_hat);

end

function [newAssignment,costSum] = optimalAssignment(ggiw_mbm,ggiw_mb_hat,n,num_mb,Wmbm)
totalCost = zeros(num_mb,1);
cost = zeros(n,n);
newAssignment = zeros(num_mb,n);
for j = 1:num_mb
    for i = 1:n
        bernoulli2.r = ggiw_mb_hat.r(i);
        bernoulli2.alpha = ggiw_mb_hat.alpha(i);
        bernoulli2.beta = ggiw_mb_hat.beta(i);
        bernoulli2.x = ggiw_mb_hat.x(:,i);
        bernoulli2.P = ggiw_mb_hat.P(:,:,i);
        bernoulli2.v = ggiw_mb_hat.v(i);
        bernoulli2.V = ggiw_mb_hat.V(:,:,i);
        for ii = 1:n
            bernoulli1.r = ggiw_mbm{j}.r(ii);
            bernoulli1.alpha = ggiw_mbm{j}.alpha(ii);
            bernoulli1.beta = ggiw_mbm{j}.beta(ii);
            bernoulli1.x = ggiw_mbm{j}.x(:,ii);
            bernoulli1.P = ggiw_mbm{j}.P(:,:,ii);
            bernoulli1.v = ggiw_mbm{j}.v(ii);
            bernoulli1.V = ggiw_mbm{j}.V(:,:,ii);
            cost(i,ii) = BernoulliEntropy(bernoulli1,bernoulli2);
        end
    end
    x = min(cost(:));
    cost = cost - x;
    [newAssignment(j,:),totalCost(j)] = assignmentoptimal(cost);
    totalCost(j) = totalCost(j) + (x.*sum(newAssignment(j,:)>0,2))';
end
costSum = totalCost'*Wmbm;

end

function cost = BernoulliEntropy(bernoulli1,bernoulli2)
    % mb1: true
    % mb2: approx
    r = bernoulli1.r; alpha = bernoulli1.alpha; beta = bernoulli1.beta;
    x = bernoulli1.x; P = bernoulli1.P; v = bernoulli1.v; V = bernoulli1.V;
    r_hat = bernoulli2.r; alpha_hat = bernoulli2.alpha; beta_hat = bernoulli2.beta;
    x_hat = bernoulli2.x; P_hat = bernoulli2.P; v_hat = bernoulli2.v; V_hat = bernoulli2.V;
    r_hat(r_hat>1) = 1; r(r>1) = 1; % avoid numerical issue
    if r_hat == 1 && r == 1
        C1 = 0;
    else
        C1 = -(1-r)*log(1-r_hat) - r*log(r_hat);
    end
    C4 = crossEntropyIW(v,V,v_hat,V_hat);
    C2 = crossEntropyGaussian(x(1:2),P(1:2,1:2),x_hat(1:2),P_hat(1:2,1:2));
    C3 = crossEntropyGamma(alpha,beta,alpha_hat,beta_hat);
    cost = C1 + r*(C2+C3+C4);
end
