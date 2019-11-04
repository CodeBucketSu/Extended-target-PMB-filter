clc;clear
dbstop if error

% VMBtype
% Type 1: Track-oriented merging w/ greedy new tracking merging
% Type 2: Efficient approximation using Linear Assignment
% Type 3: Optimal assignment
% Type 4: Track-oriented merging w/o greedy new tracking merging

load('ScenarioWith27targets.mat')

para.VMBtype = 2;           % choose an MB merging algorithm to use
para.false_alarm_rate = 60; % choose false alarm rate
para.detection_prob = 0.9;  % choose detection probability

% Model paramters and ground truth
[groundTruth,Z,model] = gen_data_many_targets(para,targetTracks);

% Total time step
K = length(groundTruth); 

% (G)OSPA parameters
gospa_p = 1;
gospa_c = 10;
gospa_vals = zeros(K,4);
ospa_vals = zeros(K,3);

estimationResults = cell(K,1);

% Initialise multi-Bernoulli and PPP parameters
[ggiw_mb,ggiw_ppp] = Initialisation(model);

fprintf('Current time step: ')
for t = 1:K
    
    fprintf('%g ',t);
    
    % Prediction step
    [ggiw_mb,ggiw_ppp] = predicting(ggiw_mb,ggiw_ppp,model);
    
    % Update step
    [ggiw_mb,ggiw_ppp,est] = updating(ggiw_mb,ggiw_ppp,Z{t},model);
    
    % Performance evaluation using GOSPA
    gospa_vals(t,:) = GOSPAmetric(est,groundTruth{t},gospa_c,gospa_p);
    % Performance evaluation using OSPA
    ospa_vals(t,:) = OSPAmetric(est,groundTruth{t},gospa_c,gospa_p);
    
    % Store estimates
    estimationResults{t} = est;
end

% Compute (G)OSPA error, averaged per MB run and per time step
averGospa = mean(mean(gospa_vals,3));
averOspa = mean(mean(ospa_vals,3));

fprintf('\nGOSPA error: %g, Localisation error: %g, Missed detection: %g, False detection: %g', ...
    averGospa(1), averGospa(2), averGospa(3), averGospa(4));

fprintf('\nOSPA error: %g, Localisation error: %g, Cardinality error: %g\n', ...
    averOspa(1), averOspa(2), averOspa(3));
