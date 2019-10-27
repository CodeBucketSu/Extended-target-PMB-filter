function [ggiw_mb,ggiw_ppp_update,est] = updating(ggiw_mb,ggiw_ppp,Z,model)

% Gating and clustering for pre-existing targets
[gatingGroup,gatingGGIW,idx_out] = groupGating(Z,ggiw_mb,model);

% Gating for unknown targets; find the measurements that are outside the
% gates of pre-existing targets but are inside the gates of unknown targets
Z_gate = gate_meas_gms(Z(:,idx_out),ggiw_ppp,model);

% Generate measurement partitions using db-scan
Partitions = genPartitions_dbscan(Z_gate,model);

% Find the most likely measurement partition for unknown targets
bestPartition = findBestPartition(ggiw_ppp,Partitions,Z_gate,model);

% Create the new GGIW components for unknown targets
ggiw_new = newBernoulli(ggiw_ppp,bestPartition,Z_gate,model);

% Update PPP for unknown targets
ggiw_ppp_update = ppp_update(ggiw_ppp,model);

G = length(gatingGroup);
if G==0
    % Pruning
    ggiw_mb = trackPruning(ggiw_new,model);
    % Extract target states
    est = state_extract(ggiw_new);
    return;
end

ggiw_mbm_gate = cell(G,1);
Wg = cell(G,1);
ggiw_vb_upd = ggiw_new;

% Perform GGIW update independently for each clustering group
for g = 1:G
    
    % Generate measurement partitions using db-scan
    MPs = genPartitions_dbscan(Z(:,gatingGroup{g}),model);
    if isempty(MPs) % this corresponds to a cluster w/o measurement
        nP = 1;
        MPs = cell(1);
    else
        nP = length(MPs);
    end
    
    log_Wp = zeros(0,1);
    ggiw_mbm_gate_hat = cell(0,1);
    
    % Perform GGIW update for each measurement partition
    for i = 1:nP
        % Create Bernoulli components for PPP
        [ggiw_mb_new,Lnew] = newBernoulli(ggiw_ppp,MPs{i},Z(:,gatingGroup{g}),model);
        
        % Create Bernoulli components for MB
        [ggiw_mb_miss,ggiw_mb_upd,Lmiss,Lupd] = mbm_update(gatingGGIW{g},model,Z(:,gatingGroup{g}),MPs{i});
        
        % Data association
        [ggiw_mbm_murty,log_Wmurty] = dataAssociation(ggiw_mb_new,Lnew,...
            ggiw_mb_upd,Lupd,ggiw_mb_miss,Lmiss,model);

        % Concatenate MBM obtained for each MB
        log_Wp = [log_Wp;log_Wmurty];
        ggiw_mbm_gate_hat = cat(1,ggiw_mbm_gate_hat,ggiw_mbm_murty);
    end

    % Pruning by only keeping MBs with total weights that correspond to
    % 1-model.threshold_w 
    Wp_normalized = exp(normalizeLogWeights(log_Wp));
    [Wp_sorted,order] = sort(Wp_normalized,'descend');
    pos = find(cumsum(Wp_sorted) >= 1-model.threshold_w,1);
    Wg{g} = exp(normalizeLogWeights(log_Wp(order(1:pos))));
    ggiw_mbm_gate{g} = ggiw_mbm_gate_hat(order(1:pos));
    
    % MBM merging
    N = length(gatingGGIW{g}.r);
    switch(model.VMBtype)
        case 1
            ggiw_mb_hat = extendedTO_onlyExist(ggiw_mbm_gate{g},Wg{g},N,model);
        case 2
            ggiw_mb_hat = extendedVMB_onlyExist(ggiw_mbm_gate{g},Wg{g},N,model);
        case 3
            ggiw_mb_hat = extendedOA_onlyExist(ggiw_mbm_gate{g},Wg{g},N,model);
        case 4
            ggiw_mb_hat = extendedTO_onlyExist_noNewMerge(ggiw_mbm_gate{g},Wg{g},N,model);
    end
    % Concatenate merged MB for each cluster
    ggiw_vb_upd = catenate(ggiw_mb_hat,ggiw_vb_upd);
end

% Pruning unlikely tracks
ggiw_vb_upd = trackPruning(ggiw_vb_upd,model);

% Recycling
[ggiw_mb,ggiw_ppp_update] = recycling_ett(ggiw_vb_upd,ggiw_ppp_update,model);

% Target states extraction
est = state_extract(ggiw_mb);

end

