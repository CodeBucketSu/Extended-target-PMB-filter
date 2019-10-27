function ggiw_mb_new_hat = newTrackMerging(ggiw_mbm,Wmbm,n_exist,model)
%Funtion to merge new tracks using KLD-based greedy approach

num_mb = length(ggiw_mbm);
new_tracks = cell(num_mb,1);
for i = 1:num_mb
    if length(ggiw_mbm{i}.r) == n_exist
        new_tracks{i}.r = zeros(0,1);
        new_tracks{i}.alpha = zeros(0,1);
        new_tracks{i}.beta = zeros(0,1);
        new_tracks{i}.x = zeros(4,0);
        new_tracks{i}.P = zeros(4,4,0);
        new_tracks{i}.v = zeros(0,1);
        new_tracks{i}.V = zeros(2,2,0);
        new_tracks{i}.l = zeros(0,1);
    else
        new_tracks{i}.r = ggiw_mbm{i}.r(n_exist+1:end);
        new_tracks{i}.alpha = ggiw_mbm{i}.alpha(n_exist+1:end);
        new_tracks{i}.beta = ggiw_mbm{i}.beta(n_exist+1:end);
        new_tracks{i}.x = ggiw_mbm{i}.x(:,n_exist+1:end);
        new_tracks{i}.P = ggiw_mbm{i}.P(:,:,n_exist+1:end);
        new_tracks{i}.v = ggiw_mbm{i}.v(n_exist+1:end);
        new_tracks{i}.V = ggiw_mbm{i}.V(:,:,n_exist+1:end);
        new_tracks{i}.l = zeros(length(new_tracks{i}.r),1);
    end
end

% get labels
label = 0;
for i = 1:num_mb
    for ii = 1:length(new_tracks{i}.r)
        if new_tracks{i}.l(ii) == 0
            label = label + 1;
            new_tracks{i}.l(ii) = label;
        end
        for j = i+1:num_mb
            nj = length(new_tracks{j}.r);
            temp2 = inf*ones(nj,1);
            for jj = 1:nj
                if new_tracks{j}.l(jj) == 0
                    temp2(jj,1) = GGIW_KLdiff2(new_tracks{i}.alpha(ii),...
                        new_tracks{i}.beta(ii),new_tracks{i}.x(:,ii),new_tracks{i}.P(:,:,ii),...
                        new_tracks{i}.v(ii),new_tracks{i}.V(:,:,ii),...
                        new_tracks{j}.alpha(jj),new_tracks{j}.beta(jj),new_tracks{j}.x(:,jj),...
                        new_tracks{j}.P(:,:,jj),new_tracks{j}.v(jj),new_tracks{j}.V(:,:,jj));
                end
            end
            [~,idx] = min(temp2);
            if temp2(idx) < model.newTrackMergingThreshold
                new_tracks{j}.l(idx) = label;
            end
        end
    end
end

% reconstruct new_tracks
r_new = zeros(num_mb,label);
alpha_new = zeros(num_mb,label);
beta_new = zeros(num_mb,label);
x_new = zeros(4,num_mb,label);
P_new = zeros(4,4,num_mb,label);
v_new = zeros(num_mb,label);
V_new = zeros(2,2,num_mb,label);

for i = 1:num_mb
    r_new(i,new_tracks{i}.l) = new_tracks{i}.r;
    alpha_new(i,new_tracks{i}.l) = new_tracks{i}.alpha;
    beta_new(i,new_tracks{i}.l) = new_tracks{i}.beta;
    x_new(:,i,new_tracks{i}.l) = new_tracks{i}.x;
    P_new(:,:,i,new_tracks{i}.l) = new_tracks{i}.P;
    v_new(i,new_tracks{i}.l) = new_tracks{i}.v;
    V_new(:,:,i,new_tracks{i}.l) = new_tracks{i}.V;
end

ggiw_mb_new_hat.r = zeros(label,1);
ggiw_mb_new_hat.alpha = zeros(label,1);
ggiw_mb_new_hat.beta = zeros(label,1);
ggiw_mb_new_hat.x = zeros(4,label);
ggiw_mb_new_hat.P = zeros(4,4,label);
ggiw_mb_new_hat.v = zeros(label,1);
ggiw_mb_new_hat.V = zeros(2,2,label);

% Merge Bernoullis with the same label
for i = 1:label
    nonzero_idx = r_new(:,i)~=0;
    w = Wmbm(nonzero_idx).*r_new(nonzero_idx,i);
    ggiw_mb_new_hat.r(i) = sum(w);
    [~,ggiw_mb_new_hat.alpha(i),ggiw_mb_new_hat.beta(i),ggiw_mb_new_hat.x(:,i),ggiw_mb_new_hat.P(:,:,i)...
        ,ggiw_mb_new_hat.v(i),ggiw_mb_new_hat.V(:,:,i)] = ...
            GGIW_merge(w,alpha_new(nonzero_idx,i),beta_new(nonzero_idx,i),...
            x_new(:,nonzero_idx,i),P_new(:,:,nonzero_idx,i),v_new(nonzero_idx,i),V_new(:,:,nonzero_idx,i));
end

end

