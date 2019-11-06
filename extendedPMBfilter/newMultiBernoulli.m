function ggiw_mb_hat = newMultiBernoulli(n_exist,Wmbm,r_exist,alpha_exist,...
    beta_exist,x_exist,P_exist,v_exist,V_exist,model)

ggiw_mb_hat.r = zeros(n_exist,1);
ggiw_mb_hat.alpha = zeros(n_exist,1);
ggiw_mb_hat.beta = zeros(n_exist,1);
ggiw_mb_hat.x = zeros(4,n_exist);
ggiw_mb_hat.P = zeros(4,4,n_exist);
ggiw_mb_hat.v = zeros(n_exist,1);
ggiw_mb_hat.V = zeros(2,2,n_exist);

for i = 1:n_exist
    w = Wmbm.*r_exist(:,i);
    % select the Bernoulli with highest (valid) existence probability
    [~,maxIdx] = max(w);    
    I = maxIdx;
    % find other Bernoullis that has small enough KL divergence to the
    % selected Bernoulli
    for l = 1:length(w)
        if l~= maxIdx
            KLD = GGIW_KLdiff2(alpha_exist(maxIdx,i),beta_exist(maxIdx,i),x_exist(:,maxIdx,i),...
                P_exist(:,:,maxIdx,i),v_exist(maxIdx,i),V_exist(:,:,maxIdx,i),alpha_exist(l,i),...
                beta_exist(l,i),x_exist(:,l,i),P_exist(:,:,l,i),v_exist(l,i),V_exist(:,:,l,i));
            if KLD < model.GGIWMergingThreshold
                I = [I l];
            end
        end
    end

    ggiw_mb_hat.r(i) = sum(w);
    
    % GGIW merging
    [~,ggiw_mb_hat.alpha(i),ggiw_mb_hat.beta(i),ggiw_mb_hat.x(:,i),ggiw_mb_hat.P(:,:,i)...
        ,ggiw_mb_hat.v(i),ggiw_mb_hat.V(:,:,i)] = GGIW_merge(w(I),alpha_exist(I,i),...
        beta_exist(I,i),x_exist(:,I,i),P_exist(:,:,I,i),v_exist(I,i),V_exist(:,:,I,i));
end

end
