function ggiw_mb_new_hat = newTrackMerging_to(ggiw_mbm,Wmbm,n_exist)
%Funtion to merge new tracks

h_w = zeros(0,1);
h_r = zeros(0,1);
h_alpha = zeros(0,1);
h_beta = zeros(0,1);
h_x = zeros(4,0);
h_P = zeros(4,4,0);
h_v = zeros(0,1);
h_V = zeros(2,2,0);

num_mb = length(ggiw_mbm);
for i = 1:num_mb
    num_new_b = length(ggiw_mbm{i}.r) - n_exist;
    if num_new_b > 0
        h_w = [h_w;Wmbm(i)*ones(num_new_b,1)];
        h_r = [h_r;ggiw_mbm{i}.r(n_exist+1:end)];
        h_alpha = [h_alpha;ggiw_mbm{i}.alpha(n_exist+1:end)];
        h_beta = [h_beta;ggiw_mbm{i}.beta(n_exist+1:end)];
        h_x = [h_x ggiw_mbm{i}.x(:,n_exist+1:end)];
        h_P = cat(3,h_P,ggiw_mbm{i}.P(:,:,n_exist+1:end));
        h_v = [h_v;ggiw_mbm{i}.v(n_exist+1:end)];
        h_V = cat(3,h_V,ggiw_mbm{i}.V(:,:,n_exist+1:end));
    end
end

[~,IA,IC] = unique(h_x(1,:)'+h_alpha+h_v);
uniqueIC = unique(IC);
for i = 1:length(IA)
   h_w(IA(i)) =  sum(h_w(find(IC)==uniqueIC(i)));
end

ggiw_mb_new_hat.r = h_r(IA).*h_w(IA);
ggiw_mb_new_hat.alpha = h_alpha(IA);
ggiw_mb_new_hat.beta = h_beta(IA);
ggiw_mb_new_hat.x = h_x(:,IA);
ggiw_mb_new_hat.P = h_P(:,:,IA);
ggiw_mb_new_hat.v = h_v(IA);
ggiw_mb_new_hat.V = h_V(:,:,IA);
