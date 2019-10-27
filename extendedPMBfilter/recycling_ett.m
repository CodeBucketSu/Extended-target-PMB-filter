function [ ggiw_mb, ggiw_ppp_update ] = recycling_ett( ggiw_mb, ggiw_ppp_update,model )

% Recycling: approximate Bernoulli components with very small existence
% probability as PPP

idx1 = ggiw_mb.r >= model.threshold_recycle;
idx2 = ggiw_mb.r < model.threshold_recycle & ggiw_mb.r > model.threshold_r;

ggiw_ppp_update.wu = [ggiw_ppp_update.wu;ggiw_mb.r(idx2)];
ggiw_ppp_update.alpha_u = [ggiw_ppp_update.alpha_u;ggiw_mb.alpha(idx2)];
ggiw_ppp_update.beta_u = [ggiw_ppp_update.beta_u;ggiw_mb.beta(idx2)];
ggiw_ppp_update.xu = [ggiw_ppp_update.xu ggiw_mb.x(:,idx2)];
ggiw_ppp_update.Pu = cat(3,ggiw_ppp_update.Pu,ggiw_mb.P(:,:,idx2));
ggiw_ppp_update.vu = [ggiw_ppp_update.vu;ggiw_mb.v(idx2)];
ggiw_ppp_update.Vu = cat(3,ggiw_ppp_update.Vu,ggiw_mb.V(:,:,idx2));

ggiw_mb.r = ggiw_mb.r(idx1);
ggiw_mb.alpha = ggiw_mb.alpha(idx1);
ggiw_mb.beta = ggiw_mb.beta(idx1);
ggiw_mb.x = ggiw_mb.x(:,idx1);
ggiw_mb.P = ggiw_mb.P(:,:,idx1);
ggiw_mb.v = ggiw_mb.v(idx1);
ggiw_mb.V = ggiw_mb.V(:,:,idx1);

end

