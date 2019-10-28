function [ ggiw_ppp ] = ppp_update( ggiw_ppp,model )
% GGIW PPP UPDATE: MISSED DETECTION

Pd = model.Pd;
% Update undetected PPP
temp1 = (ggiw_ppp.beta_u./(ggiw_ppp.beta_u+1)).^ggiw_ppp.alpha_u;

wu1 = (1-Pd)*ggiw_ppp.wu;
wu2 = Pd*temp1.*ggiw_ppp.wu;

% Gamma mixture reduction
beta = ggiw_ppp.beta_u;
n = length(ggiw_ppp.wu);
for i = 1:n
    w_bar = wu1(i) + wu2(i);
    temp = 1/w_bar*ggiw_ppp.alpha_u(i)*(wu1(i)/beta(i)+wu2(i)/(beta(i)+1));
    beta_hat = ggiw_ppp.alpha_u(i)/temp;
    ggiw_ppp.beta_u(i) = beta_hat;
end

ggiw_ppp.wu = wu1 + wu2;

% Prune undetected PPP intensity
idx1 = ggiw_ppp.wu > model.threshold_u;
idx2 = -log((ggiw_ppp.beta_u./(ggiw_ppp.beta_u+1)).^ggiw_ppp.alpha_u)>1;
idx3 = ggiw_ppp.vu > 6;

idx = idx1&idx2&idx3;

ggiw_ppp.wu = ggiw_ppp.wu(idx);
ggiw_ppp.alpha_u = ggiw_ppp.alpha_u(idx);
ggiw_ppp.beta_u = ggiw_ppp.beta_u(idx);
ggiw_ppp.xu = ggiw_ppp.xu(:,idx);
ggiw_ppp.Pu = ggiw_ppp.Pu(:,:,idx);
ggiw_ppp.vu = ggiw_ppp.vu(idx);
ggiw_ppp.Vu = ggiw_ppp.Vu(:,:,idx);

end
