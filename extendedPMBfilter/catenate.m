function ggiw_mb_catenate = catenate(ggiw_mb1,ggiw_mb2)
% catenate ggiw_mb1 with ggiw_mb2
if isempty(ggiw_mb2) || isempty(ggiw_mb2.r)
    ggiw_mb_catenate = ggiw_mb1;
elseif isempty(ggiw_mb1) || isempty(ggiw_mb1.r)
    ggiw_mb_catenate = ggiw_mb2;
else
    ggiw_mb_catenate.r = [ggiw_mb1.r;ggiw_mb2.r];
    ggiw_mb_catenate.alpha = [ggiw_mb1.alpha;ggiw_mb2.alpha];
    ggiw_mb_catenate.beta = [ggiw_mb1.beta;ggiw_mb2.beta];
    ggiw_mb_catenate.x = [ggiw_mb1.x ggiw_mb2.x];
    ggiw_mb_catenate.P = cat(3,ggiw_mb1.P,ggiw_mb2.P);
    ggiw_mb_catenate.v = [ggiw_mb1.v;ggiw_mb2.v];
    ggiw_mb_catenate.V = cat(3,ggiw_mb1.V,ggiw_mb2.V);
end
end