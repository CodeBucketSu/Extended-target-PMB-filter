function [ est ] = state_extract( ggiw_mb )

d = 2;
% select Bernoulli density with existence probability larger than 0.5
ss = ggiw_mb.r >= 0.5;
est.x = ggiw_mb.x(:,ss);
v = ggiw_mb.v(ss);
V = ggiw_mb.V(:,:,ss);

for i = 1:length(v)
    est.X(:,:,i) = V(:,:,i)./(v(i)-2*d-2);
end
if isempty(v)
    est.X = zeros(2,2,0);
end

end