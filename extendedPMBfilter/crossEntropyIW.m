function [crossEntropy] = crossEntropyIW(vh,Vh,vi,Vi)
% Compute cross entropy bewteen two inverse Wisharts
d = 2;
crossEntropy = (vi-d-1)*d/2*log(2)-(vi-d-1)/2*log(det(Vi))+...
    logamma2((vi-d-1)/2)+vi/2*(log(det(Vh))-d*log(2)-psi(0,(vh-d-1)/2)...
    -psi(0,(vh-d-2)/2))+trace((vh-d-1)*inv(Vh)*Vi)/2;
end