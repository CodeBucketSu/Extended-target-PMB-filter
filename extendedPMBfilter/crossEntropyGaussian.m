function [crossEntropy] = crossEntropyGaussian(xh,Ph,xi,Pi)
% Compute cross entropy between two Gaussians
temp = xh - xi;
inv_Pi = inv(Pi);
crossEntropy = (trace(inv_Pi*Ph)+temp'*inv_Pi*temp+log(det(2*pi*Pi)))/2;
end