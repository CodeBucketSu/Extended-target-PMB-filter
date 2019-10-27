function [crossEntropy] = crossEntropyGamma(ah,bh,ai,bi)
% Compute cross entropy between two Gammas
crossEntropy = -ai*log(bi) + gammaln(ai) - (ai-1)*(psi(0,ah)-log(bh)) + bi*ah/bh;
end