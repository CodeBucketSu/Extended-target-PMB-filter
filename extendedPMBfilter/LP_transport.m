function [ Cmin,q ] = LP_transport( C,ph )

[H,N] = size(C);
c = reshape(C',H*N,1);
beq = [ph;ones(N,1)];
Aeq1 = kron(eye(H),ones(1,N));
Aeq2 = repmat(eye(N),1,H);
% replace inf with a very large number; needed by the LP solver
c(isinf(c)) = max(c(~isinf(c)))+1e12;
Aeq = [Aeq1;Aeq2];

options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
[x,Cmin] = linprog(c,[],[],Aeq,beq,zeros(1,H*N),ones(1,H*N),options);
q = reshape(x,N,H)';

end

