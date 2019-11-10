function [ Cmin,q ] = LP_transport( C,ph )

[H,N] = size(C);
c = reshape(C',H*N,1);
beq = [ph;ones(N,1)];
Aeq1 = kron(eye(H),ones(1,N));
Aeq2 = repmat(eye(N),1,H);

Aeq = [Aeq1;Aeq2];

% can be largely speed up by using the Gurobi solver
options = optimoptions('linprog','Algorithm','dual-simplex','Display','off');
[x,Cmin] = linprog(c,[],[],sparse(Aeq),beq,zeros(1,H*N),ones(1,H*N),options);

% params.outputflag = 0;          % Silence gurobi
% params.method     = 1;          % Use dual simplex method
% model_gurobi.A = sparse(Aeq);
% model_gurobi.obj = c;
% model_gurobi.sense = '=';
% model_gurobi.rhs = beq;
% result = gurobi(model_gurobi, params);
% x = result.x;
% Cmin = result.objval;
q = reshape(x,N,H)';

end

