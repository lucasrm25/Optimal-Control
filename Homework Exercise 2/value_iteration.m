% -------------------------------------------------------------------------
% Finite-time infinite-horizon optimal control - Value Function Iteration and LP
%
% Programmed by:
%   Lucas Rath
%   Julius Hiller
%
% -------------------------------------------------------------------------

clear all; close all; clc;


%% Value Function Iteration

fprintf('------------------------------------------------\n')
fprintf('               Value Function Iteration:\n')
fprintf('------------------------------------------------\n\n')

% ---------------------------------------------------------------------
% Define discrete-time infinite-horizon DP problem
% ---------------------------------------------------------------------

% number of inputs
m = 3;
% number of states
n = 8;

% cost matrix: fo = @(s,u)
fo = [3   3 1
      5   3 1
      6   6 5
      1   0 1
      3   3 2
      2.5 2 4
      1   1 1
      0   0 0];

% motion model: s_{k+1} = f(s_{k},u)
f = [2 2 3
     7 5 4
     4 6 5
     7 8 2
     4 4 6
     1 7 8
     8 8 8
     8 8 8];


% ---------------------------------------------------------------------
% Value function algorithm
% ---------------------------------------------------------------------

% max number of iterations
kmax = 100;
% residue
epsilon = 1e-3;
% discount factor
alpha = 0.9;

% allocate memory
V = nan(kmax,n);
u = zeros(1,n);

% value function initial guess
V(1,:) = 0;

k=1;
while k<kmax
    for s=1:n
        % V = min{TuV(s)}
        [V(k+1,s), idx] = min( fo(s,:) + alpha * V(k,f(s,:)) );
        % u = argmin{TuV(s)}
        u(s) = idx-1;
    end
    if norm(V(k+1,:)-V(k,:)) < epsilon
        break;
    end
    k = k+1;
end

fprintf('Optimal Value function (after %d iterations):\n',k);
disp(V(k+1,:))
fprintf('Optimal policy:\n');
disp(u)

% plot iterations
figure('Color','w'); hold on; grid on;
plot(0:kmax-1,V')
xlabel('iteration')
ylabel('Value function')
legend({'V(\xi_1)','V(\xi_2)','V(\xi_3)','V(\xi_4)','V(\xi_5)','V(\xi_6)','V(\xi_7)','V(\xi_8)'}, 'Location', 'Northwest')


% ---------------------------------------------------------------------
%  Calculate optimal input sequence when starting at s=1
% ---------------------------------------------------------------------
s = 1;
fprintf('Optimal input sequence when starting in state s=%d:\n',s);
while s~=8
    fprintf('\tState:%d, Optimal input:%d, Cost:%.1f\n',s,u(s),fo(s,u(s)+1))
    s  = f(s,u(s)+1);
end
fprintf('\tState:%d, Optimal input:%d, Cost:%.1f\n',s,u(s),fo(s,u(s)+1))




%% Value Function as Linear Programming

fprintf('\n\n')
fprintf('------------------------------------------------\n')
fprintf('     Value Function as Linear Programming:\n')
fprintf('------------------------------------------------\n\n')

% build system of linear inequalities
Vx = sym('Vx',[n,1],'real');
constraints = sym(nan(n*m,1));
i=1;
for xi=1:n
    for im=1:m
        constraints(i) = Vx(xi) == fo(xi,im) + alpha * Vx(f(xi,im));
        i=i+1;
    end
end
[A,b] = equationsToMatrix(constraints);
A = double(A);
b = double(b);
c = -ones(n,1);


% find optimal value function
[V_LP,fval,exitflag,output,lambda] = linprog(c,A,b); % solves min c'*x such that A*x ≤ b.
V_LP = V_LP';

% Calculate optimal policy, given optimal value function
for s=1:n
    [~, idx] = min( fo(s,:) + alpha * V_LP(f(s,:)) );
    % u = argmin{TuV(s)}
    u_LP(s) = idx-1;
end


fprintf('Optimal Value function (Linear Programming approach):\n');
disp(V_LP)
fprintf('Optimal policy:\n');
disp(u_LP)

