% Lucas Rath
% Julius Hiller

clear all; close all; clc;


%% Question 1

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
gamma = 0.9;

% allocate memory
V = nan(kmax,n);
u = zeros(n,1);

% init iteration variables
V(1,:) = 0;


k=1;
while k<kmax
    for s=1:n
        % V = min{.}
        [V(k+1,s), idx] = min( fo(s,:) + gamma * V(k,f(s,:)) );
        % u = argmin{.}
        u(s) = idx-1;
    end
    if norm(V(k+1,:)-V(k,:)) < epsilon
        break;
    end
    k = k+1;
end
disp(V)
disp(u')

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
while s~=8
    fprintf('State:%d, Optimal input:%d, Cost:%.1f\n',s,u(s),fo(s,u(s)+1))
    s  = f(s,u(s)+1);
end
fprintf('State:%d, Optimal input:%d, Cost:%.1f\n',s,u(s),fo(s,u(s)+1))




