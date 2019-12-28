% Lucas Rath
% Julius Hiller

clear all; close all; clc;

%% config

saveimages = false;
showtitles = true;


if saveimages
    saveimg = @(fig,name,format) fp.savefig(fig,name,format);
else, saveimg = @(fig,name,format) 0; end
if showtitles
    stitle = @(text) title(text);
else, stitle = @(text) 0; end


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
kmax = 10;
% discount factor
gamma = 0.9;


V = zeros(kmax,n);
for k=1:kmax   
    % V(k+1,:) = min( fo + gamma * reshape(V(k,f),n,m), [],2);
    for s=1:n
        % V = min{.}
        [V(k+1,s), idx] = min( fo(s,:) + gamma * V(k,f(s,:)) );
        % u = argmin{.}
        u(s) = idx-1;
    end
end

% show converged value function
V(end,:)
figure('Color','w'); hold on; grid on;
plot(0:kmax,V)

% ---------------------------------------------------------------------
%  Calculate optimal input sequence when starting at s=1
% ---------------------------------------------------------------------
s = 1;
while s~=8
    fprintf('State:%d, Optimal input:%d, Cost:%.1f\n',s,u(s),fo(s,u(s)+1))
    s  = f(s,u(s)+1);
end
fprintf('State:%d, Optimal input:%d, Cost:%.1f\n',s,u(s),fo(s,u(s)+1))




