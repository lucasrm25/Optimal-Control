% -------------------------------------------------------------------------
% Model Predictive Control example
%
% Programmed by:
%   Lucas Rath
%   Julius Hiller
%
% -------------------------------------------------------------------------

clear all; close all; clc;


% -------------------------------------------------------------------------
% Analyse discrete-time system
% -------------------------------------------------------------------------


% xk+1 = A*xk + B*uk
A = [1 2; -1 3];
B = [0;1];

n = size(A,1);
m = size(B,2);

% determine equilibrium of unforced system (xk = A xk)
null(eye(2)-A)

% check if system is stable
eig(A)
abs(eig(A))




% -------------------------------------------------------------------------
% Formulate MPC problem
% -------------------------------------------------------------------------

Q = eye(2);
R = 1.5;
K = [-0.91,2.85];
P = [10.65 16.02
     16.02 67.01];
c = min(eig(P))/norm(K)^2;


% [P, ~, K] = dare(A,B,Q,R);
Mdare = (A-B*K)'* P *(A-B*K) - P + Q + K'*R*K
eig(Mdare)


% time horizon
N = 2;


% -------------------------------------------------------------------------
% Define MPC problem (self developed general MPC class)
% -------------------------------------------------------------------------


% define cost functions
fo   = @(t,mu_x,var_x,u,e,r) mu_x'*Q*mu_x + u'*R*u;
fend = @(t,mu_x,var_x,e,r)   mu_x'*P*mu_x;   % end cost function

% define dynamics
f  = @(mu_x,var_x,u) deal( A*mu_x + B*u, A*var_x*A');
% define additional constraints
h    = @(x,u,e) [];
hend = @(x,u,e) [];
g    = @(x,u,e) [];
gend = @(x,u,e) x'*P*x - c;
u_lb = -1;
u_ub =  1;
dt = 0.1;

% Initialize NMPC object;
mpc = NMPC(f, h, hend, g, gend, u_lb, u_ub, n, m, 0, fo, fend, N, dt);
mpc.tol     = 1e-2;
mpc.maxiter = 200;



% -------------------------------------------------------------------------
% Simulate MPC
% -------------------------------------------------------------------------

% initial state
x0 = [0.3; -0.25];
% x0 = [1; -0.9];

% define simulation time
kmax = 30;                  % steps to simulate

% initialize variables to store simulation results
out.x_mpc = [x0 NaN(n, kmax)];     % true states
out.u_mpc =     NaN(m, kmax);      % applied input

for k=1:kmax
    % calculate optimal input
    [~,u_opt,~] = mpc.optimize(out.x_mpc(:,k), 0,0);
    out.u_mpc(:,k) = u_opt(:,1);
    
    % check constraints
    % [x_pred,~] = mpc.predictStateSequence(out.x(:,k), 0, u_opt);
    % if gend(x_pred(:,end),0,0) > 0 
    %     error('Problem is unfeasible.. .aborting');
    % end
    
    % simulate real system
    [out.x_mpc(:,k+1),~] = f(out.x_mpc(:,k),zeros(n),out.u_mpc(:,k));
end

% -------------------------------------------------------------------------
% Simulate DLQR
% -------------------------------------------------------------------------

% initial state
x0 = [0.3; -0.25];
% x0 = [1; -0.9];

[K_lqr,~,~] = dlqr(A,B,Q,R,[]); 

% define simulation time
kmax = 30;                  % steps to simulate

% initialize variables to store simulation results
out.x_lqr = [x0 NaN(n, kmax)];     % true states
out.u_lqr =     NaN(m, kmax);      % applied input

for k=1:kmax
    % calculate optimal input
    out.u_lqr(:,k) = - K_lqr * out.x_lqr(:,k);
    
    % simulate real system
    [out.x_lqr(:,k+1),~] = f(out.x_lqr(:,k),zeros(n),out.u_lqr(:,k));
end


% -------------------------------------------------------------------------
% Plot results
% -------------------------------------------------------------------------

close all;
figure('Color','white','Position',[449  493  424  292]); hold on; grid on;
plot(0:kmax, out.x_mpc', '-', 'LineWidth',1.5)
plot(0:kmax, out.x_lqr', '--','LineWidth',1.5)
legend({'x_1^{mpc}','x_2^{mpc}','x_1^{lqr}','x_2^{lqr}'})
xlabel('time step k')
xlim([-1,kmax+1])
fp.savefig(gcf,'mpc_x','jpg')

figure('Color','white','Position',[807  519  424  292]); hold on; grid on;
plot(0:kmax-1, out.u_mpc', '-',  'LineWidth',1.5)
plot(0:kmax-1, out.u_lqr', '--', 'LineWidth',1.5)
legend({'u_1^{mpc}', 'u_1^{lqr}'})
xlabel('time step k')
xlim([-1,kmax+1])
fp.savefig(gcf,'mpc_u','jpg')

figure('Color','white','Position',[807  519  424  292]); hold on; grid on;
axis square
plot(out.x_mpc(1,:), out.x_mpc(2,:),'-o', 'LineWidth',1.5, 'DisplayName','x_k^{mpc}')
plot(out.x_lqr(1,:), out.x_lqr(2,:),'-o', 'LineWidth',1.5, 'DisplayName','x_k^{lqr}')
plot(out.x_mpc(1,1), out.x_mpc(2,1), '*', 'MarkerFaceColor','r','LineWidth',1.5, 'DisplayName','x_0')
Xf = sigmaEllipse2D( 0, inv(P), sqrt(c), 100 ); % i=50; Xf(:,i)'*P*Xf(:,i)-c
p1 = patch(Xf(1,:),Xf(2,:),'blue','FaceAlpha',.2, 'EdgeColor','none');     % plot(Xf(1,:),Xf(2,:));
p1.DisplayName = 'X_f';
legend('Location','bestoutside')
xlabel('x_1')
ylabel('x_2')
fp.savefig(gcf,'mpc_x-xy','jpg')



% -------------------------------------------------------------------------
% 
% -------------------------------------------------------------------------

% function [nonl_ineq, nonl_eq] = nonlcon(z)
%     
%     d = 1.3325;
%     P = [4.2 7;7 36.1];
%     T = blkdiag(zeros(2),zeros(2),P,0,0,0);
% 
%     nonl_ineq = z'*T*z-d; %nonlinear ineqality
%     nonl_eq = 0; %nonlinear equality
% 
% end


% % % % % -------------------------------------------------------------------------
% % % % % Define MPC problem
% % % % % -------------------------------------------------------------------------
% % % % 
% % % % 
% % % % 
% % % % gamma = kron(eye(N),B);
% % % % omega = A;
% % % % for i=1:N-1
% % % %     gamma = gamma + kron(diag(ones(N-i,1),-i),A^i*B);
% % % %     omega = [omega; A^(i+1)];
% % % % end
% % % % Qb = blkdiag( kron(eye(N-1),Q), P );
% % % % Rb = kron(eye(N),R);
% % % % 
% % % % gamma = zeros(N*n,N*m);
% % % % for i=0:N-1
% % % %     gamma = gamma + kron(diag(ones(N-i,1),-i),A^i*B);     % fill lower-diagonals with Ad^i*Bd
% % % % end
% % % % gamma = [zeros(n,N*m); gamma];
% % % % 
% % % % omega = eye(n);
% % % % for i=1:N
% % % %     omega = [omega; A^i];                              % extend matrix
% % % % end
% % % % 
% % % % Qb = blkdiag(kron(eye(N),Q), P);
% % % % Rb = kron(eye(N),R);
% % % % 
% % % % 
% % % % 
% % % % % Define quadratic minimization problem
% % % % 
% % % % % J = (xb'*Qb*xb + ub'*Rb*ub + xN'*Qf*xN)    =   0.5*[xb,ub]'*H*[xb,ub] + f'*[xb,ub]' + d
% % % % %   where:
% % % % H = blkdiag(Qb,Rb);
% % % % f = zeros(n+N*(n+m),1);
% % % % 
% % % % % equality constraints:  
% % % % %  xb = omega*x0 + gamma*ub    =>   Aeq*[xb,ub] = beq
% % % % %   where:
% % % % Aeq = [ eye(n*(N+1)), -gamma ];
% % % % beq = omega * x0;
% % % % 
% % % % 
% % % % % inequality constraints: |uk|<1
% % % % uvar = sym('u',[N,m])
% % % % xvar = sym('x',[N+1,n])
% % % % zvar = [xvar(:); uvar(:)]
% % % % ineqconstr = [];
% % % % for i=1:N
% % % %     %   |uk| <= 1
% % % %     ineqconstr = [ineqconstr; uvar(i,:) < 1; -uvar(i,:) < 1 ];
% % % % end
% % % % [Aineq,bineq] = equationsToMatrix( lhs(ineqconstr)-rhs(ineqconstr),zvar);
% % % % 
% % % % 
% % % % % define quadratic inequality constraint
% % % % T = zeros((N+1)*n+N*m);
% % % % T(N+n+1:N+n+n,N+n+1:N+n+n) = P;
% % % % d = c;



% % % % % -------------------------------------------------------------------------
% % % % % Simulate
% % % % % -------------------------------------------------------------------------
% % % % 
% % % % % user input
% % % % tol  = 1e-12;
% % % % maxiter = 200;
% % % % 
% % % % options = optimoptions('quadprog',...
% % % %                        'Algorithm','interior-point-convex',...
% % % %                        'ConstraintTolerance',tol,...
% % % %                        'MaxIterations', maxiter,...
% % % %                        'Display','iter');
% % % % [Y,~,~] = fmincon( @(z)z'*H*z ,x0,A,b,Aeq,beq,lb,ub,nonlcon,options)
% % % % 
% % % % % check constraints error: (Aeq*Y-beq)'
% % % % 
% % % % % reshape solution
% % % % x_opt = reshape( Y(1:N*n), n, []);
% % % % u_opt = reshape( Y(N*n+1:end), m, []);
