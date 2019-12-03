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




%% User input


% cost function parameters
Q = eye(2);
R = 1;

% horizon length
tf = 10;            % horizon length
N = 50;             % number of steps
h = tf/N;           % step size

% continuous dynamic equation parameters
x0 = [-1 1]';
Ac = [1 2; 
      1 1];
Bc = [1;
      2];

n = length(Ac);     % state dimension
m = size(Bc,2);     % input dimension


%% Discretize model

% exact discretization      - c2d(ss(Ac,Bc,eye(n),0), h, 'zoh')
Ad_1 = expm(Ac*h);
syms t
Bd_1 = double( int(expm(Ac*t), t, 0, h) * Bc );


% Euler discretization
Ad_2 = eye(n) + h*Ac;
Bd_2 = Bc*h;


% compare results
Ad_1, Ad_2
Bd_1, Bd_2


%% Question F) Lets verify the difference in the solution when using the two different discretization methods

% Exact discretization (Ad_1,Bd_1)
[H, f, d, Aeq, beq]        = define_opt_control_problem (Ad_1, Bd_1, Q, R, x0, N, n, m);
[x_opt_1, u_opt_1]         = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);
[x_opt_1_KKT, u_opt_1_KKT] = solve_opt_control_KKT      (H, f, d, Aeq, beq, N, n, m);

% Euler discretization (Ad_2,Bd_2)
[H, f, d, Aeq, beq]        = define_opt_control_problem (Ad_2, Bd_2, Q, R, x0, N, n, m);
[x_opt_2, u_opt_2]         = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);
[x_opt_2_KKT, u_opt_2_KKT] = solve_opt_control_KKT      (H, f, d, Aeq, beq, N, n, m);


norm(x_opt_1-x_opt_1_KKT)
norm(x_opt_1-x_opt_2)



fig = figure('Color','white','Position',[262  317  612  420]); hold on; grid on;
clr = lines(20);
stairs(0:h:tf-h, x_opt_1(1,:)', '-','Color',clr(1,:),'LineWidth',2,'DisplayName','x(1)_{quadprog}');
stairs(0:h:tf-h, x_opt_1(2,:)', '-','Color',clr(2,:),'LineWidth',2,'DisplayName','x(2)_{quadprog}');
stairs(0:h:tf-h, u_opt_1(1,:)', '-','Color',clr(3,:),'LineWidth',2,'DisplayName','u_{quadprog}'); 
stairs(0:h:tf-h, x_opt_1_KKT(1,:)', '--','Color',clr(4,:),'LineWidth',2,'DisplayName','x(1)_{KKT}');
stairs(0:h:tf-h, x_opt_1_KKT(2,:)', '--','Color',clr(5,:),'LineWidth',2,'DisplayName','x(2)_{KKT}');
stairs(0:h:tf-h, u_opt_1_KKT(1,:)', '--','Color',clr(6,:),'LineWidth',2,'DisplayName','u_{KKT}'); 
legend
xlabel('time [s]')
stitle('Comparison quadprog solver and KKT solution (for exact discretization)')
saveimg(fig, 'question-f-quadprog-KKT', 'jpg');


fig = figure('Color','white','Position',[877  305  612  420]); hold on; grid on;
clr = lines(20);
stairs(0:h:tf-h, x_opt_1(1,:)', '-','Color',clr(1,:),'LineWidth',2,'DisplayName','x(1)_{exact}');
stairs(0:h:tf-h, x_opt_1(2,:)', '-','Color',clr(2,:),'LineWidth',2,'DisplayName','x(2)_{exact}');
stairs(0:h:tf-h, u_opt_1(1,:)', '-','Color',clr(3,:),'LineWidth',2,'DisplayName','u_{exact}'); 
stairs(0:h:tf-h, x_opt_2(1,:)', '--','Color',clr(4,:),'LineWidth',2,'DisplayName','x(1)_{Euler}');
stairs(0:h:tf-h, x_opt_2(2,:)', '--','Color',clr(5,:),'LineWidth',2,'DisplayName','x(2)_{Euler}');
stairs(0:h:tf-h, u_opt_2(1,:)', '--','Color',clr(6,:),'LineWidth',2,'DisplayName','u_{Euler}'); 
legend
xlabel('time [s]')
stitle('Comparison exact and Euler discretizations (using quadprog solver)')
saveimg(fig, 'question-f-disc', 'jpg');


fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0*0, x_opt_1-x_opt_2, u_opt_1-u_opt_2, h, tf)
title({'Difference of optimal states and control inputs',
       'for exact and euler discretization methods',
       'using quadprog solver'});
saveimg(fig, 'question-f-error', 'jpg');


% h_e = 0.02:0.02:1;
% err = 0*h_e;
% for i=1:numel(h_e)
%     % Euler discretization
%     N = floor(tf/h_e(i));
%     Ad_e = eye(n) + h_e(i)*Ac;
%     Bd_e = Bc*h_e(i);
%     [H, f, d, Aeq, beq] = define_opt_control_problem (Ad_e, Bd_e, Q, R, x0, N, n, m);
%     [x_opt_e, u_opt_e]  = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);
%     err(i) = norm(x_opt_1' - interp1(0:h_e(i):tf-h_e(i), x_opt_e', 0:h:tf-h) );
% end
% 
% fig = figure('Color','white','Position',[877  305  612  420]); hold on; grid on;
% clr = lines(20);
% plot(h_e, err, '-','Color',clr(1,:),'LineWidth',2);
% xlabel('Time step size h [s]')
% ylabel('Error norm')
% saveimg(fig, 'question-f-err', 'jpg');



%% Question G) What happens if we change Q?

[H, f, d, Aeq, beq] = define_opt_control_problem (Ad_1, Bd_1, 0.2*Q, R, x0, N, n, m);
[x_opt, u_opt]      = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt, u_opt, h, tf)
stitle('Opt. control solution using exact discretization (Q = 0.2*I_2)');
saveimg(fig, 'question-g-alpha-0p2', 'jpg');


[H, f, d, Aeq, beq] = define_opt_control_problem (Ad_1, Bd_1, 1*Q, R, x0, N, n, m);
[x_opt, u_opt]      = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt, u_opt, h, tf)
stitle('Opt. control solution using exact discretization (Q = 1*I_2)');
saveimg(fig, 'question-g-alpha-1p0', 'jpg');


[H, f, d, Aeq, beq] = define_opt_control_problem (Ad_1, Bd_1, 40*Q, R, x0, N, n, m);
[x_opt, u_opt]      = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m);

fig = figure('Color','w'); hold on; grid on;
plot_sim ( fig, x0, x_opt, u_opt, h, tf)
stitle('Opt. control solution using exact discretization (Q = 40*I_2)');
saveimg(fig, 'question-g-alpha-40p0', 'jpg');


%% Question H) Compare continuous-time and discrete-time solution using the optimal control input

u     = u_opt_1;
x_opt = x_opt_1;

sys = ss(Ac,Bc,eye(n),0);
[y,t,x] = lsim(sys, u_opt_1, 0:h:tf-h, x0, 'zoh');

% figure('Color','w'); hold on; grid on;
% plot(0:h:tf-h, x_opt, '-', 'LineWidth',2)
% plot(0:h:tf-h, x', '-.', 'LineWidth',2)
% saveimg(fig, 'question-h-discretization-error', 'jpg');


fig = figure('Color','white','Position',[322  392  725  282]);
hold on; grid on;
plot ( 0:h:tf-h, x_opt-x', 'o-','LineWidth',2)
ylabel('Error')
xlabel('Discretization time points t_k')
legend('error x(1)','error x(2)')
stitle('Error between continuous time and discrete time simulation');
saveimg(fig, 'question-h-discretization-error', 'jpg');


% close all;



%% Helper functions


function [H, f, d, Aeq, beq] = define_opt_control_problem ( Ad, Bd, Q, R, x0, N, n, m)
    %% Define Optimal control matrices

    gamma = zeros(N*n,N*m);
    for i=0:N-2
        gamma = gamma + kron(diag(ones(N-i-1,1),-i-1),Ad^i*Bd);     % fill lower-diagonals with Ad^i*Bd
    end
    
    omega = eye(n);
    for i=1:N-1
        omega = [omega; Ad^i];                              % extend matrix
    end
    
    Qb = kron(eye(N),Q);
    Rb = kron(eye(N),R);                        % Rb = blkdiag(R, R, ..., R)
    
    %% Define quadratic minimization problem

    % J = (xb'*Qb*xb + ub'*Rb*ub + x0'*Q*x0)    =   0.5*[xb,ub]'*H*[xb,ub] + f'*[xb,ub]' + d
    %   where:
    H = blkdiag(Qb,Rb);
    f = zeros(N*(n+m),1);
    d = 0;

    % equality constraints:  
    %  xb = omega*x0 + gamma*ub    =>   Aeq*[xb,ub] = beq
    %   where:
    Aeq = [ eye(n*N), -gamma ];  % Aeq = [ eye(n*N), -gamma ];
    beq = omega * x0;
end


function [x_opt, u_opt] = solve_opt_control_quadprog (H, f, d, Aeq, beq, N, n, m)
    %% Solve optimal control problem numerically

    % user input
    tol  = 1e-12;
    maxiter = 200;
    
    options = optimoptions('quadprog',...
                           'Algorithm','interior-point-convex',...
                           'ConstraintTolerance',tol,...
                           'MaxIterations', maxiter,...
                           'Display','iter');
    [Y,~,~] = quadprog(H, f, [],[], Aeq,beq, [],[], [],options);
    
    % check constraints error: (Aeq*Y-beq)'

    % reshape solution
    x_opt = reshape( Y(1:N*n), n, []);
    u_opt = reshape( Y(N*n+1:end), m, []);
end

function [x_opt, u_opt] = solve_opt_control_KKT (H, f, d, Aeq, beq, N, n, m)
    %% Solve optimal control problem using KKT conditions
    
    % 1 step solution
    if rank(H) == size(H,1)
        nu = -(Aeq/H*Aeq')\(Aeq/H*f + beq);
        Y  = -H\(Aeq'*nu + f);
        
    % 1 iterative solution - Newton's method
    else    
        % user input
        tol  = 1e-12;
        maxiter = 20;

        Y = zeros((n+m)*N,1);
        nu= zeros(n*N,1);

        iter  = 1;
        delta = Inf;
        while iter<=maxiter && norm(delta) >= tol
            A = [ 0.5*(H'+H)  Aeq';
                  Aeq         zeros(n*N,n*N)];
            b = -[0.5*(H'+H)*Y+f+Aeq'*nu;
                  Aeq*Y-beq];
            delta = A \ b;
            Y = Y  + delta(1:N*(n+m));
            nu= nu + delta(N*(n+m)+1:end);
            fprintf('iter: %4d, tolerance:%f\n',iter,norm(delta));
            iter = iter+1;
        end
    end
    
    % reshape solution
    x_opt = reshape( Y(1:N*n), n, []);
    u_opt = reshape( Y(N*n+1:end), m, []);
end


function plot_sim ( fig, x0, x_opt, u_opt, h, tf )
    %% Show results
    figure(fig);
    clr = lines(20);
    % subplot(3,1,1)
    stairs(0:h:tf-h, x_opt(1,:)', '-','Color',clr(1,:),'LineWidth',2,'DisplayName','x(1)');
    stairs(0:h:tf-h, x_opt(2,:)', '-','Color',clr(2,:),'LineWidth',2,'DisplayName','x(2)');
    stairs(0:h:tf-h, u_opt(1,:)', '-','Color',clr(3,:),'LineWidth',2,'DisplayName','u'); 

    legend
    xlabel('time [s]')
end

