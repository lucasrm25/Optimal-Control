clear all; close all; clc;



syms x1 x2 real

f1 = 0.5 * (x1^2 + x2^2);
g1 = x1 + x2 + 2;           % <=0

f2 = x1 + x2;
h2 = [ (x1-1)^2+x2^2-1;
       (x1-2)^2+x2^2-4 ];  % == 0 
   
   
%% 1) A1

fig = figure(1); hold on; grid on; fig.Color = 'w';
c1 = plot([0 -2],[-2 0]);
c1.DisplayName = 'x1 + x2 + 2 <= 0';
c1.LineWidth = 2;
c1.Color = 'g';

p1 = patch([-3,-2,0,0,-3],[0,0,-2,-3,-3],'green','FaceAlpha',.3, 'EdgeColor','none');
p1.DisplayName = 'x1 + x2 + 2 <= 0';


[X1,X2] = meshgrid(-3:0.05:0 , -3:0.05:0);
F = 0.5 * (X1.^2 + X2.^2);
% s = surf(X1,X2,F,'FaceAlpha',0.5, 'EdgeColor', 'none', 'DisplayName','f1 = 0.5 * (x1^2 + x2^2)');
% c = colorbar;
% c.Label.String = 'f2 = x1 + x2'
[~,cnt] = contour(X1,X2,F,'ShowText','on');
cnt.DisplayName  = '0.5 * (x1^2 + x2^2)';
axis equal
legend([cnt,p1])
xlabel('x1');
ylabel('x2');
saveas(fig, 'q1-a1', 'jpg')

%% 1) A2

fig = figure(2); hold on; grid on; fig.Color = 'w';
c1 = circle(1,0,1);
c1.DisplayName = '(x1-1)^2+x2^2-1 = 0';
c1.LineWidth = 2;
c2 = circle(2,0,2);
c2.DisplayName = '(x1-2)^2+x2^2-4 = 0';
c2.LineWidth = 2;
xlabel('x1'); ylabel('x2');

[X1,X2] = meshgrid(-0.5:0.05:4.5 , -2.5:0.05:2.5);
F = X1 + X2;
% s = surf(X1,X2,F,'FaceAlpha',0.5, 'EdgeColor', 'none', 'DisplayName','f2 = x1 + x2');
% c = colorbar;
% c.Label.String = 'f2 = x1 + x2'
[~,cnt] = contour(X1,X2,F,'ShowText','on');
cnt.DisplayName  = 'x1 + x2';
axis equal
legend([cnt,c1,c2])
xlabel('x1');
ylabel('x2');
saveas(fig, 'q1-a2', 'jpg')



%% 1) C1

f = @(x) 0.5 * (x(1)^2 + x(2)^2);
% g1 = x1 + x2 + 2;    % <=0       
A = [1 1];
B = -2;
x0 = [-2 -2]';
[xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(f,x0,A,B)



%% 1) C2
f = @(x) x(1) + x(2);
h = @(x) [ (x(1)-1)^2+x(2)^2-1;
           (x(1)-2)^2+x(2)^2-4 ];  % == 0 
x0 = [-2 -2]';
[xmin,fval,exitflag,output,lambda,grad,hessian] = fmincon(f,x0,[],[],[],[],[],[],@hfun)


%% 2) C)

load('dataset.mat')

xfun = @(x) [1, x(1), x(1)^2, x(1)^3];
X = cell2mat( arrayfun( xfun, x , 'UniformOutput', false) );

beta = (X'*X)\X'*y;

xfit = (0:0.1:10)';
Xfit = cell2mat( arrayfun( xfun, xfit , 'UniformOutput', false) );
yfit = Xfit * beta;

fig = figure(3); hold on; grid on; fig.Color = 'w';
p1 = scatter(x,y)
p2 = plot(xfit,yfit)

p2.LineWidth = 2;
p1.DisplayName = 'Data';
p2.DisplayName = 'fitted curve: f(x) = \beta_1 + \beta_2*x + \beta_3*x^2 + \beta_4*x^3'
legend([p1,p2])
xlabel('x'); ylabel('y');
saveas(fig, 'q2', 'jpg')


%% helper


function [hineq, heq] = hfun(x)
    hineq = [];
    heq = [ (x(1)-1)^2+x(2)^2-1;
            (x(1)-2)^2+x(2)^2-4 ];
end

function p = circle(x,y,r)
    ang=0:0.01:2*pi; 
    xp=r*cos(ang);
    yp=r*sin(ang);
    p = plot(x+xp,y+yp);
end