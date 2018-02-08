% Main.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% The main file for the FEM 1D method to solve Problem 4 of Homework 1 for
% MTH 655.

clear all

n = [2,4,8];
N = length(n);
x = [0:0.001:1]';
ExactSol = Exact(x);

figure(1)
plot(x, ExactSol, 'k')
xlabel('x')
ylabel('u(x)')
hold on

for i = 1:N
    
    [FemSol, x] = SimpleFEM1D(n(i));
    ExactSol = Exact(x');
    error(i) = norm(ExactSol-FemSol, inf);
    plot(x,FemSol,'-o')
    
end

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)./n; % uniform mesh size

legend('Exact', 'h = 0.5', 'h = 0.25', 'h = 0.125')
hold off

figure(2)
loglog(h,h,'b-',h,h.^2,'k-',h,error,'*-r')
xlabel('Mesh size, h')
ylabel('Max Norm Error, e(h)')
legend('Linear', 'Quadratic', 'Linear FEM Error')

% SimpleFEM1D.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% A simple FEM 1D method to solve Problem 4 of Homework 1 for MTH 655.

function [FemSol, x] = SimpleFEM1D(N)

a = 0; % left endpoint
b = 1; % right endpoint
h = (b-a)/N; % uniform mesh size
x = a:h:b; % mesh nodes

A = GStiff(x); % Global Stiffness matrix
F = GLoad(x); % Load vector

FemSol = zeros(N+1, 1); % Initialize the FEM solution
FemSol(2:N) = A(2:N,2:N)\F(2:N); % Solve the linear system
                                 % for interior nodes

end

% GStiff.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the 1D global stiffness matrix for FEM.

function A = GStiff(x)

% ===Input: vector x of mesh nodes===
% ===Output: Assembled stiffness matrix A===

N = length(x)-1; % number of elements
A = zeros(N+1,N+1); % initialize the stiffness matrix to zero

for i = 1:N % loop over elements
    h = x(i+1)-x(i); % element length
    n = [i i+1];
    A(n,n) = A(n,n) + [1 -1; -1 1]/h; % incorporate local stiffness
                                       % matrix into global matrix
end

end

% GLoad.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the 1D load vector for FEM.

function F = GLoad(x)

% ===Input: vector x of mesh nodes===
% ===Output: Load vector F using Trapezoidal Rule===

N = length(x)-1;
F = zeros(N+1,1);

for i = 1:N
    h = x(i+1) - x(i);
    n = [i, i+1];
    F(n) = F(n) + [feval(@Loadf,x(i)); feval(@Loadf,x(i+1))]*(h/2);
    
end

end

% Loadf.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the local 1D load vector for FEM.

function f = Loadf(x)

% ===Input = x, mesh nodes at which to evaluate the load function f
% ===Output = f, value of load f at x

f = (pi^2)*sin(pi*x);

end

% Exact.m
% Peter Ferrero, Oregon State University, MTH 655, 1/31/2018
% A function to compute the exact solution for the simple FEM method.

function y = Exact(x)

% ===Input = x, mesh nodes at which to evaluate the exact solution y
% ===Output = y, the exact solution at x

y = sin(pi*x);

end