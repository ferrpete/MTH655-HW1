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