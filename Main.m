% Main.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% The main file for the FEM 1D method to solve Problem 4 of Homework 1 for
% MTH 655.

h = 2.^(-[1;2;3]);
N = length(h);
x = [0:0.001:1]';
ExactSol = Exact(x);

figure(1)
plot(x, ExactSol, 'k')
xlabel('x')
ylabel('u(x)')
hold on

for i = 1:N
    
    [FemSol, x] = SimpleFEM1D(h(i));
    ExactSol = Exact(x');
    error(i) = norm(FemSol-ExactSol, inf);
    plot(x,FemSol,'-o')
    
end

legend('Exact', 'h = 0.5', 'h = 0.25', 'h = 0.125')
hold off

figure(2)
loglog(h,h.^2,'k-',h,error,'*-r')
xlabel('Mesh size, h')
ylabel('Max Norm Error, e(h)')
legend('Quadratic', 'Linear FEM Error')