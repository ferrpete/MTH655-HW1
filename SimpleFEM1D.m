% SimpleFEM1D.m
% Peter Ferrero, Oregon State University, MTH655, 1/31/2018
% A simple FEM 1D method to solve Problem 4 of Homework 1 for MTH 655.

function [FemSol, x] = SimpleFEM1D(h)

a = 0; % left endpoint
b = 1; % right endpoint
x = a:h:b; % mesh nodes
N = length(x); % number of nodes

A = GStiff(x); % Global Stiffness matrix
F = GLoad(x); % Load vector

FemSol = zeros(N, 1); % Initialize the FEM solution
FemSol(2:N-1) = A(2:N-1,2:N-1)\F(2:N-1); % Solve the linear system
                                 % for interior nodes

end