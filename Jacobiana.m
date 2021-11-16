function [ J ] = Jacobiana( x,f,c )
%UNTITLED Summary of this function goes here
%   Creating the Jacobian matrix 
%   parameters 'C' and the delta 'h'

h = mean(abs(diff(x)))/1000;
J = zeros(length(x),length(c));

for k = 1:length(c)
    delta = 0*c;
    delta(k) = h;
    J(:,k) = f(x,c+delta) - f(x,c-delta);
end
J = J/2/h;

