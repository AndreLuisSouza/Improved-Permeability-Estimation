function y=multiweibullPDF(x,c) 
% Multi-Weibull probability density function (First Derivative)
% This function allows for variable number of Weibull functions

if mod(length(c),4)>0 % Module after division
    error('length(c) might be multiple of 4')
end

y=0*x; % Zeros 
for k=1:length(c)/4
    y = y + c(4*k-3) * c(4*k)/c(4*k-1) *  ((x-c(4*k-2))/c(4*k-1)).^(c(4*k)-1) .* ...
        exp(-( (x-c(4*k-2))./c(4*k-1)).^c(4*k)) .*((x-c(4*k-2))>0); 
    % c(4*k-3) = Amplitude (c1); 
    % c(4*k-2) = Entry (c2); 
    % c(4*k-1)) = Scale (c3);
    % c(4*k) = obliquity/shape (c4); 
    % ((x-c(4*k-2))>0) = eliminates negative valuesfrom each multiweibull    
end