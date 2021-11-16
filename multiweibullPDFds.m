function y=multiweibullPDFds(x,c) 
% Multi-Weibull probability density function (Second Derivative)
% This function allows for variable number of Weibull functions

if mod(length(c),4)>0 % Module after division
    error('length(c) might be multiple of 4')
end

y=0*x; % Zeros
for k=1:length(c)/4
    c1 = c(4*k-3);
    c2 = c(4*k-2);
    c3 = c(4*k-1);
    c4 = c(4*k);
    y = y + c1*c4/c3/c3*( ...
        (c4-1)*(((x-c2)/c3).^(c4-2)) - c4*(((x-c2)/c3).^(2*c4-2))) .* exp(-((x-c2)/c3).^c4) .* ((x-c2)>0); % PDF

    % c(4*k-3) = Amplitude (c1);
    % c(4*k-2) = Entry (c2);
    % c(4*k-1)) = Scale (c3);
    % c(4*k) = obliquity/shape (c4);
    % ((x-c(4*k-2))>0) = eliminates negative valuesfrom each multiweibull
end