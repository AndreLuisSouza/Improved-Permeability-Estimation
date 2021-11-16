function c=chuteinicialweibullCDF(x,y)  
% Define the number of summed Multi-Weibulls by clicking with the mouse left button
% When finished the sellection, click with the right button to start the
% optimization process

plot_fit_PcHg(x, y)
dydx = [0; diff(y)./diff(x)];
dydx(dydx>200)=NaN;
xclick=[];
while 1    
    [x0,~,but]=ginput(1); %Click at inflexion points with the left mouse button 
    if but==1
        xclick=[xclick;x0];
        yclick=interp1(x,y,xclick); % Interpolate values of "y" = (SHg) to each "x" = (Pc) selected
        c=zeros(4*length(xclick),1);
        for k=1:length(xclick)
            
               
            c(4*k)=3; % Obliquity (Beta) parameter specializing Weibull function to Gaussian distribution
            
            c(4*k-1)=std(x); % Scale (Eta) parameter (tortuosity)
            
            c(4*k-2)=xclick(k)-c(4*k-1)*(((c(4*k)-1)/c(4*k))^(1/c(4*k)));  % Position (ommega) parameter (Entry pressure)
            
            c(4*k-3)= yclick(k)*c(3)/c(4) ...
                      .* ((xclick(k)-c(4*k-2))/c(4*k-1)).^(1-c(4*k)) ...
                      .* exp(( (xclick(k)-c(4*k-2))./c(4*k-1)).^c(4*k)); % Volume (alpha) parameter 
        end
        c(1:4:end)=c(1:4:end)/sum(c(1:4:end))*100; % Normalizing volume to 100% 
        plot_fit_PcHg(x, y, c)      
    else 
        break
    end
end
