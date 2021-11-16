function plot_fit_PcHg(x, y, c)
% Plot Raw Data + First Guess + Fit Multi-Weibull

%r = (-((((2*480*cos(2*pi*140/360))./x)*1450.377)/10000));

if nargin<3
    c = [];
end
carrega_cores
clf

% Pc x SHg (raw data)
subplot(211)
plot(x, y, 'o', 'color', cicloPcHg(1,:))
hold on
box on
grid on
% try 
%     set(gca,'gridalpha',1) % Old Matlab version do not support this
%     command
% end
%set(gca,'gridcolor',[0.75, 0.75, 0.75], 'ticklabelinterpreter','latex')
xlabel('$\log_{~10}(P_c, \textrm{psi})$','interpreter','latex','fontsize',10)
ylabel('$S_{Hg}$ (\%)','interpreter','latex','fontsize',10)
title('Weibull','interpreter','latex','fontsize',10)
xlim([min(x) max(x)])
ylim([0 100])
set(gca, 'FontName', 'Times New Roman')

%Pc x derivada
subplot(212)
dydx = [0; diff(y)./diff(x)];
dydx(dydx>200)=NaN;
plot(x,dydx,'o','color', cicloPcHg(1,:))
hold on
box on
grid on
% axis tight
% try 
%     set(gca,'gridalpha',1) %Versões antigas não tem esse comando
% end
%set(gca,'gridcolor',[0.75, 0.75, 0.75], 'ticklabelinterpreter','latex')
%set(gca,'gridcolor',[0.75, 0.75, 0.75])
xlabel('$\log_{~10}(P_c, \textrm{psi})$','interpreter','latex','fontsize',10)
ylabel('$\displaystyle\frac{d S_{Hg}}{d\log_{~10}P_{c}}$ (\%/psi)','interpreter','latex','fontsize',10)
xlim([min(x) max(x)])
% ylim([0 100])
set(gca, 'FontName', 'Times New Roman')
%Curvas do fit
xfit = interp(x, 20);
for k=1:length(c)/4
    subplot(211)
    plot(xfit, multiweibullCDF(xfit, c(4*k-3:4*k)), ...
        'linewidth', 1.5, 'color', cicloPcHg(k+2,:));
    subplot(212)
    plot(xfit, multiweibullPDF(xfit, c(4*k-3:4*k)), ...
        'linewidth', 1.5, 'color', cicloPcHg(k+2,:));
end

if ~isempty(c)
    subplot(211)
    plot(xfit, multiweibullCDF(xfit, c(c~=0)), '--', 'linewidth', 2, 'color', cicloPcHg(2,:));
    subplot(212)
    plot(xfit, multiweibullPDF(xfit, c(c~=0)), '--', 'linewidth', 2, 'color', cicloPcHg(2,:));
    pause(0.1)
end


