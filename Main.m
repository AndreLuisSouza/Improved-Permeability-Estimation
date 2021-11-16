%% Multi-Weibull-Gaussian first guess followed by fit Multi-Weibull (4 free parameters) Pc x Shg;
clc
clear

% Cumulative Weibull function (cdf_wb):
% cdf_wb = @(x,c) c(1) * ( 1-exp(-((x-c(2))/c(3)).^c(4))) .* ((x-c(2))>0);  --->>> ((x-c(2))>0) eliminates negative values from each Weibull;

% Probability Weibull distribution (pdf_wb'):
% pdf_wb' = @(x,c) c(1) * c(4)/c(3) * ((x-c(2))/c(3)).^(c(4)-1).*exp(-(x-c(2))/c(3)).^c(4) .* ((x-c(2))>0); 

% Probability Weibull distribution (pdf_wb''): 
% pdf_wb'' = @(x,c) c1*c4/c3/c3*((c4-1)*(((x-c2)/c3).^(c4-2)) - c4*(((x-c2)/c3).^(2*c4-2))) .* exp(-((x-c2)/c3).^c4) .* ((x-c2)>0); 

% Parameterss:
% A (c(1)) = partitioned volume fractions (%);
% B (c(2)) = position factor (entry pressure - log10(psi));
% C (c(3)) = dispersion factor (log10(psi));
% D (c(4)) = obliquity/shape parameter (geometry (adm));

pchgm = importdata('MICP_data.txt'); % Dados PcHg de MERO
%Poço: 2-ANP-2A-RJS
depth = pchgm(:,1);
perm = pchgm(:,2);
phi = pchgm(:,3)/100;
pc = pchgm(:,4); %*********************************************************
pc_oa = pc.* -(35*cos(35*2*pi/360))./(480*cos(140*2*pi/360)); % System convertion O/A
sat = pchgm(:,5); %******************************************************
shg = pchgm(:,10);
fzi = (0.0314*(perm./(phi)).^(0.5))./(phi./(1-phi)); % Flow Zone Indicator

%%
maxit = 1000; % Max number of iterations
tol = 1e-12; % Error tolerance
a=1; % Gradient first guess (Steepest Descent)

ID = unique(depth); % Identify the samples
npar_wb = 28; % 4 (param) x N Weibull's - determines the resulting fit vector length Ex.: 7 Weibull x 4 parameters = 28;
c_wb = zeros(length(ID),npar_wb); % Pre alocating the vector of parameters
pfit_wb = zeros(length(ID),npar_wb); % Pre alocating a matrix of multiple adjusted parameters
R2_wb = zeros(length(ID),1); % Pre alocating the goodness-of-fit R2
vet = logspace(-3,5,300); % Sinthetic vector of pressure (psi)
rad = (-((((2*480*cos(2*pi*140/360))./pc)*1450.377)/10000)); % Vector of pore throats (mum)
phi_samples = zeros(length(ID),1); % Pre alocating a vector of samples porosity
perm_samples = zeros(length(ID),1); % Pre alocating a vector of samples permeability

ind=0;
while ind<=length(ID) 
    ind = ind+1;
    
    posicoes = find(depth==ID(ind)); % Selects the interval of the data
    
    x = pc(posicoes); % Pressure vector
    x1 = rad(posicoes); % Pore throat vector
    y = shg(posicoes); % Hg saturation vector
    dydx = [0;diff(y)./diff(log10(x))]; % Finite difference Hg saturation vector
    
    [phi_samples(ind,:)] = unique(phi(posicoes)); 
    [perm_samples(ind,:)] = unique(perm(posicoes)); 
    
    m_wb = zeros(length(ID),length(x));
    m_data = zeros(length(ID),length(x));
     
    figure(1)
    clf
    set(gcf, 'name',['Sample ', num2str(ind), '_First Guess'])
    c0_wb=chuteinicialweibullCDF(log10(x),y); 
    % Plot Raw curve and the finite difference first derivative
    % Ginput Function:
    % Select first-guess parameters by clicking left
    % button at disired positions. Click the right button when finshed to start the optimization process. 
    % Algorithm: Levemberg-Marquadt; Convergence/Inversion: Steepest Descent + Gauss-Newton
    % Constraints on c(1,2,3,4) 
    % Visit code: fitmultiweibullCDF

    c1_wb=fitmultiweibullCDF(log10(x),y,c0_wb,maxit,tol); 
    % Plot Fit Multi-Weibull
    % Creates figure 2 with the resulting fit
    
    figure(2)
    clf
    set(gcf, 'name',['Sample ', num2str(ind), 'Fit Multi-Weibull'])
    plot_fit_PcHg(log10(x), y, c1_wb)   
    subplot(211)
    [R2_wb(ind)] = 1 - (norm(multiweibullCDF(log10(x),c1_wb) - y) / norm(y - mean(y)))^2;
    Erro = ['R2 =' num2str(R2_wb(ind))];
    text(0,90,Erro,'FontSize',12,'interpreter','latex');
        
    %*********************************************************************************************************************
    c_wb = [c1_wb;zeros(npar_wb-length(c1_wb),1)]; % Vector of adjusted parameters;
   % Number of parmeters variable acording with the number of Multi-Weibulls
           
   % [pfit_wb(ind,:)] = c_wb'; % Saving a matrix of multiple parameters
    
    choice = questdlg('Approve?','Fit Multi-Weibull','Yes','No','Cancel','Yes');
    if strcmp(choice, 'No')
        ind = ind-1;
    elseif strcmp(choice, 'Cancel')
        error('Cancelled by User')
    end
       
end