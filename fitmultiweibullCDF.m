function c=fitmultiweibullCDF(x,y,c0,maxit,tol)

if nargin<4 % Number of input arguments
    tol=1e-6;
end

if nargin<4
    maxit=1000;
end

a=100;
c=c0;
f=@(x,c) multiweibullCDF(x,c);

for k=1:maxit
    F = f(x, c) - y; % Inicial Error   
    e1=norm(F); % norm of Inicial Error
    Jac = Jacobiana(x, f, c); % Partial derivatives matrix (central limit)    
    dc = -(Jac'*Jac + a*a*eye(size(Jac,2)))\(Jac'*F); % Fit = non-linear Least Squares: "Levemberg-Marquadt" --> "Gauss-Newton"
    e2=norm(f(x,c+dc) - y); % norm of the updated Error

    if norm(dc)/norm(c)<tol % Stop criteria
        break
    end
    
    if e2<e1 && isreal(dc)
        c=c+dc; % Step  forward solution
        a=0.9*a; % Reduced step
        
        for m=1:length(c)/4 % Parameters
                        
            if c(4*m-3)<0 % Constraint of Amplitude <0%
                c(4*m-3)=c(4*m-3)-dc(4*m-3);
                a=1.1*a;
            end
            if c(4*m-3)>100 % Constraint of Amplitude >100%
                c(4*m-3)=c(4*m-3)-dc(4*m-3);
                a=1.1*a;
            end
            
%             if c(4*m-2)<-1 % Constraint of Entry pressure <0.1 psi
%                 c(4*m-2)=c(4*m-2)-dc(4*m-2);
%                 a=1.1*a;
%             end
%             
%              if c(4*m-2)>1000 % Constraint of Entry pressure >0.5 psi
%                 c(4*m-2)=c(4*m-2)-dc(4*m-2);
%                 a=1.1*a;
%             end
            
            if 10^c(4*m-1)<=0 % Constraint of scale parameter (tortuosity)
                c(4*m-1)=abs(c(4*m-1));
                a=1.1*a;
            end
            
            if 10^c(4*m-1)>100 % Constraint of scale parameter (tortuosity)
                c(4*m-1)=abs(c(4*m-1));
                a=1.1*a;
            end
            
            if c(4*m)<2 || c(4*m)>20% Constraint of obliquity parameter (Beta)
                c(4*m)=c(4*m)-dc(4*m);
                a=1.1*a;
            end
        end
        
        if sum(c(1:4:end))>100 % Constraint of amplitude >100% or <0
            c(1:4:end) = c(1:4:end)/sum(c(1:4:end))*100;
        end
        
    else
        a=1.1*a;  % Aumenta o passo        
    end
    
    fprintf('k=%d, alpha=%f, corr=%e, erro=%e\n',k,a,norm(dc)/norm(c),e2)
    
end

% Sort the Weibull functions by increasing size of modal pore throat 
cl = reshape(c,[],length(c)/4);
xmoda = cl(2,:) + cl(3,:).*(((cl(4,:)-1)./cl(4,:)).^(1./cl(4,:)));
[foo, I] = sort(xmoda,'ascend');
c=reshape(cl(:,I),[],1);



