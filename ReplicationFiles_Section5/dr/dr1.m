function res = dr1(Y,D,P,hat,Mu1_,Mu0_)
% This funciton calculates DRI ATE estimate and the corresponding
% standard error reported in Appendix B.3.
% INPUT:
% Y - observed outcome
% D - observed treatment indicator
% P - Estimated/true propensity score
% hat: 1 if analysis is based on estimated propensity scores
% Mu_1 and Mu_0 -  known mean functions in case for the known propensities

% OUTPUT: collected in an object with
% .Ym - DRI weighted outcome
% .ate - ATE estimate
% .se -  standard error
n = size(Y,1); B = size(Y,2);
options = optimoptions('fsolve','Display','none');
f1 = @(y1,p1,par) [mean((y1-log(par(1) + par(2).*p1)));...
    mean((y1-log(par(1) + par(2).*p1)).*p1)];
x0 = [0.03;0.03]; % starting values for optimizer for GMM on mean funciton
res.ate = zeros(B,1); res.se = zeros(B,1); res.Ym = zeros(n,B);
for b = 1:B 
    y = Y(:,b);
    d = D(:,b);
    if hat==0 % for known propensity scores
        p=P(:,b);
    else
        p=P;
    end
    if hat ==1 % for estimated propensity scores
        % estimating mean functions
        y1 = y(d==1); y0 = y(d==0); p1 = p(d==1); p0 = p(d==0);
        
        del1 = fsolve(@(par)f1(y1,p1,par),x0,options);
        del0 = fsolve(@(par)f1(y0,p0,par),x0,options);
        
        del1 = real(del1); del0 = real(del0);
        mu1 = log(del1(1) + del1(2)*p); mu0 = log(del0(1) + del0(2)*p);
        mu1 = real(mu1); mu0 = real(mu0);
        % doubly robust I
        res.Ym(:,b) = (y-mu1).*d./p - (y-mu0).*(1-d)./(1-p) + (mu1-mu0);
        res.ate(b) = mean(res.Ym(:,b));
        
        res.se(b) = sqrt(mean(((d.*(y-mu1))./p).^2) + mean((((1-d).*(y-mu0))./(1-p)).^2)...
            + mean((mu1-mu0-res.ate(b)).^2));
        res.muhat1(:,b)= mu1;  res.muhat0(:,b)= mu0;
    else
        % for known propensity scores
        mu1_ = Mu1_(:,b); mu0_ = Mu0_(:,b);
        res.Ym(:,b) = (y-mu1_).*d./p - (y-mu0_).*(1-d)./(1-p) + (mu1_-mu0_);
        res.ate(b) = mean(res.Ym(:,b));
        res.se(b) = std(res.Ym(:,b)-res.ate(b));  
    end
end

end

