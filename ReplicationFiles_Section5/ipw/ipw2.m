function res = ipw2(Y,D,P,X,hat)
% This funciton calculates IPWII ATE estimate and the corresponding
% standard error reported in Appendix B.3.
% INPUT:
% Y - observed outcome
% D - observed treatment indicator
% P - estimated/true propensity score
% X - confounding variables
% hat: 1 if analysis is based on estimated propensity scores

% OUTPUT: collected in an object with
% .Ym - IPWII weighted outcome
% .ate - ATE estimate
% .se -  standard error
n = size(Y,1); B = size(Y,2);

res.ate = zeros(B,1); res.se = zeros(B,1); res.Ym = zeros(n,B);
for b = 1:B
    y = Y(:,b);
    d = D(:,b);
    p = P(:,b);
    x = [ones(n,1),X(:,b)];
    % Inverse Probability Weighting II
    res.Ym(:,b) = (y.*d./p)./mean(d./p) - (y.*(1-d)./(1-p))./mean((1-d)./(1-p));
    res.ate(b) = mean(res.Ym(:,b));
    if hat ==1
        % for estimated propesnity scores, see Appenix B.3:
        dd= d.*(1-p).*x - (1-d).*p.*x;
        S=inv(dd'*dd./n);
        % potential outcome 1
        mu1h = mean((y.*d./p)./mean(d./p));
        A = mean(d./p);
        k1 = (d./p).*(y-mu1h);
        u = (k1- (mean(k1.*dd)*S*dd')');
        D = mean(u.^2);
        v1 = D/(A^2);
        % potential outcome 0
        mu0h = mean((y.*(1-d)./(1-p))./mean((1-d)./(1-p)));
        A = mean((1-d)./(1-p));
        k0 = ((1-d)./(1-p)).*(y-mu0h);
        u = k0- (mean(k0.*dd)*S*dd')';
        D = mean(u.^2);
        v2 = D/(A^2);
        % covariance
        alpha = mean(k1.*dd)*S*mean(k0.*dd)';
        alpha=alpha/(mean(d./p)*mean((1-d)./(1-p)));
        res.se(b)= sqrt(v1+v2-2*alpha);
    else
        % for known propensity score
        mh1=mean((y.*d./p)./mean(d./p));
        mh0=mean((y.*(1-d)./(1-p))./mean((1-d)./(1-p)));
        res.se(b) = sqrt(mean(d./p.^2.*(y-mh1).^2)/mean(d./p)^2 + mean((1-d)./(1-p).^2.*(y-mh0).^2)/mean((1-d)./(1-p))^2);
    end
end

end

