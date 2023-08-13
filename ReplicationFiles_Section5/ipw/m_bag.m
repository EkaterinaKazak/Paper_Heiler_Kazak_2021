function [ mboot2] = m_bag(y_,d_,p_,x,B)
% bootstrap aggregation of m-choice selection.
% INPUT: 
% y_ - observed outcome
% d_ - observed treatment indicator
% p_ - estimated/true propensity score
% x  - confounding variables (for the score)
% B  - B_1 number of bootstrap iterations
% INPUT: 
% mboot2 - optimal choice of m for the m-out-of-n bootstrap

q = 0.75;  K=100; % = S in the paper
n = length(y_);
mgrid = ceil(n.*q.^(0:1:40));mgrid = mgrid(mgrid>(log(n))^2);mgrid(diff(mgrid)==0)=[];
lm = length(mgrid);

mopt2=zeros(K,1);
for k=1:K
    tb2 = zeros(B,lm);
    for j = 1:lm
        m = mgrid(j);
        % m out of n bootstrap (sampling with replacement)
        J = ceil(rand(m*B,1)*n);
        xss = zeros(m,B); dss = zeros(m,B); 
        xss(:) = x(J); dss(:) = d_(J);
        pss = zeros(m,B);pss(:) = p_(J);
        sed=zeros(B,1); ated=zeros(B,1);
        for b = 1:B
            d = dss(:,b); x_ = [ones(m,1),xss(:,b)];
            p = pss(:,b); %
            % compute self-normalization
            dd= d.*(1-p).*x_ - (1-d).*p.*x_;
            S=inv(dd'*dd./m);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % for the D/p ratio
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % D/p
            Dm = (d-p)./(p.*(1-p));
            ated(b) = mean(Dm);
            % potential outcome 1
            k1 = ((d-p)./p);
            u1 = (k1- (mean(k1.*dd)*S*dd')');
            v1 = mean(u1.^2);
            % potential outcome 0
            k0 = ((p-d)./(1-p));
            u0 = k0- (mean(k0.*dd)*S*dd')';
            v2 = mean(u0.^2);
            % covariance
            alpha = mean(u1.*u0);
            % standard error d/p (see derivations in Appendix B.3):
            sed(b) = sqrt(v1+v2-2*alpha);
        end
    tb2(:,j)=(sqrt(m).*(ated)./(sed))';
    
    end
RHO2 = zeros(lm-1,1);
for i =1:(lm-1)
    [~,~,RHO2(i)] = kstest2(tb2(:,i+1),tb2(:,i));
end

[~,iopt2] = min(RHO2);
mopt2(k) = mgrid(iopt2);
end
mboot2= floor(mean(mopt2));
end