close all; clear all; clc; warning off;
% This code allows for replication of results reported in Tables 5.1 - 5.3
% for the DRI estimator. It generates empirical null rejection
% probabilities as the share of simulation draws (MC) where the the
% two-sided hypothesis discussed in Section 4 was rejected.

% Complementary functions:
% 1) dr1.m % generates ATE estimates, standard errors (also for bootstrap)
% 2) m_bag.m %bootstrap aggregation of m-choice selection, described in
% Appendix B.3 in detail

% Global parameters:
MC = 6000; % monte-carlo repetitions
Bh = 500; % B_2 bootstrap samples
Bbag = 10000; % B_1 bootstrap samples
n = 5e2; % sample size


la = 50; a = linspace(0,0.5,la); % parameters for Crump(2009)
tc = 0.04;  % parameter for Huber(2013)

% Regularity of the problem
c = 0.5; c0 = 0; mux = 0;   %regular identification
%c = 1.5; c0 = 0; mux = 0;  %irregular identification
%c = 0.9; c0 = 0; mux = 0;  %irregular identification through misspecification

% Note, for irregular identification through misspecification the following
% changes are to be implemeted in the MC loop:
% 1. New DGP:
% x = random('uniform',0,1,n,1);
% bx = c0 + c*x;
% d = (bx-random('uniform',0,1,n,1)>0);
% p = cdf('uniform',bx,0,1);
% 2: x=log(x); to be inserted in line 67

% Potential outcomes:
a1 = 0.01; b1 = 0.25; a2 = 0.02; b2 = 0.05;
fg = @(x) (exp(-x+mux)./((1+exp(-x+mux)).^2)).*log((a1*(1+exp(-c0-c.*x))+b1)./(a2*(1+exp(-c0-c.*x))+b2));
x = random('logistic',mux,1,1e6,1);
del = integral(fg,min(x),max(x)); % true ATE

% space
i0 = zeros(MC,1); ate = i0;
drh = i0; tnh = i0; blh = i0; ahath = i0; bhath = i0;
pv1h = i0; pv2h = i0; pv3h = i0; pv4h = i0; pv5h = i0; pv6h = i0; pv7h = i0;
for m = 1:MC
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DGP
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rng(m) % control seed for replicability
    x = random('logistic',mux,1,n,1);
    bx = c0 + c*x;
    d = (bx-random('logistic',0,1,n,1)>0);
    p = cdf('logistic',bx,0,1);
    
    mx0 = log(a2 + b2*p);
    mx1 = log(a1 + b1*p);
    y0 = mx0 + randn(n,1);
    y1 = mx1 + randn(n,1);
    y = (1-d).*y0 + d.*y1;
    
    ate(m) = mean(mx1-mx0);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Part II: inference for the estimated phat
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % DR I: phat
    b = glmfit(x,d,'binomial','link','logit','constant','on');
    X = [ones(n,1),x];
    ph = 1./(1+exp(-X*b));
    
    f1 = dr1(y,d,ph,1,0,0);
    muhat1=f1.muhat1;muhat0=f1.muhat0;
    Yh = f1.Ym;
    drh(m) = f1.ate;
    tnh(m) = sqrt(n).*(drh(m)-del)./f1.se;
    
    %2. Fixed Trimming p in [0.1 0.9];
    idx = (ph<=0.9 & ph>=0.1);
    yt = y(idx==1); dt = d(idx==1); pt = ph(idx==1);  nt = length(yt);
    f1 = dr1(yt,dt,pt,1,0,0);
    drfh = f1.ate;
    tntfh = sqrt(nt).*(drfh-del)./f1.se;
    
    %3. Trimming Crump(2009)
    xx = ph.*(1-ph);
    I = 2*sum((repmat(xx,1,la)>=repmat(a.*(1-a),n,1))./repmat(xx,1,la))./sum((repmat(xx,1,la)>=repmat(a.*(1-a),n,1)));
    yy = 1./(a.*(1-a))-I;
    id = (yy<=0);
    as = min(a(id));
    idx = (ph<=(1-as) & ph>=as);
    yt = y(idx==1); dt = d(idx==1); pt = ph(idx==1);nt = length(yt);
    f1 = dr1(yt,dt,pt,1,0,0);
    drth = f1.ate;
    tnth = sqrt(nt).*(drth-del)./f1.se;
    
    %4. Trimming Huber(2013)
    i1 = 1./ph<=tc*sum(1./ph); i2 = 1./(1-ph)<=tc*sum(1./(1-ph));
    idx = i1.*i2;
    yt = y(idx==1); dt = d(idx==1); pt = ph(idx==1);  nt = length(yt);
    f1 = dr1(yt,dt,pt,1,0,0);
    drtrh = f1.ate;
    tntrh = sqrt(nt).*(drtrh-del)./f1.se;
    
    % 5.  m out of n bootstrap
    blh(m) = m_bag(y,d,ph,x,Bbag);
    mb = blh(m);
    
    % Bootstrapping with phat
    J = ceil(rand(n*Bh,1)*n);
    xss = zeros(n,Bh); dss = zeros(n,Bh); yss = zeros(n,Bh);
    xss(:) = x(J); dss(:) = d(J); yss(:) = y(J);
    
    J_m = ceil(rand(mb*Bh,1)*n);
    xss_m = zeros(mb,Bh); dss_m = zeros(mb,Bh); yss_m = zeros(mb,Bh);
    xss_m(:) = x(J_m); dss_m(:) = d(J_m); yss_m(:) = y(J_m);
    
    tntfhx = zeros(Bh,1); tnthx = zeros(Bh,1); tntrhx = zeros(Bh,1);
    tb = zeros(Bh,1); tbm = zeros(Bh,1);
    for bb = 1:Bh
        b = glmfit(xss(:,bb),dss(:,bb),'binomial','link','logit','constant','on');
        Xss = [ones(n,1),xss(:,bb)];
        ph = 1./(1+exp(-Xss*b));
        f1h = dr1(yss(:,bb),dss(:,bb),ph,1,0,0);
        
        u = f1h.Ym;
        %1.iid
        tb(bb) = sqrt(n).*( mean(u)-drh(m))./f1h.se;
        %2. Fixed Trimming p in [0.1 0.9];
        idx = (ph<=0.9 & ph>=0.1);
        yt = yss(idx==1,bb); dt = dss(idx==1,bb); pt = ph(idx==1); xt = xss(idx==1,bb); nt = length(yt);
        f1bh = dr1(yt,dt,pt,1,0,0);
        tntfhx(bb) = sqrt(nt)*(f1bh.ate-drfh)./f1bh.se;
        
        %3. Trimming Crump(2009)
        xx = ph.*(1-ph);
        I = 2*sum((repmat(xx,1,la)>=repmat(a.*(1-a),n,1))./repmat(xx,1,la))./sum((repmat(xx,1,la)>=repmat(a.*(1-a),n,1)));
        yy = 1./(a.*(1-a))-I;
        id = (yy<=0);
        as = min(a(id));
        idx = (ph<=(1-as) & ph>=as);
        yt = yss(idx==1,bb); dt = dss(idx==1,bb); pt = ph(idx==1);  xt = xss(idx==1,bb); nt = length(yt);
        f1bh = dr1(yt,dt,pt,1,0,0);
        tnthx(bb) = sqrt(nt)*(f1bh.ate-drth)./f1bh.se;
        
        %4. Trimming Huber(2013)
        i1 = 1./ph<=tc*sum(1./ph); i2 = 1./(1-ph)<=tc*sum(1./(1-ph));
        idx = i1.*i2;
        yt = yss(idx==1,bb); dt = dss(idx==1,bb); pt = ph(idx==1);  xt = xss(idx==1,bb); nt = length(yt);
        f1bh = dr1(yt,dt,pt,1,0,0);
        tntrhx(bb) = sqrt(nt)*(f1bh.ate-drtrh)./f1bh.se;
        
        %5. moon
        b_m = glmfit(xss_m(:,bb),dss_m(:,bb),'binomial','link','logit','constant','on');
        Xm = [ones(mb,1),xss_m(:,bb)];
        ph_m = 1./(1+exp(-Xm*b_m));
        
        f1h = dr1(yss_m(:,bb),dss_m(:,bb),ph_m,1,0,0);
        tbm(bb) = sqrt(mb).*( f1h.ate-drh(m))./f1h.se;
        
    end
    % Calculate p-values
    %1. iid bootstrap
    pv1h(m) = mean(abs(tb)>abs(tnh(m)));
    %2. Fixed Trimming p in [0.1 0.9];
    pv2h(m) = mean(abs(tntfhx)>abs(tntfh));
    %3. Trimming Crump(2009)
    pv3h(m) = mean(abs(tnthx)>abs(tnth));
    %4. Trimming Huber(2013)
    pv4h(m) = mean(abs(tntrhx)>abs(tntrh));
    %5.  moon bootstrap
    pv5h(m) = mean(abs(tbm)>abs(tnh(m)));
end

PVh = [pv1h,pv2h,pv3h,pv4h,pv5h];
% Empirical null rejection probabilities:
%1. iid bootstrap
%2. Fixed Trimming p in [0.1 0.9];
%3. Trimming Crump(2009)
%4. Trimming Huber(2013)
%5. m-out-of-m bootstrap
mean(PVh<=0.05)'




