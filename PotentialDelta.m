% clear;
% close all;

%% Functions
% Integral of phi
Phi_g = @(g) @(x) log(cosh(g*x))/g;

% Gauss-Hermite method
[~,r,w] = hermipol(65);

% Potential function
V_D = @(D0,Phi) @(Ds) -Ds.^2/2+arrayfun(@(D) intGauss(@(zs) arrayfun(@(z) intGauss(@(x) ...
    Phi((D0-abs(D))^0.5*x+abs(D)^0.5*z),r,w).^2,zs),r,w),Ds);

%% Evaluate
Case = 1;
nd = 20;
switch Case
    case 0	% g<1
        g = 0.5;
        nd0 = 3;
        D0 = linspace(0.5,1.5,nd0); % Delta(0)
    case 1	% g>1
        g = 2;
        nd0 = 6;
        D0 = linspace(0.44,0.485,6);	% Delta(0)
end

d = linspace(-1,1,2*nd+1);
v = zeros(nd0,2*nd+1);
for i = 1:nd0
    d0 = D0(i);
    V = V_D(d0,Phi_g(g));
    v(i,:) = V(d0*d);
    v(i,:) = v(i,:)-v(i,nd+1);
end

cmap = jet(nd0);
figure; hold on;
for i = 1:nd0
    plot(d,v(i,:),'color',cmap(i,:),'LineWidth',2);
end

switch Case
    case 0	% g<1
        xlim(D0(ceil(nd0/2))*[-1,1]);
        ylim([min(v(ceil(nd0/2),1)),0]);
    case 1	% g>1
        axis tight;
        yl = get(gca,'YLim');   ylim(yl(1)*[1,-1]);
end
xlabel('\Delta/\Delta(0)');
title('V(\Delta)-V(0)');
leg = legend(arrayfun(@(x) num2str(x,'%.3f'),D0,'UniformOutput',0),'Location','South');
set(leg,'EdgeColor',[1,1,1]);
set(get(leg,'Title'),'String','\Delta(0)');

%% Solution
if Case==1
    err = @(D0) arrayfun(@(d0) diff(feval(V_D(d0,Phi_g(g)),[0,d0])),D0);
    D00 = fminsearch(@(d0) err(d0).^2,D0(ceil(nd0/2)));
    V = V_D(D00,Phi_g(g));
    
    phi = @(x) tanh(g*x);
    dV = @(Ds) -Ds+arrayfun(@(D) intGauss(@(zs) arrayfun(@(z) intGauss(@(x) ...
        phi((D00-abs(D))^0.5*x+abs(D)^0.5*z),r,w).^2,zs),r,w),Ds);
    
    dt = 0.005;
    maxlag = 20;
    nlag = ceil(maxlag/dt);
    C = zeros(1,nlag+1);
    C1 = C;
    C(1) = D00;
    for t = 1:nlag
        C(t+1) = C(t) + C1(t)*dt;
        C1(t+1) = C1(t) - dV(C(t))*dt;
        if C(t+1)<=0
            nlag = t-1;
            C(t+1:end) = [];
            C1(t+1:end) = [];
            break;
        end
    end
    ttlags = (0:nlag)*dt;
    
    figure;
    subplot(211);
    plot(D00*d,V(D00*d),'LineWidth',2);
    axis tight;
    xlabel('\Delta');
    title('V(\Delta)-V(0)');
    subplot(212);
    plot(ttlags,C,'LineWidth',2);
    xlim(ttlags([1,end]));
    ylim([0,D00]);
    xlabel('\tau');
    title('\Delta(\tau)');
    
    figure(1);	subplot(212);   hold on;    % simulated results
    tlagss = tlags(tlags>ttlags(end));
    plot([ttlags,tlagss],[C,zeros(size(tlagss))],'r','LineWidth',2);
    legend({'Estimated','Predicted'});
end

