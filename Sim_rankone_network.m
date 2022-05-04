clear;
close all;

%% Parameters
N = 1000;	% number of neurons
g = 0.9;	% random strength
Mm = 1;     % structure mean m
Mn = 1.2;	% structure mean n
Sig_m = 0.4;	% sigma m
Sig_n = 0.8;	% sigma n
rho = 0.25;     % correlation coefficient m, n
phi = @(x) tanh(x);	% activation function

dt = 0.02;	% time step
T = 30;	% simulation time
ntrials = 1;

ntrials = 100;
dt = 0.1;
T = 20;
rng(0);

%% Theoretical stationary solution
dphi = @(x) sech(x).^2; % derivative of phi
phi2 = @(x) tanh(x).^2; % square of phi
[~,r,w] = hermipol(65); % Gauss-Hermite method
IntG = @(f,mu,d0) intGauss(@(x) f(mu+d0^0.5*x),r,w);
% update rules
D0 = @(mu,d0) g^2*IntG(phi2,mu,d0) + (Sig_m/Mm*mu)^2;
MU = @(mu,d0) Mm*Mn*IntG(phi,mu,d0) + rho*Sig_m*Sig_n*mu*IntG(dphi,mu,d0);

% phase plane
ng = 15;	% number of grid points
bm = 0.8; bd = 0.5;	% bound of mu and d0
[gm,gd] = meshgrid(linspace(-bm,bm,ng),linspace(0,bd,ng));
vm = zeros(size(gm));
vd = zeros(size(gd));
for i = 1:ng^2
    vm(i) = MU(gm(i),gd(i))-gm(i);
    vd(i) = D0(gm(i),gd(i))-gd(i);
end
figure(1);
quiver(gm,gd,vm,vd);
axis tight;
xlabel('\mu');	ylabel('\Delta_0');

iter = 1000;	d = 0.2;
mu = zeros(2,iter);
d0 = zeros(2,iter);
d0(:,1) = bd*rand(2,1);
mu(:,1) = bm*[-1;1].*rand(2,1);
for i = 1:2
    for t = 1:iter-1
        mu(i,t+1) = mu(i,t)+(MU(mu(i,t),d0(i,t))-mu(i,t))*d;
        d0(i,t+1) = d0(i,t)+(D0(mu(i,t),d0(i,t))-d0(i,t))*d;
    end
end

figure(1);  hold on;
for i = 1:2
    plot(mu(i,:),d0(i,:),'r');
    plot(mu(i,end),d0(i,end),'ro');
end
xlabel('\mu');  ylabel('\Delta_0');
% figure;
% for i = 1:2
%     subplot(2,1,i); hold on;
%     plot(mu(i,:),'b');	plot(d0(i,:),'r');
% end
% xlabel('iterations');	legend({'\mu','\Delta_0'});

%% Initialization
% structured connection
Sigma = [Sig_m^2,rho*Sig_m*Sig_n;rho*Sig_m*Sig_n,Sig_n^2];
mn = mvnrnd([Mm,Mn],Sigma,N);
P = mn(:,1)*mn(:,2)'/N;

figure;
plot(mn(:,1),mn(:,2),'.');
xlabel('m');    ylabel('n');


hs = zeros(N,ntrials);
for j = 1:ntrials
% random connection
chi = g*randn(N,N)/N^0.5;
J = chi + P;	% connectivity matrix
h0 = 0 + 1.0*randn(N,1);	% initial state

%% Simulation
[h,t] = sim_net( J, phi, h0, T, dt );

hs(:,j) = h(:,end);
end

%% Calculate sample first/second moments
mu_t = mean(hs,1);
sgn = sign(mu_t);
hs = bsxfun(@times,hs,sgn);
mu_i = mean(hs,2);
d0_i = var(hs,[],2);
mu_e = mean(mu_i);
d0_e = var(hs(:));
kappa = mean(mn(:,2).*phi(mu_i));

figure(1);
plot(mu_e*[-1,1],d0_e*[1,1],'go');

c_mu = corrcoef(mn(:,1),mu_i);
c_mu = c_mu(2);
disp(c_mu);
figure; hold on;
plot(mu(2,end)/Mm*mn(:,1),mu_i,'b.');
axis tight;
xl = get(gca,'XLim');
plot(xl,xl,'r');
xlabel('m_i\kappa');	ylabel('\mu _i');

%% Plot
% Activity
nshow = 1:10;
figure;
plot(t,h(nshow,:));
xlim(t([1,end]));
xlabel('time'); ylabel('h');
title(num2str(numel(nshow),'Activities of %g units'));

