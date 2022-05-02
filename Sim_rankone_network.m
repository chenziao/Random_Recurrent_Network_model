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
figure;
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

hold on;
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

% random connection
chi = g*randn(N,N)/N^0.5;
J = chi + P;	% connectivity matrix
h0 = 0 + 1.0*randn(N,1);	% initial state

%% Simulation
[h,t] = sim_net( J, phi, h0, T, dt );

%% Plot
% Activity
nshow = 1:10;
figure;
plot(t,h(nshow,:));
xlim(t([1,end]));
xlabel('time'); ylabel('h');
title(num2str(numel(nshow),'Activities of %g units'));

figure;
plot(mn(:,1),h(:,end),'.');
xlabel('m');    ylabel('h');




