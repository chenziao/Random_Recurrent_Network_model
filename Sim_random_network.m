clear;
close all;

%% Parameters
dt = 0.05;	% time step
N = 1000;	% number of neurons

Case = 1;
switch Case
    case 0	% g<1
        g = 0.5;    % activation gain phi(g*x)
        ntrials = 1;	% number of trials
        T = 30;	% simulation time
        Tmin = 10;	% delay for calculating autocorrelation
    case 1	% g>1
        g = 2;
        ntrials = 20;
        T = 300;	% simulation time
        Tmin = 100;
end
phi = @(x) tanh(g*x);	% activation function
rng(10);

maxlag = 50;
nlag = ceil(maxlag/dt);
tmin = ceil(Tmin/dt);
C = zeros(N,2*nlag+1);
Cm = zeros(ntrials,2*nlag+1);
for j = 1:ntrials
%% Initialization
J = randn(N,N)/N^0.5;   % connectivity matrix
h0 = 0 + 1.0*randn(N,1);    % initial state

%% Simulation
[h,t] = sim_net( J, phi, h0, T, dt );

%% Autocorrelation
for i = 1:N
    [C(i,:),lags] = xcorr(h(i,tmin:end),nlag,'unbiased');
end
Cm(j,:) = mean(C,1);

end

%% Plot
% Autocorrelation
tlags = lags*dt;
figure;
subplot(211);
plot(tlags,Cm);
xlim([0,tlags(end)]);	ylim([-0.1,0.5])
ylabel('\Delta');
title(num2str(ntrials,'Estimated autocorrelation of %d trials'));
subplot(212);
plot(tlags,mean(Cm,1),'b','LineWidth',2);
xlim([0,tlags(end)]);	ylim([-0.1,0.5])
xlabel('\tau');	ylabel('\Delta');
title(num2str(ntrials,'Averaged autocorrelation over %d trials'));

% Activity
nshow = 1:5;
figure;
plot(t,h(nshow,:));
xlim(t([1,end]));
xlabel('time'); ylabel('h');
title(num2str(numel(nshow),'Activities of %g units'));
