function [h,t] = sim_net( J, phi, h0, T, dt )
% Simulate recurrent neural network
% J - connection matrix
% phi - activation function
% h0 - initial state
% T - simulation time
% dt - time step size
h0 = h0(:);
nt = floor(T/dt)+1;
h = zeros(size(h0,1),nt);

h(:,1) = h0;
prog = 0;   quant = 100;
fmt = '%3.0f%%';
nb = fprintf(fmt,0);
bksp = repmat('\b',[1,nb]);
for t = 1:nt-1
    h(:,t+1) = h(:,t)+(-h(:,t)+J*phi(h(:,t)))*dt;
    if floor(quant*t/nt)>prog
        prog = floor(quant*t/nt);
        fprintf([bksp,fmt],100*prog/quant);
    end
end
fprintf([bksp,fmt,'\n'],100);
t = (0:nt-1)*dt;
end

