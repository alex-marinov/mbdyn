% generate inputs for direct.mbd

N = 1000;
dt = 1e-3;
t = [0:N]'*dt;

% base displacement
x = zeros(size(t));
x(501:end) = 0.01*(1 - cos(5*2*pi*(t(501:end) - t(501))))/2;

save -ascii 'prescribed_displacement.drv' x

% force
f = zeros(size(t));
f(101) = 2000;

save -ascii 'prescribed_force.drv' f

