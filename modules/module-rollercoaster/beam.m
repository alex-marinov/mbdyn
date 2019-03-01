% $Header$

L = 1.;
EJ = 1.e3;

N = 20;
d = [0:N]'/N*L;

n = 10;
for i = 1:n,
	GA = 1.e3*i;

	% esatto
	vexact(:,i) = d./GA + 1./6./EJ*(d.^2).*(3.*L-d);

	% FEM 2 nodi rilassato
	vfem(:,i) = d./GA + 1./4./EJ*L^2*d;

	% FV 2 nodi
	vfv(:,i) = (L/GA+1./4.*L^2/EJ*(2.*d-L)).*[zeros(N/2,1);ones(N/2+1,1)];
end

