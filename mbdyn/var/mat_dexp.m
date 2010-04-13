function Gamma = mat_dexp(theta)

t = sqrt(theta'*theta);

if (t > 0.),
	a = sin(t)/t;
	b = (1 - cos(t))/t^2;
	c = (1 - a)/t^2;

	theta_cross = cross(theta);

	Gamma = eye(3) + b*theta_cross + c*theta_cross*theta_cross;

else
	Gamma = eye(3);
end

