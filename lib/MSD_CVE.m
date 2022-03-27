function [D, err_var, n] = MSD_CVE(x,y,t)

dx = diff(x);
dy = diff(y);
d2x = dx(1:end-1).*dx(2:end);
d2y = dy(1:end-1).*dy(2:end);
dt = diff(t);

dx = dx(isfinite(dx));
dy = dy(isfinite(dy));
d2x = d2x(isfinite(d2x));
d2y = d2y(isfinite(d2y));
dt = dt(isfinite(dt));

dt = mean(dt);
R = 1/6;

Dx = mean(dx.^2)/(2*dt) + mean(d2x)/dt;
Dy = mean(dy.^2)/(2*dt) + mean(d2y)/dt;

D = mean([Dx Dy]);

err_varx = R * mean(dx.^2) + (2*R-1) * mean(d2x);
err_vary = R * mean(dy.^2) + (2*R-1) * mean(d2y);

err_var = mean([err_varx err_vary]);

n = sum(~isnan(d2x));

end