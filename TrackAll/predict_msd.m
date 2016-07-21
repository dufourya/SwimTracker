function F = predict_msd(x,xdata,v2)


b = x(1);
c = x(2);
a = v2;

t = xdata;
%a = x(3);

F = (2*a*b*(b.*cos(c.*t) - b.*exp(t./b) + t.*exp(t./b) - 2*b^2*c.*sin(c.*t) - b^3*c^2.*cos(c.*t) + b^3*c^2.*exp(t./b) + b^2*c^2.*t.*exp(t./b)))./(exp(t./b)*(b^2*c^2 + 1)^2);

%d = x(3);
%F =(2*(a-d)*b^4*c^2 + 2*(a-d)*b^3*c^2*t - 2*(a-d)*b^2 + 2*(a-d)*b*t + d*b^4*c^4*t.^2 + 2*d*b^2*c^2*t.^2 + d*t.^2)./(b^4*c^4 + 2*b^2*c^2 + 1) - (2*(a-d)*b^4*c^2*cos(c*t) - 2*(a-d)*b^2*cos(c*t) + 4*(a-d)*b^3*c*sin(c*t))./(exp(t/b) + 2*b^2*c^2*exp(t/b) + b^4*c^4*exp(t/b));
end