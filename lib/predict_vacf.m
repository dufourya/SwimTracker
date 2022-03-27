function C = predict_vacf(x,xdata)

    b = x(1);
    c = x(2);
    a = x(3);
    
    C = a*exp(-xdata./b).*cos(c.*xdata);
end
