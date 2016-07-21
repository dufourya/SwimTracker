function F = predict_vacf(x,xdata,v2)
    a = v2;
    b = x(1);
    c = x(2);
    %a = x(3);
    F = a*exp(-xdata./b).*cos(c.*xdata);
    
    %d = x(3);
    %F = (a-d)*exp(-xdata./b).*cos(c.*xdata)+d;
end
