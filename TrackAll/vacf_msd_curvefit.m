function err = vacf_msd_curvefit(x, xdata, yvacf, ymsd)
    
    fitvacf = predict_vacf(x, xdata, yvacf(1));
    fitmsd = predict_msd(x, xdata, yvacf(1));
    
    err_vacf = fitvacf - yvacf;
    err_msd = fitmsd - ymsd;
    
    weights = sqrt(numel(xdata)-1:-1:0);
    
    err = [err_vacf.*weights'./fitvacf(1) err_msd.*weights'./max(fitmsd)];
    
end