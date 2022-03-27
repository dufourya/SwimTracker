function err = vacf_msd_curvefit(x, delay, vacf_all, msd_all, n)
    
    fitvacf = predict_vacf(x, delay);
    
    err_vacf = (fitvacf - vacf_all).*n;
    
    msd2vacf = diff(diff(msd_all/2))/(mean(diff(delay))^2);
    err_msd = (fitvacf(2:end-1) - msd2vacf).*n(2:end-1);
    
    err = [err_vacf' err_msd'];
end