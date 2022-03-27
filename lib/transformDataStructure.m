function track = transformDataStructure(traj,file_info, id)
% warning off;
%% change matrix to data structure
track.metadata = file_info;
track.x = traj(:,1); % in um
track.y = traj(:,2); % in um
% track.objarea = traj(:,3);
track.objintensity = traj(:,4);
track.objgap = traj(:,5);
track.xvec = traj(:,7);
track.yvec = traj(:,8);
track.speed = traj(:,9);
track.angle = traj(:,10);
track.time = traj(:,11);
track.frozen = traj(:,12);
track.trajID = id;

n = numel(track.x);
frozen = sum(track.frozen);

if (n-frozen)<3
    track =[];
    return;
end
%% calculate mean square displacement and kurtosis.
% msd(x) returns the msd x frames away.
track.msd = zeros(n-frozen,1);
track.msd_n = zeros(n-frozen,1);
track.msd_n(1) = sum(~track.objgap);
track.kurtosis = zeros(n-frozen,1);
% msd_all = [];
% delay_msd = [];
for i = 1:(n-1-frozen)
    x_sq_disp = (track.x(i+1:end-frozen) - track.x(1:end-i-frozen)).^2;
    y_sq_disp = (track.y(i+1:end-frozen) - track.y(1:end-i-frozen)).^2;
    track.msd(i+1) = nanmean(x_sq_disp + y_sq_disp);
    track.msd_n(i+1) = sum(~isnan(x_sq_disp));
    track.kurtosis(i+1) = kurtosis(sqrt(x_sq_disp + y_sq_disp));
%     msd_all = [msd_all; x_sq_disp + y_sq_disp];
%     delay_msd = [delay_msd; (track.time(i+1)-track.time(1))*ones(numel(x_sq_disp),1)];
end

track.alphamsd = nan(n-frozen,1);
track.alphamsd(1:end-1) = diff(log(track.msd))./diff(log(track.time(1:end-frozen)-track.time(1)));
track.msd_delay = track.time(1:end-frozen)-track.time(1);

%% calculate some parameters
track.acceleration = [NaN; track.speed(2:end-frozen) - track.speed(1:(end-1-frozen))];
track.deltaangle = [NaN; track.angle(2:end-frozen) - track.angle(1:(end-1-frozen))];
track.trajtime = track.time(max(1,end-frozen)) - track.time(1);
track.meanspeed = nanmean(track.speed(1:end-frozen));
track.meanacceleration = nanmean(track.acceleration(1:end-1));
track.varspeed = nanvar(track.speed(1:end-frozen));
track.cvspeed = sqrt(track.varspeed) / track.meanspeed;
track.meanangle = abs(nanmean(track.angle(1:end-frozen)));
track.meanobjintensity = nanmean(track.objintensity(1:end-frozen));
track.varobjintensity = nanvar(track.objintensity(1:end-frozen));
track.cvobjintensity = sqrt(track.varobjintensity) / track.meanobjintensity;

[f, xi] = ksdensity(track.speed);
[~,idx] = max(f);
track.modespeed = xi(idx);
track.rand_speed = track.speed - track.modespeed;
track.rel10_speed = log10(track.speed / track.modespeed);

[f, xi] = ksdensity(track.angle);
[~,idx] = max(f);
track.modeangle = xi(idx);
track.rand_angle = track.angle - track.modeangle;

logit_angle = (wrapToPi(track.rand_angle * mean(diff(track.time))) + pi)/(2*pi);
track.logit_angle = log(logit_angle./(1-logit_angle));

%% calculate diffusion coefficient using CVE method
diff_cve = nan(n-frozen,1);
err_cve = nan(n-frozen,1);
n_cve = nan(n-frozen,1);

for j = 1:((n-frozen)-1)
    [diff_cve(j+1), err_cve(j+1), n_cve(j+1)] = MSD_CVE(track.x(1:j:end-frozen),track.y(1:j:end-frozen),track.time(1:j:end-frozen));
end

track.diffcoeff_cve = diff_cve;
track.diffcoeff_cve_err = err_cve;
track.diffcoeff_cve_n = n_cve;
track.diffcoeff_cve_mean = sum(diff_cve(n_cve>0).*n_cve(n_cve>0))/sum(n_cve(n_cve>0));
track.diffcoeff_cve_runtime = 2*track.diffcoeff_cve_mean/(track.meanspeed^2);

if track.diffcoeff_cve_mean <= 0
    track.diffcoeff_cve_mean = NaN;
    track.diffcoeff_cve_runtime = NaN;
end

%% calculate the velocity auto-correlation function
track.vacf=zeros(n-frozen,1);
% vacf_all = [];
% delay_vacf = [];
for i = 1:n-frozen
    vacf = dot([track.xvec(1:(end-i+1-frozen)), ...
        track.yvec(1:(end-i+1-frozen))], ...
        [track.xvec(i:(end-frozen)), ...
        track.yvec(i:(end-frozen))],2);
    track.vacf(i) = nanmean(vacf);
%     vacf_all = [vacf_all; vacf];
%     delay_vacf = [delay_vacf; (track.time(i)-track.time(1))*ones(numel(vacf),1)];
end

% %% calculate the angle auto-correlation function
% track.aacf=zeros(n,1);
% for i = 1:n
%     temp = dot([track.xvec(1:end-i+1),track.yvec(1:end-i+1)],[track.xvec(i:end),track.yvec(i:end)],2)./(sqrt(sum(([track.xvec(1:end-i+1) track.yvec(1:end-i+1)]').^2)).*sqrt(sum(([track.xvec(i:end) track.yvec(i:end)]').^2)))';
%     temp(isnan(temp) | isinf(temp)) = [];
%     if numel(temp) > 0
%         track.aacf(i) = mean(temp);
%     end
% end

% try
% if n > 9
%% fit model to VACF and MSD

f = track.diffcoeff_cve_runtime;
w = 0;

x0 = [f w track.vacf(1)];
lb = [0 0 0];
ub = [Inf Inf Inf];

if numel(track.msd_delay)<4 || sum(~isfinite(vacf_msd_curvefit(x0,track.msd_delay(1:end-1),track.vacf(1:end-1),track.msd(1:end-1),track.msd_n(1:end-1))))>0
    x = [NaN NaN NaN];
else
    options = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off');
    x = lsqnonlin(@vacf_msd_curvefit,x0,lb,ub,options,...
        track.msd_delay(1:end-1), track.vacf(1:end-1), track.msd(1:end-1), track.msd_n(1:end-1));
end
% if track.trajtime>30
%     plot(track.msd_delay, track.vacf,'o')
%     hold on,
%     msd2vacf = diff(diff(track.msd/2))/(mean(diff(track.msd_delay))^2);
%     plot(track.msd_delay(2:end-1), msd2vacf,'+');
%     plot(track.msd_delay,predict_vacf(x,track.msd_delay));
%     hold off
%     ylim([-track.vacf(1) track.vacf(1)]);
%     % set(gca,'XScale','log');
%     pause;
% end

% track.fit_resnorm = resnorm;
track.fit_runtime = x(1);
track.fit_angle = x(2);
track.fit_speed = sqrt(x(3));
track.fit_diffcoeff = x(3)/(2*x(1)*(1/x(1)^2+x(2)^2));
track.fit_diffcoeff_nocurve = x(3)*x(1)/2;

%% create structure for tumble detection
track.tumble_d2 = [];
track.tumble = [];
track.tumblebias = [];
track.revfreq = [];
track.meanrunspeed = [];
track.runtime = [];
track.meanruntime = [];
track.runlength = [];
track.tumbletime = [];
track.tumbleangle = [];
track.runspeed = [];
track.diffcoeff = [];

% close all,
% figure('Position',[100 500 2100 700])
% subplot(1,3,1),
% plot(track.x,track.y);
% axis equal
% subplot(1,3,2),hold on,
% plot(track.msd_delay, track.vacf,'o');
% plot(unique(delay_vacf),predict_vacf(x,unique(delay_vacf)),'r-');
% plot(unique(delay_vacf),predict_vacf(x0,unique(delay_vacf)),'g-');
% axis square
% subplot(1,3,3),hold on,
% plot(track.msd_delay, track.msd./track.msd_delay,'o');
% plot(unique(delay_vacf),predict_msd(x,unique(delay_vacf))./unique(delay_vacf),'r-');
% plot(unique(delay_vacf),predict_msd(x0,unique(delay_vacf))./unique(delay_vacf),'g-');
% axis square
% pause;

% %% calculate the powerspectrum of vacf
% [power,freq,~] = getAvgSpectrum(track.vacf,track.time(2)-track.time(1),1);
% [~,ind] = max(smooth(power));
% track.vacf_powerangle = freq(ind);

% else
%     track.vacf_powerangle = [];
%     track.fit_resnorm = [];
%     track.fit_runtime = [];
%     track.fit_angle = [];
%     track.fit_speed = [];
%     % track.fit_drift = x(3);
%     track.fit_diffcoeff = [];
%     track.fit_diffcoeff_nocurve = [];
% end

% %% fit model to AACF = exp(-x/b)*cos(c*x)
% w = max(2*pi*abs(track.meanangle),2*pi*track.aacf_powerangle);
% N = numel(track.vacf);
%
% s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[1 w],'Display','off','Weights',N:-2:2,'Lower',[0 0],'Upper',[Inf 2*pi]);
% f = fittype('exp(-x/b)*cos(c*x)','independent','x','options',s);
%
% track.aacf_fit = fit(track.time(1:floor(N/2))-track.time(1), track.aacf(1:floor(N/2)), f);
%
% track.aacf_fit_angle = track.vacf_fit.c/(2*pi);
% track.aacf_fit_rotdiff = 1/track.vacf_fit.b;

% figure,
% scatter(track.time(1:floor(N/2))-track.time(1), track.msd(1:floor(N/2)));
% hold on
% plot(track.time(1:floor(N/2))-track.time(1),p(1).*(track.time(1:floor(N/2))-track.time(1))+p(2));
%
% a = track.vacf_fit.a;
% b = 1/track.vacf_fit.b;
% c = track.vacf_fit.c;
%
% subplot(2,2,1)
% scatter(track.time-track.time(1), track.vacf);
% hold on
% plot(track.time-track.time(1),predict_vacf(x,track.time-track.time(1),ydatavacf(1)),'r');
% hold off;
% subplot(2,2,[2 4])
% plot(track.x, track.y);
% axis equal;
% subplot(2,2,3)
% scatter(track.time-track.time(1), track.msd);
% hold on
% plot(track.time-track.time(1),predict_msd(x,track.time-track.time(1),ydatavacf(1)),'r');hold off;
% pause();

% %% fit model to MSD = 2*a*b*(x-b*(1-exp(-x/b)))
% s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[track.meanspeed^2 0.001],'Display','off','Weights',N:-6:6,'Lower',[0 0],'Upper',[Inf Inf]);
% f = fittype('2*a*b*(x-b*(1-exp(-x/b)))','independent','x','options',s);
%
% [track.msd_fit, track.msd_gof] = fit(track.time(1:floor(N/6))-track.time(1), track.msd(1:floor(N/6)), f);
%
% track.msd_fit_speed = sqrt(track.msd_fit.a);
% track.msd_fit_runtime = track.msd_fit.b;
% track.msd_fit_diffcoeff = 2*track.msd_fit.b*track.msd_fit.a/4;

% if track.msd_fit_diffcoeff>1
%     subplot(2,2,1)
%     plot(track.time(1:floor(N/4))-track.time(1), track.msd(1:floor(N/4)));
%     hold on
%     plot(track.msd_fit);
%     title(track.msd_gof.rmse);
%     hold off;
%     subplot(2,2,[2 4])
%     plot(track.x, track.y);
%     axis equal;
%     subplot(2,2,3)
%     plot(track.time(1:floor(N/4))-track.time(1), track.vacf(1:floor(N/4)));
%     hold on
%     plot(track.vacf_fit);
%     title(track.vacf_gof.rmse);
%     hold off;
%     pause();
% end
end
