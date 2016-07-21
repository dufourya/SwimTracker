function track = transformDataStructure(traj,file_info)
warning off;

%% change matrix to data structure
track.metadata = file_info;
track.x = traj(:,1); % in um
track.y = traj(:,2); % in um
track.objarea = traj(:,3);
track.objintensity = traj(:,4);
track.objgap = traj(:,5);
track.xvec = traj(:,7);
track.yvec = traj(:,8);
track.speed = traj(:,9);
track.angle = traj(:,10);
track.time = traj(:,11);
track.frozen = traj(:,12);

frozen = sum(track.frozen);
%% calculate mean square displacement and kersosis.
track.msd = zeros(numel(track.x)-frozen,1);
track.kurtosis = zeros(numel(track.x)-frozen,1);
for i = 1:(length(track.x)-1-frozen)
    x_sq_disp = (track.x(i+1:end-frozen) - track.x(1:end-i-frozen)).^2;
    y_sq_disp = (track.y(i+1:end-frozen) - track.y(1:end-i-frozen)).^2;
    track.msd(i+1) = mean(x_sq_disp + y_sq_disp);
    track.kurtosis(i+1) = kurtosis(sqrt(x_sq_disp + y_sq_disp));
end

%% calculate some parameters
track.acceleration = [track.speed(2:end-frozen) - track.speed(1:(end-1-frozen)); 0];
track.deltaangle = [track.angle(2:end-frozen) - track.angle(1:(end-1-frozen)); 0];
track.trajtime = track.time(max(1,end-frozen)) - track.time(1);
track.meanspeed = mean(track.speed(1:end-frozen));
track.meanangle = abs(mean(track.angle(1:end-frozen)));

%% create structure for tumble detection
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
track.tumbleModel = [];
track.ppState = [];
track.diffcoeff = [];


%% calculate the velocity auto-correlation function
track.vacf=zeros(numel(track.x)-frozen,1);
for i = 1:numel(track.x)-frozen
    track.vacf(i) = mean(dot([track.xvec(1:end-i+1-frozen), ...
                              track.yvec(1:end-i+1-frozen)], ...
                             [track.xvec(i:end-frozen), ...
                              track.yvec(i:end-frozen)],2));
end


try
    %% calculate the powerspectrum of vacf
    [power,freq,~] = getAvgSpectrum(track.vacf,track.time(2)-track.time(1),1);
    [~,ind] = max(power(2:end));
    track.vacf_powerangle = freq(ind+1);

    %% fit model to VACF and MSD
    w = max(2*pi*abs(track.meanangle),2*pi*track.vacf_powerangle);
    
    xdata = track.time(1:end-frozen)-track.time(1);
    ydatavacf = track.vacf(1:floor(numel(xdata)/2));
    ydatamsd = track.msd(1:floor(numel(xdata)/2));
    xdata = xdata(1:floor(numel(xdata)/2));
    x0 = [5 w];
    
    options = optimset('MaxFunEvals',10000,'MaxIter',10000,'Display','off');
    
    [x,resnorm] = lsqnonlin(@vacf_msd_curvefit,x0,[0 0],[Inf Inf],options,xdata,ydatavacf,ydatamsd);
    
    track.fit_resnorm = resnorm/floor(numel(xdata)/2);
    track.fit_runtime = x(1);
    track.fit_angle = x(2)/(2*pi);
    track.fit_diffcoeff = track.meanspeed^2/(2*x(1)*(1/x(1)^2+x(2)^2));
    track.fit_diffcoeff_nocurve = track.meanspeed^2*x(1)/2;
catch exception
    track.vacf_powerangle = [];
    track.fit_resnorm = [];
    track.fit_runtime = [];
    track.fit_angle = [];
    track.fit_diffcoeff = [];
    track.fit_diffcoeff_nocurve = [];
end

end
