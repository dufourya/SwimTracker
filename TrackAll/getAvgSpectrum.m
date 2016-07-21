function [ power, freq, segLength ] = getAvgSpectrum( series, dt, segments )
%GETAVGSPECTRUM Takes a time series and generates the normalized average 
%   power spectrum based on the given dt.
%   
%   [power, freq, segLength] = getAvgSpectrum(series, dt, segments) - Allows you to
%   specify the exact number of segments you wish to split the series into.
%   Series is your series of data.  dt is the time in seconds between each 
%   event in the series.
%
%   5/03/2007
%   Michael Sneddon

%Split up the series into equal parts and do the math
series = series(1:(end-mod(length(series),segments)));
stateCount = length(series)/segments;

for i = 0:(segments-1)
    
    startIndex = round(i*stateCount)+1;
    endIndex = round((i+1)*stateCount);
    
    transformedSeries = fft(series(startIndex:endIndex));
    N = length(transformedSeries);
    transformedSeries(1) = [];
    segpower = abs(transformedSeries(1:round(N/2))).^2;
    if(size(segpower,2) > 1), segpower = segpower'; end;
    
    if(i==0) 
        power = segpower; 
    else
        power = segpower + power;
    end;
    
    
end;

segLength = (endIndex - startIndex + 1)*dt;
power = power./segments;

%Finish up and return the results
freq = (1:round(N/2))/(round(N/2))*(1/dt)*1/2;   %added *1/2
power = power*dt/stateCount;     