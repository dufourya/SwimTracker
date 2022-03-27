function normalizedVals = imquantnorm(values)

tiedFlag = true;
avgFcn = @nanmean;

% allocate some space for the normalized values
normalizedVals = values;
valSize = size(values);
rankedVals = NaN(valSize);

% find nans
nanvals = isnan(values);
numNans = sum(nanvals);
ndx = ones(valSize);
N = valSize(1);

% create space for output
if tiedFlag
    rr = cell(size(values,2),1);
end

% for each column we want to ordered values and the ranks with ties
for col = 1:valSize(2)
    [sortedVals,ndx(:,col)] = sort(values(:,col));
    if(tiedFlag)
        rr{col} = sort(tiedrank(values(~nanvals(:,col),col)));
    end
    M = N-numNans(col);
    % interpolate over the non-NaN data to get ranked data for all points
    rankedVals(:,col) = interp1(1:(N-1)/(M-1):N,sortedVals(1:M),1:N);
end

% take the mean of the ranked values
mean_vals = feval(avgFcn,rankedVals,2);

% Extract the values from the normalized distribution
for col = 1:size(values,2)
    M = N-numNans(col);
    if tiedFlag
        normalizedVals(ndx(1:M,col),col) = interp1(1:N,mean_vals,1+((N-1)*(rr{col}-1)/(M-1)));
    else
        normalizedVals(ndx(1:M,col),col) = interp1(1:N,mean_vals,1:(N-1)/(M-1):N);
    end
end

