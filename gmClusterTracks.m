function clusters = gmClusterTracks(tracks,maxComponents,minComponents)

if nargin <3 || isempty(minComponents) || minComponents < 1
    minComponents = 1;
end
warning('off');
meanspeed = log10([tracks.meanspeed]);
meanaccel = [tracks.meanacceleration]./[tracks.meanspeed];
varspeed  = sqrt([tracks.varspeed])./[tracks.meanspeed];
% meaninten = real(log10([tracks.meanobjintensity]));
% varintens = sqrt([tracks.varobjintensity])./[tracks.meanobjintensity];
diffcoeff = log10([tracks.diffcoeff_cve_mean]);
meanrunti = log10([tracks.diffcoeff_cve_runtime]);
trajtime  = [tracks.trajtime];

options = statset('MaxIter',1000);
BIC = ones(1,maxComponents)*Inf;
GMModels = cell(1,maxComponents);

% X = vertcat(meanspeed,meanaccel,varspeed,meaninten,varintens,diffcoeff,...
%     meanrunti);
if isrow(meanspeed)
    X = vertcat(meanspeed,meanaccel,varspeed,diffcoeff,meanrunti)';
else
    X = horzcat(meanspeed,meanaccel,varspeed,diffcoeff,meanrunti);
end

if sum(~isnan(diffcoeff)) > 5
    k = minComponents;
    nBIC = 0;
    
    while nBIC<1 && k<=maxComponents
        fprintf('Clustering with %d components...', k);
        GMModels{k} = fitgmdist(X,k,'Options',options,'CovarianceType','full',...
            'RegularizationValue',10^-10,'Replicates',5);
        
        if GMModels{k}.BIC >= min(BIC)
            nBIC = nBIC + 1;
        else
            nBIC = 0;
        end
        
        BIC(k)= GMModels{k}.BIC;
        fprintf('BIC = %d\n',BIC(k));
        
        k = k+1;
    end
    
    [minBIC, numComponents] = min(BIC);
    BestModel = GMModels{numComponents};
    
    %%
    % unsortedclusters = ones(size(X',1),1);
    % posteriors = BestModel.posterior(X');
    % r = rand(numel(unsortedclusters),1);
    % unsortedclusters = sum(r>cumsum(posteriors,2),2)+1;
    % unsortedclusters(isnan(posteriors(:,1))) = 0;
    %%
    unsortedclusters = BestModel.cluster(X);
    unsortedclusters(isnan(unsortedclusters))=0;
    %%
    
    cumulTime = zeros(1,numComponents);
    
    for k = 1:numComponents
        cumulTime(k) = sum(trajtime(unsortedclusters==k));
    end
    
    [~, order] = sort(cumulTime,'descend');
    clusters = unsortedclusters;
    for k = 1:numComponents
        clusters(unsortedclusters==order(k)) = k;
    end
    fprintf('Final clusters with %d components (BIC = %d)\n', numComponents, minBIC);
else
    clusters = ones(numel(tracks),1);
    clusters(isnan(diffcoeff)) = 0;
    fprintf('Not enough tracks to cluster!\n');
end

warning('on');