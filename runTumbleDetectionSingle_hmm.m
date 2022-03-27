function [tracks, tumble_model] = runTumbleDetectionSingle_hmm(tracks, params, learn_model_hmm)

%     startMatlabPool();
    if nargin < 3 || learn_model_hmm <0
       learn_model_hmm = 0;
    end
    
%     learn_model_hmm  = params.runTumbleDetectionSingle.learn_model_hmm;
    
%     if ~isfield(tracks,'ppState')
%         tracks(1).ppState = [];
%     end
    
    if ~isfield(tracks, 'diffcoeff')
        tracks(1).diffcoeff = [];
    end
    
    if learn_model_hmm
        fprintf('Learning tumble model...');
        %%
        X = cell(numel(tracks),1);
        for i = 1:numel(tracks)
%             X{i} = [tracks(i).acceleration(1:end-2)./tracks(i).speed(2:end-1) tracks(i).deltaangle(1:end-2)]';
              X{i} = (tracks(i).acceleration(2:end-2).*tracks(i).deltaangle(2:end-2)./(tracks(i).speed(2:end-2).^2))';
        end
%         x= cell2mat(X');
%         figure,
%         hist(x(1,:),100);
%         figure,
%         plot(x(1,:),x(2,:),'.')
        %%
        %'gauss', 'mixGaussTied', 'discrete', 'student'
        X = X(cellfun(@numel,X)>10);
        X = X(cellfun(@(x) sum(isnan(x))==0,X));
        tumble_model = hmmFit(X,learn_model_hmm,'gauss');
        fprintf('done.\n');
    else
        tumble_model = params.tumble_model;
    end

    fprintf('Detecting tumbles...');
    TUMBLE = find(squeeze(tumble_model.emission.Sigma)>0.1);
    for i = 1:numel(tracks)
        if numel(tracks(i).x)<5
            tracks(i).tumble = tracks(i).x * NaN;
%             tracks(i).tumbleModel = tumble_model;
            tracks(i).meanrunspeed = NaN;
%             tracks(i).ppState = nan(numel(tracks(i).acceleration),3);    
        else
%             X = [tracks(i).acceleration(1:end-2)./tracks(i).speed(2:end-1) tracks(i).deltaangle(1:end-2)]';
              X = (tracks(i).acceleration(2:end-1).*tracks(i).deltaangle(2:end-1)./(tracks(i).speed(2:(find([tracks(i).frozen; 1],1)-2))).^2)';
            tracks(i).tumble = [NaN ismember(hmmMap(tumble_model, X),TUMBLE) NaN]';
%             tracks(i).tumbleModel = tumble_model;
            tracks(i).meanrunspeed = mean(tracks(i).speed(tracks(i).tumble==0));
%             tracks(i).ppState = nan(numel(tracks(i).acceleration),3);    
        end
    end
    fprintf('done.\n');
        
    fprintf('Calculating track statistics...');
    for j = 1:numel(tracks)
        times = [tracks(j).time(2:end)-tracks(j).time(1:(end-1)); ...
                 tracks(j).time(2)-tracks(j).time(1)];
        
        tracks(j).tumblebias = sum(tracks(j).tumble==1) / ...
                                   sum(~isnan(tracks(j).tumble));
        
        tracks(j).tumbletime = [];
        tracks(j).runtime = [];
        tracks(j).runlength = [];
        tracks(j).tumbleangle = [];
        
        tumbleind = find(tracks(j).tumble==1);
        tumbletime = 0;
        
        if length(tumbleind)>1
            k=1;
            for i=2:length(tumbleind)
                tumbletime = tumbletime + (tracks(j).time(tumbleind(i-1)) - ...
                    tracks(j).time(tumbleind(i-1)-1));
                if (tumbleind(i)-tumbleind(i-1))>1
                    tracks(j).tumbletime = [tracks(j).tumbletime tumbletime];
                    tumbletime = 0;
                    runtime = tracks(j).time(tumbleind(i)) - ...
                        tracks(j).time(tumbleind(i-1)+1);
                    tracks(j).runtime = [tracks(j).runtime runtime];
                    runlength = ...
                        sum(tracks(j).speed(tumbleind(i-1):tumbleind(i)-1) .* ...
                        times(tumbleind(i-1):tumbleind(i)-1));
                    tracks(j).runlength = [tracks(j).runlength, runlength];
                    tracks(j).tumbleangle = [tracks(j).tumbleangle, ...
                        180/pi*atan2(tracks(j).xvec(tumbleind(k)-1)* ...
                        tracks(j).yvec(tumbleind(i-1)+1)- ...
                        tracks(j).yvec(tumbleind(k)-1)* ...
                        tracks(j).xvec(tumbleind(i-1)+1), ...
                        tracks(j).xvec(tumbleind(k)-1)* ...
                        tracks(j).xvec(tumbleind(i-1)+1)+ ...
                        tracks(j).yvec(tumbleind(k)-1)* ...
                        tracks(j).yvec(tumbleind(i-1)+1))];
                    k=i;
                end
            end
            tracks(j).tumbletime = [tracks(j).tumbletime tumbletime];
        elseif length(tumbleind)==1
            tracks(j).tumbletime = (tracks(j).time(tumbleind(1)+1)- ...
                tracks(j).time(tumbleind(1)));
            tracks(j).tumbleangle = [tracks(j).tumbleangle, ...
                180/pi*atan2(tracks(j).xvec(tumbleind(1)-1)* ...
                tracks(j).yvec(tumbleind(1)+1)- ...
                tracks(j).yvec(tumbleind(1)-1)* ...
                tracks(j).xvec(tumbleind(1)+1), ...
                tracks(j).xvec(tumbleind(1)-1)* ...
                tracks(j).xvec(tumbleind(1)+1)+ ...
                tracks(j).yvec(tumbleind(1)-1)* ...
                tracks(j).yvec(tumbleind(1)+1))];
        end
        
        tracks(j).runspeed = tracks(j).runlength./tracks(j).runtime;
        tracks(j).meanruntime = mean(tracks(j).runtime);
        tracks(j).revfreq = 2*(length(tracks(j).tumbletime))/(tracks(j).time(end-1)-tracks(j).time(2));
        tracks(j).diffcoeff = 0.5.*tracks(j).meanruntime.*tracks(j).meanspeed.^2;
    end
    fprintf('done.\n');
end
