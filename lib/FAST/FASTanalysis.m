function FASTanalysis(file, varargin)

fprintf('FAST data analysis...');
visfig = 'On';
screensize = get( groot, 'Screensize' );

if contains(varargin,'hpcc')
    visfig = 'Off';
end

%% load data file
fileInfo = RadFileInfo(file);

load(strcat(fileInfo.name, RadFileInfo.swimtrackerExt),'tracks');
metadataTracks = Metadata(fileInfo.meta_file);
metadataTracks.read();
metadataPics = Metadata(strcat(fileInfo.name, FileInfo.picsMetaExt));
metadataPics.read();
load(strcat(fileInfo.name, FileInfo.picturesExt),'pictures');

pixelSizeTracks = metadataTracks.getPixelSize();
pixelSizePics = metadataPics.getPixelSize();

imageHeight = metadataTracks.getImageHeightPixel();
imageWidth = metadataTracks.getImageWidthPixel();

% xcorrection = metadataPics.getXCorrection();
% ycorrection = metadataPics.getYCorrection();

%% find trajectories that finish at the freezing time
end_times = cellfun(@max,{tracks.time},'UniformOutput', false);
t = [end_times{:}] == max([end_times{:}]);
tracks0 = tracks(~t);
tracks = tracks(t);
traj_coordinates = zeros(numel(tracks),2);

for i = 1:numel(tracks)
    traj_coordinates(i,:) = [tracks(i).x(end), tracks(i).y(end)];
end

coordinates = cat(2,pictures.XYcoord);
coordinates = coordinates';

coordinates(:,1) = (imageWidth/2*pixelSizeTracks) + coordinates(:,1);
coordinates(:,2) = (imageHeight/2*pixelSizeTracks) + coordinates(:,2);

correction = (max(traj_coordinates)+min(traj_coordinates))/2 - (max(coordinates)+min(coordinates))/2;

coordinates(:,1) = coordinates(:,1) + correction(1);
coordinates(:,2) = coordinates(:,2) + correction(2);

h1 = figure('color','w','Name',file,'Position',[round((screensize(3)-1600)/2)...
    round((screensize(4)-600)/2) 1600 600],'visible',visfig);
subplot(1,4,1),
hold on
scatter(traj_coordinates(:,1),traj_coordinates(:,2),'bo');
scatter(coordinates(:,1),coordinates(:,2),'*','red');
axis image;
axis off;
title('All coordinates');

% subplot(1,4,2),
% hold on
% scatter(traj_coordinates(:,1),traj_coordinates(:,2),'bo');
% scatter(coordinates(:,1),coordinates(:,2),'*','red');
% axis equal;
% axis off;
% title('All coordinates');
%% match frozen objects with trajectories
tol = 1;
resdist = 0;
iter = 1;

while tol>10^-10 && iter<20

    dx = repmat(coordinates(:,1),1,size(traj_coordinates,1))-repmat(traj_coordinates(:,1),1,size(coordinates,1))';
    dy = repmat(coordinates(:,2),1,size(traj_coordinates,1))-repmat(traj_coordinates(:,2),1,size(coordinates,1))';

    D = dx.^2 + dy.^2;
    [~, idx] = min(D,[],1);
    [~, idy] = min(D,[],2);

    %     ex = find(idx(idy') == 1:length(idy));
    ey = find(idy(idx)' == 1:length(idx));

    %     figure,
    %     hold on
    %     scatter(traj_coordinates(idy(ex),1),traj_coordinates(idy(ex),2));
    %     scatter(coordinates(idx(ey),1),coordinates(idx(ey),2),'*','red');
    %     axis square;
    %     hold off

    ii = sub2ind(size(D),idx(ey)',idy(idx(ey)));
    [f, xi] = ksdensity([dx(ii),dy(ii)]);
    [~, fi] = max(f);

    xcorrection = xi(fi,1);
    ycorrection = xi(fi,2);

    coordinates(:,1) = coordinates(:,1) - xcorrection;
    coordinates(:,2) = coordinates(:,2) - ycorrection;

    tol = (resdist - sqrt(xcorrection^2+ycorrection^2))^2;
    resdist = sqrt(xcorrection^2+ycorrection^2);
    iter = iter+1;

end

%%
selected_matches = idx(ey);
X = [dx(ii), dy(ii)];
opts = statset('MaxIter',1000);
gm = fitgmdist(X,10,'RegularizationValue',10^-5,'Options',opts,'Replicates',5);
[~, gi] = max(gm.ComponentProportion);
X = X - gm.mu(gi,:);
X = (X(:,1)/sqrt(gm.Sigma(1,1,gi))).^2 + (X(:,2)/sqrt(gm.Sigma(2,2,gi))).^2;
thresh = chi2inv(0.99,2);
indgm = X<thresh;
% CDF = mvncdf(X,gm.mu(gi,:),gm.Sigma(:,:,gi));
% Xpost = gm.posterior(X);
% indgm = Xpost(:,gi)>0.8;
% indgm = CDF > 0.025 & CDF < 0.975;
selected_matches = selected_matches(indgm);

% sum(indgm)/numel(indgm)
%%

a1 = idy(idx(ey));
a2 = idx(ey);

coordinates = coordinates + gm.mu(gi,:);

subplot(1,4,2),
hold on
scatter(traj_coordinates(a1(~indgm),1),traj_coordinates(a1(~indgm),2),'bo');
scatter(coordinates(a2(~indgm),1),coordinates(a2(~indgm),2),'*','red');
axis image;
axis off;
title(sprintf('Rejected matches (%0.0f%%)',sum(~indgm)/numel(indgm)*100));

subplot(1,4,3),
hold on
scatter(traj_coordinates(a1(indgm),1),traj_coordinates(a1(indgm),2),'bo');
scatter(coordinates(a2(indgm),1),coordinates(a2(indgm),2),'*','red');
axis image;
axis off;
title(sprintf('Accepted matches (%0.0f%%)',sum(indgm)/numel(indgm)*100));

subplot(1,4,4),
hold on,
for i = a1(indgm)'
    plot(tracks(i).x,tracks(i).y);
end
scatter(coordinates(a2(indgm),1),coordinates(a2(indgm),2),5,'*');
axis image;
axis off;
title(sprintf('Matched trajectories (%d)',sum(indgm)));

% h2 = figure('color','w','Name',file,'Position',[round((screensize(3)-1800)/2)...
%     round((screensize(4)-1200)/2) 1800 1200],'visible',visfig);

%% analyze object fluorescence levels
fi = 0;
images = {};
p = [pictures.Phase];
p = prctile(p(:),[0.1 99.9]);

for i = 1:size(pictures,2)

    if sum(selected_matches==i)>0

        sumYFP = NaN;
        sumCFP = NaN;
        sumRFP = NaN;
        bgYFP = NaN;
        bgCFP = NaN;
        bgRFP = NaN;
        spotYFP = NaN;
        spotCFP = NaN;
        spotRFP = NaN;
        spotintYFP = NaN;
        spotintCFP = NaN;
        spotintRFP = NaN;
%         spotCoordCFP = NaN;
%         spotCoordYFP = NaN;
%         spotCoordRFP = NaN;
        volume = NaN;

        mask = FASTcellmask(pictures(i).Phase);

        if sum(mask(:))>0

            ang = regionprops(mask,'Orientation');
            hcell = imrotate(mask,-ang.Orientation);
            volume = sum((sum(hcell,1)./2.*pixelSizePics).^2.*pi.*pixelSizePics);
            b = bwboundaries(mask);
            boundary = bwperim(mask);
            b = b{1};

            Phase = cat(3, pictures(i).Phase, pictures(i).Phase, pictures(i).Phase);
            Phase = (Phase-p(1))/(p(2)-p(1));
            %             m = cat(3, boundary, boundary, boundary);
            %             Phase(m) = 0;
            m = cat(3, boundary, false & boundary, false & boundary);
            Phase(m) = 1;
            %             imshow(Phase)
            CFP = 0*Phase;
            YFP = 0*Phase;
            RFP = 0*Phase;

            if ~isempty(pictures(i).YFP)
                [sumYFP, bgYFP, bgSTD, spotYFP, spotintYFP, spotCoordYFP] = FASTcellFluorescence(mask, pictures(i).YFP);
                YFP = pictures(i).YFP;
                YFP = (YFP-bgYFP)/(bgSTD*10);
                YFP = cat(3, YFP, YFP, YFP);
                YFP(m) = 0.6;
                if ~isnan(spotRFP) && spotRFP>0
                    overlay = [spotCoordYFP(:,1)+1, spotCoordYFP(:,2)];
                    overlay = cat(1,overlay,[spotCoordYFP(:,1)-1, spotCoordYFP(:,2)]);
                    overlay = cat(1,overlay,[spotCoordYFP(:,1), spotCoordYFP(:,2)+1]);
                    overlay = cat(1,overlay,[spotCoordYFP(:,1), spotCoordYFP(:,2)-1]);
                    spotidx = sub2ind(size(YFP),overlay(:,1),overlay(:,2),2*ones(size(overlay,1),1));
                    YFP(spotidx) = 1;
                    spotidx = sub2ind(size(YFP),overlay(:,1),overlay(:,2),1*ones(size(overlay,1),1));
                    YFP(spotidx) = 0;
                    spotidx = sub2ind(size(YFP),overlay(:,1),overlay(:,2),3*ones(size(overlay,1),1));
                    YFP(spotidx) = 0;
                end
            end

            if ~isempty(pictures(i).CFP)
                [sumCFP, bgCFP, bgSTD, spotCFP, spotintCFP, spotCoordCFP] = FASTcellFluorescence(mask, pictures(i).CFP);
                CFP = pictures(i).CFP;
                CFP = (CFP-bgCFP)/(bgSTD*10);
                CFP = cat(3, CFP, CFP, CFP);
                CFP(m) = 0.6;
                if ~isnan(spotRFP) && spotRFP>0
                    overlay = [spotCoordCFP(:,1)+1, spotCoordCFP(:,2)];
                    overlay = cat(1,overlay,[spotCoordCFP(:,1)-1, spotCoordCFP(:,2)]);
                    overlay = cat(1,overlay,[spotCoordCFP(:,1), spotCoordCFP(:,2)+1]);
                    overlay = cat(1,overlay,[spotCoordCFP(:,1), spotCoordCFP(:,2)-1]);
                    spotidx = sub2ind(size(CFP),overlay(:,1),overlay(:,2),2*ones(size(overlay,1),1));
                    CFP(spotidx) = 1;
                    spotidx = sub2ind(size(CFP),overlay(:,1),overlay(:,2),1*ones(size(overlay,1),1));
                    CFP(spotidx) = 0;
                    spotidx = sub2ind(size(CFP),overlay(:,1),overlay(:,2),3*ones(size(overlay,1),1));
                    CFP(spotidx) = 0;
                end
            end

            if ~isempty(pictures(i).RFP)
                [sumRFP, bgRFP, bgSTD, spotRFP, spotintRFP, spotCoordRFP] = FASTcellFluorescence(mask, pictures(i).RFP);
                RFP = pictures(i).RFP;
                RFP = (RFP-bgRFP)/(bgSTD*10);
                RFP = cat(3, RFP, RFP, RFP);
                RFP(m) = 0.6;
                if ~isnan(spotRFP) && spotRFP>0
                    overlay = [spotCoordRFP(:,1)+1, spotCoordRFP(:,2)];
                    overlay = cat(1,overlay,[spotCoordRFP(:,1)-1, spotCoordRFP(:,2)]);
                    overlay = cat(1,overlay,[spotCoordRFP(:,1), spotCoordRFP(:,2)+1]);
                    overlay = cat(1,overlay,[spotCoordRFP(:,1), spotCoordRFP(:,2)-1]);
                    spotidx = sub2ind(size(RFP),overlay(:,1),overlay(:,2),2*ones(size(overlay,1),1));
                    RFP(spotidx) = 1;
                    spotidx = sub2ind(size(RFP),overlay(:,1),overlay(:,2),1*ones(size(overlay,1),1));
                    RFP(spotidx) = 0;
                    spotidx = sub2ind(size(RFP),overlay(:,1),overlay(:,2),3*ones(size(overlay,1),1));
                    RFP(spotidx) = 0;
                end
            end

            if (~isnan(sumCFP) || ~isnan(sumYFP) || ~isnan(sumRFP)) && fi < 108

                %                 width = 30 + max(max(b(:,1))-min(b(:,1)), max(b(:,2))-min(b(:,2)));
                width = 60;
                xmin = (max(b(:,2))+min(b(:,2))-width)/2;
                ymin = (max(b(:,1))+min(b(:,1))-width)/2;
                rect = [xmin ymin width width];

                images = [images {imcrop(Phase,rect)}];
                images = [images {imcrop(CFP,rect)}];
                images = [images {imcrop(YFP,rect)}];
                images = [images {imcrop(RFP,rect)}];

                %                 subplot(8,12,mod(fi*4,96)+1,'replace'),
                %                 imshow(pictures(i).Phase',[]);
                %                 hold on
                %                 plot(b(:,1),b(:,2),':r','linewidth',1.5)
                % %                 xlabel(volume);
                %                 xlim([xlim1 xlim2]);
                %                 ylim([ylim1 ylim2]);
                % %                 title('Phase');
                %                 subplot(8,12,mod(fi*4,96)+2,'replace'),
                %                 imshow(pictures(i).CFP',[]);
                %                 hold on
                %                 plot(b(:,1),b(:,2),':r','linewidth',1.5)
                %                 xlim([xlim1 xlim2]);
                %                 ylim([ylim1 ylim2]);
                % %                 xlabel(spotCFP);
                % %                 title('CFP');
                %                 subplot(8,12,mod(fi*4,96)+3,'replace'),
                %                 imshow(pictures(i).YFP',[]);
                %                 hold on
                %                 plot(b(:,1),b(:,2),':r','linewidth',1.5)
                %                 xlim([xlim1 xlim2]);
                %                 ylim([ylim1 ylim2]);
                % %                 xlabel(spotYFP);
                % %                 title('YFP');
                %                 subplot(8,12,mod(fi*4,96)+4,'replace'),
                %                 imshow(pictures(i).RFP',[]);
                %                 hold on
                %                 plot(b(:,1),b(:,2),':r','linewidth',1.5)
                %                 xlim([xlim1 xlim2]);
                %                 ylim([ylim1 ylim2]);
                % %                 xlabel(spotRFP);
                % %                 title('RFP');

                fi = fi + 4;
            end

            %     else
            %         figure,
            %         imshow(pictures(i).Phase,[]);
        end

        tracks(idy(i)).sumYFP = sumYFP;
        tracks(idy(i)).sumRFP = sumRFP;
        tracks(idy(i)).sumCFP = sumCFP;

        tracks(idy(i)).spotYFP = spotYFP;
        tracks(idy(i)).spotRFP = spotRFP;
        tracks(idy(i)).spotCFP = spotCFP;

        tracks(idy(i)).spotintYFP = spotintYFP;
        tracks(idy(i)).spotintRFP = spotintRFP;
        tracks(idy(i)).spotintCFP = spotintCFP;

        tracks(idy(i)).bgYFP = bgYFP;
        tracks(idy(i)).bgRFP = bgRFP;
        tracks(idy(i)).bgCFP = bgCFP;

        %         tracks(idy(i)).spotCoordCFP = spotCoordCFP;
        %         tracks(idy(i)).spotCoordYFP = spotCoordYFP;
        %         tracks(idy(i)).spotCoordRFP = spotCoordRFP;

        tracks(idy(i)).MajorAxisLength = pictures(i).MajorAxisLength;
%         tracks(idy(i)).Solidity = pictures(i).Solidity;
%         tracks(idy(i)).NormVar = pictures(i).NormVar;
%         tracks(idy(i)).LogL = pictures(i).LogL;

        tracks(idy(i)).Volume = volume;
        tracks(idy(i)).Pictime = pictures(i).time;
    end
end

h2 = figure('color','w','Name',file,'visible',visfig);
imdisp(images,'Size',[9 12],'Border',[0.01 0.01]);

ind = cellfun(@isempty,{tracks.sumCFP}) & cellfun(@isempty,{tracks.sumYFP}) & cellfun(@isempty,{tracks.sumRFP});
tracksFAST = tracks(~ind);

ind = isnan([tracksFAST.sumCFP]) & isnan([tracksFAST.sumYFP]) & isnan([tracksFAST.sumRFP]);
tracksFAST = tracksFAST(~ind);

metadata0 = getMetadataTracks(tracks0);
metadataFAST = cell(numel(tracksFAST),1);
[metadataFAST{:}] = deal('FAST');

tracks = tracksFAST;

rf = setdiff(fieldnames(tracksFAST),fieldnames(tracks0));
tracksFAST = rmfield(tracksFAST,rf);

if contains(varargin,'hpcc')
    plotTrackStats([tracks0 tracksFAST],unique([metadata0; metadataFAST]),[metadata0; metadataFAST],[],strcat(fileInfo.name_,'.FAST_tracks'),'png','hpcc');
else
    plotTrackStats([tracks0 tracksFAST],unique([metadata0; metadataFAST]),[metadata0; metadataFAST],[],strcat(fileInfo.name_,'.FAST_tracks'),'png');
end

leftover_pics = pictures(setdiff(1:numel(pictures),selected_matches));
matched_pics =  pictures(selected_matches);

set(h1,'PaperPositionMode','auto')
print(h1,strcat(file,'.FAST_matches.png'),'-dpng','-r150','-noui');
% print(h1,strcat(file,'.FAST_matches.svg'),'-dsvg','-r150','-noui');

set(h2,'PaperPositionMode','auto')
print(h2,strcat(file,'.FAST_cells.png'),'-dpng','-r150','-noui');
% print(h2,strcat(file,'.FAST_cells.svg'),'-dsvg','-r150','-noui');

save(fileInfo.FASTfile, 'tracks', 'matched_pics', 'leftover_pics');

crop(strcat(file,'.FAST_matches.png'));
crop(strcat(file,'.FAST_cells.png'));

fprintf('done.\n');

%%

% if do_plots
%     temp = tracks(~cellfun('isempty',{tracks.sumRFP}));
%     tempfluo = temp(~isnan([temp.sumRFP]));
%     tempfluo = tempfluo([tempfluo.sumRFP]./[tempfluo.sumRFP]>0);
%     tempfluo = tempfluo([tempfluo.trajtime]>30);
%     fratio = log2([tempfluo.sumRFP]./[tempfluo.sumRFP]);
%     fratio = (fratio - min(fratio))/max(fratio - min(fratio));
%
%     figure,
%     imshow(zeros(ceil(imageWidth*pixelSizeTracks),ceil(imageHeight*pixelSizeTracks)));
%     hold on;
%     for i = 1:numel(fratio)
%         plot(tempfluo(i).x,tempfluo(i).y,'color',[1-fratio(i) fratio(i) 0],'linewidth',2);
%         scatter(tempfluo(i).x(end),tempfluo(i).y(end),'MarkerFaceColor',[1-fratio(i) fratio(i) 0],'MarkerEdgeColor',[1-fratio(i) fratio(i) 0]);
%         %         scatter(tempfluo(i).x(tempfluo(i).tumble == 1),tempfluo(i).y(tempfluo(i).tumble == 1),'+','MarkerEdgeColor','w');
%     end
%     hold off;
%
%     figure,
%     scatter(log10([tempfluo.fit_runtime]) , fratio(~cellfun('isempty',{tempfluo.fit_runtime})));
%     xlabel('Log10 run length (s)');
%     ylabel('Log2 YFP/RFP');
% end
end
