function [tracks, good_pics, leftover_pics] = analyze_freezing_assay(file, do_plots)
%% parameters

cutoff = 1;
% do_plots = 1;

%% load data file

fileInfo = FileInfo(file);

load(strcat(fileInfo.name, RadFileInfo.swimtrackerExt),'tracks');
%load(strcat(fileInfo.name, RadFileInfo.swimtrackerExt),'tracks');
metadataTracks = Metadata(fileInfo.meta_file);
metadataTracks.read();
metadataPics = Metadata(strcat(fileInfo.name, FileInfo.picsMetaExt));
metadataPics.read();
load(strcat(fileInfo.name, FileInfo.picturesExt),'pictures');

pixelSizeTracks = metadataTracks.getPixelSize();
pixelSizePics = metadataPics.getPixelSize();

imageHeight = metadataTracks.getImageHeightPixel();
imageWidth = metadataTracks.getImageWidthPixel();

xcorrection = metadataPics.getXCorrection();
ycorrection = metadataPics.getYCorrection();
%%
%tracks = likelihood_tumblebias(tracks, do_plots);
%% find trajectories that finish at the freezing time
end_times = cellfun(@max,{tracks.time},'UniformOutput', false);
t = [end_times{:}] == max([end_times{:}]);
tracks = tracks(t);
traj_coordinates = zeros(numel(tracks),2);

for i = 1:numel(tracks)
    traj_coordinates(i,:) = [tracks(i).x(end), tracks(i).y(end)];
end
%
coordinates = cat(2,pictures.XYcoord);
coordinates = coordinates';

xcorrection = 38;
ycorrection = -46;
coordinates(:,1) = (imageWidth/2*pixelSizeTracks) + coordinates(:,1) + xcorrection;
coordinates(:,2) = (imageHeight/2*pixelSizeTracks) + coordinates(:,2) + ycorrection;

% figure,
% scatter(coordinates(:,1), coordinates(:,2),'o');
% hold on,
% scatter(traj_coordinates(:,1), traj_coordinates(:,2),'*');
% axis equal
%% align two sets of coordinates
% if min(coordinates(:,1))<1
%     coordinates(:,1) = coordinates(:,1) - min(coordinates(:,1)) + 1;
% end
% if min(coordinates(:,2))<1
%     coordinates(:,2) = coordinates(:,2) - min(coordinates(:,2)) + 1;
% end
%{
mat1 = zeros(ceil(pixelSizeTracks * imageWidth)+100);
mat2 = mat1;
m1 = min([traj_coordinates(:,1); coordinates(:,1)]);
m2 = min([traj_coordinates(:,2); coordinates(:,2)]);
indx = sub2ind(size(mat1),1+round(coordinates(:,1)-m1)',1+round(coordinates(:,2)-m2)');
mat1(indx) = 1;
indx = sub2ind(size(mat2),1+round(traj_coordinates(:,1)-m1)',1+round(traj_coordinates(:,2)-m2)');
mat2(indx) = 1;

se = strel('disk',5,4);
H = fspecial('disk',5);

mat1 = imfilter(imdilate(mat1,se),H,'replicate');
mat2 = imfilter(imdilate(mat2,se),H,'replicate');

figure,
imshowpair(mat1,mat2);

[optimizer, metric] = imregconfig('monomodal');
tform = imregtform(mat1,mat2,'translation',optimizer,metric,'DisplayOptimization',1);

movmat2 = imwarp(mat2,tform,'OutputView',imref2d(size(mat1)));
figure,
imshowpair(mat1,movmat2);
title(file);

return

coordinates(:,1) = coordinates(:,1) + tform.T(3,2);
coordinates(:,2) = coordinates(:,2) + tform.T(3,1);
%}
%% match frozen objects with trajectories
dx = repmat(coordinates(:,1),1,size(traj_coordinates,1))-repmat(traj_coordinates(:,1),1,size(coordinates,1))';
dy = repmat(coordinates(:,2),1,size(traj_coordinates,1))-repmat(traj_coordinates(:,2),1,size(coordinates,1))';

D = dx.^2 + dy.^2;

[minx,  idx] = min(D,[],1);
[miny,  idy] = min(D,[],2);

ex = find(idx(idy') == 1:length(idy));
ey = find(idy(idx)' == 1:length(idx));

%mean(sqrt([miny(ex)' minx(ey)]))
%%
max_dist = mean(sqrt([minx(ey) miny(ex)'])) + cutoff;

ey(sqrt(minx(ey))>max_dist) = [];
ex(sqrt(miny(ex))>max_dist) = [];

if do_plots
    %%
    figure,
    scatter(coordinates(:,1), coordinates(:,2),'o');
    axis equal
    hold on
    scatter(traj_coordinates(:,1), traj_coordinates(:,2),'*');
    hold off
    %%
    figure,
    hist(sqrt([miny(ex)' minx(ey)]),50)
    
    figure,
    hold on
    scatter(traj_coordinates(idy(ex),1),traj_coordinates(idy(ex),2));
    scatter(coordinates(idx(ey),1),coordinates(idx(ey),2),'+','red');
    hold off
    
    figure,
    imshow(zeros(ceil(imageWidth*pixelSizeTracks),ceil(imageHeight*pixelSizeTracks)));
    hold on;
    r = rand(length(ey),3);
    
    for i = 1:length(ey)
        plot(tracks(ey(i)).x,tracks(ey(i)).y,'color',r(i,:));
        scatter(coordinates(idx(ey(i)),1),coordinates(idx(ey(i)),2),'MarkerFaceColor',r(i,:),'MarkerEdgeColor',r(i,:));
    end
    hold off;
    figure, fi = 0;
end


%% analyze object fluorescence levels and sizes from 100x pictures

for i = 1:size(pictures,2)
    
    bw = edge(pictures(i).Phase,'log');
    bw = imfill(bw,'holes');
    ang = regionprops(bw, 'Orientation','Eccentricity','MajorAxisLength','Area','Centroid');
    [~, mi] = max([ang.Area]);
    bw = bwselect(bw,ang(mi).Centroid(1),ang(mi).Centroid(2));
    r = imrotate(bw,-ang(mi).Orientation);
    volume = sum((sum(r,1)./2.*pixelSizePics).^2.*pi.*pixelSizePics);
    bw1 = imdilate(bw,strel('sphere',10));
    axislength = ang(mi).MajorAxisLength * pixelSizePics;
    
    YFP = NaN;
    CFP = NaN;
    RFP = NaN;
    bgYFP = NaN;
    bgCFP = NaN;
    bgRFP = NaN;
    
    if ~isempty(pictures(i).YFP)
        bgYFP = mean(pictures(i).YFP(~bw1));
        YFP = (sum(pictures(i).YFP(bw1)) - (sum(sum(bw1)) * bgYFP));
    end
    
    if ~isempty(pictures(i).CFP)
        bgCFP = mean(pictures(i).CFP(~bw1));
        CFP = (sum(pictures(i).CFP(bw1)) - (sum(sum(bw1)) * bgCFP));
    end
    
    if ~isempty(pictures(i).RFP)
        bgRFP = mean(pictures(i).RFP(~bw1));
        RFP = (sum(pictures(i).RFP(bw1)) - (sum(sum(bw1)) * bgRFP));
    end
    
    if do_plots
        boundary1 = bwboundaries(bw1);
        boundary = bwboundaries(bw);
        if numel(boundary) > 0
            subplot(6,8,mod(fi,48)+1),
            fi = fi +1;
            imshow(pictures(i).RFP',[]);
            b = boundary{1};
            b1 = boundary1{1};
            hold on
            plot(b(:,1),b(:,2),'color','green')
            plot(b1(:,1),b1(:,2),'color','red')
            hold off
        end
        title(volume);
    end
    
    
    if sum(idx(ey)==i)>0
        tracks(idy(i)).sumYFP = YFP;
        tracks(idy(i)).sumRFP = RFP;
        tracks(idy(i)).sumCFP = CFP;
        tracks(idy(i)).bgYFP = bgYFP;
        tracks(idy(i)).bgRFP = bgRFP;
        tracks(idy(i)).bgCFP = bgCFP;
        tracks(idy(i)).axislength = axislength;
        tracks(idy(i)).Phase = pictures(i).Phase;
        tracks(idy(i)).volume = volume;
        tracks(idy(i)).pictime = pictures(i).time;
        tracks(idy(i)).match_dist = sqrt(miny(i));
    end
    
end
%[tracks(:).exp_name] = deal(fileInfo.name);


leftover_pics = pictures(setdiff(1:numel(pictures),idx(ey)));
good_pics =  pictures(idx(ey));
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