function [objs] = fo4_rp_ysd(img, imgun, mask, thresh, ste, removeHalo)

diagnostic_plot = 0;  % for debugging -- plot things.
regmax = imregionalmax(img);
regmax(~mask) = 0;
thresh = quantile(img(regmax),1-thresh);
% img = abs(img);
isbright = (img >= thresh);

% BW = imregionalmax(img) & isbright;

dimg = imdilate(img, ste);

BW = (img == dimg) & isbright & regmax & mask;

% imshow(BW,[]);

% Bandpass filter (if nsize > 0) -- use Grier et al bpass.m
% if bpfiltsize>0
%     noisesize = 1;  % size of noise, px
%     filtimg = bpass(abs(img),noisesize,bpfiltsize);
% else
% filtimg = bpass(abs(img),1,ceil(nsize)*2+1);
% filtimg = imgaussfilt(abs(img));
% imshow(img,[]);
% pause();
% imshow(filtimg,[]);
% pause;
% end

% Three options for thresholding -- see above
% For the chosen option, determine the points that "pass" the threshold.
% Move to a separate function, so it can be called by the GUI

% dimg = imdilate(filtimg, ste);
% imshow(dimg,[]);
% pause;
% Local maxima in the filtered image
% BW = (filtimg == dimg);
% imshow(BW,[]);
% pause;
% find local maxima that are also above-threshold
% isbright = (dimg > abs(thresh));

% look at each local maximum
[y, x] = find(BW);
% if numel(x)>1000
%     [~,o] = sort(img(BW),'descend');
%     x=x(o(1:1000));
%     y=y(o(1:1000));
% end

% Get rid of maxima too close to the edge
[lenx, leny] = size(ste.Neighborhood);  % 'floor' isn't really necessary, but this
% lenx = ceil(2*scalewavelet);
% leny = lenx;
% is the size of "nhood = getnhood(ste);" for a disk structuring
% element of that size
% lenx = ceil(1.5*lenx);
% leny = ceil(1.5*leny);

edgeind = x < lenx | x > (size(img,2) - lenx) | y < leny | y > (size(img,1) - leny);
x(edgeind) = [];
y(edgeind) = [];

if diagnostic_plot
    %         sc = get( groot, 'Screensize' );
    %         h = figure('color','w','Name','Diagnostics','Position',...
    %             [100 100 sc(3)-200 sc(4)-200]);
    subplot(1,2,1),
    img_comp = cat(3,img.*double(mask),img.*double(mask),img);
    img_comp = imadjust(img_comp,stretchlim(img_comp,0));
    imshow(img_comp);
    hold on,
    plot(x,y,'g.');
    hold off;
    xlabel('Frame');
    %     subplot(1,3,2),
    %     imshow(dimg,[]);
    %     hold on,
    %     plot(x,y,'go');
    %     hold off;
    %     xlabel('Dilated image');
    subplot(1,2,2),
    s = std(imgun(:));
    m = mean(imgun(:));
    imshow(imgun,[m-3*s m+3*s]);
    xlabel('Non filtered');
    hold on,
    plot(x,y,'g.');
    hold off;
    drawnow;
    %     pause;
end

% Compute "masses"
savemass = zeros(1, length(x));
rect = zeros(length(x), 4);
savecorr = ones(length(x),2);

% removeHalo = 1;

% Compute the first neighborhood of a maximum, to know the image size
% Can skip if there are no objects to find
if ~isempty(x)
    cropimg = zeros(lenx,leny,length(x));
    for k = 1:length(x)
        rect(k,:) = [round(x(k) - lenx/2) round(y(k) - leny/2) (lenx-1) (leny-1)];
        tmp = imcrop(imgun, rect(k,:));
        cropimg(:,:,k) = tmp;
        savemass(k) = max(reshape(imcrop(img, rect(k,:)),1,[]));
        if removeHalo
        rtmp = rot90(tmp);
        [savecorr(k,1), savecorr(k,2)] = corr(tmp(:),rtmp(:));
        end
        %         savecorr(k) = savecorr(k)/4;
        %                 figure;
        %                 imshow(tmp,[],'InitialMagnification',500);
        %                 disp(savecorr(k))
        %                 drawnow;
        %                 pause;
        
    end
else
    objs = [];
    return;
end


if removeHalo
    FDR = 0.05;
    p = sort(savecorr(:,2));
    p = p(find((p' .* (1:numel(p)))./(1:numel(p))<FDR,1,'last'));
    indsym = savecorr(:,2)<=p & savecorr(:,1)>0;
    cropimg = cropimg(:,:,indsym);
    savemass = savemass(indsym);
    rect = rect(indsym,:);
end

% Do refinement (find center) around each local max
% If there are many points, use radialcenter_stk.m rather than radialcenter.m
% for radial symmetry method -- even faster!
xcent = zeros(1,size(cropimg,3));
ycent = xcent;
sigma = xcent;
lsumx = 1:size(cropimg,2);
lsumy = 1:size(cropimg,1);
Lx = lsumx(end);
Ly = lsumy(end);
% for j = 1:length(x)

% Radial-symmetry based fit -- fast, accurate
% If <10 points, use radialcenter_stk.m for extra speed (avoid
% redundant grid calculations); else radialcenter.m
% ind = savemass<0;
% cropimg(:,:,ind) = -cropimg(:,:,ind);

if numel(savemass) < 10
    for j = 1:numel(savemass)
        [xcent(j), ycent(j), sigma(j)] = radialcenter(cropimg(:,:,j));
    end
else
    [xcent, ycent, sigma] = radialcenter_stk(cropimg);
end

%
% if 0
%     figure,
%     for i = 1:numel(xcent)
%     imshow(cropimg(:,:,i),[]);
%     hold on, plot(xcent(i), ycent(i),'o'), hold off;
%     pause;
%     end
% end
% Is the center within reasonable bounds?
% If not, replace with centroid
% frequency of bad cases ~ 1/100,000 !  (Extremely rare)
% See notes Oct. 26, 2011
% This conditional statement can slow things; Delete?
badcase = find(abs(xcent - Lx/2)>1.5*Lx | abs(ycent - Ly/2)>1.5*Ly);
for j = badcase
    ci = cropimg(:,:,j);
    xcent(j) = sum(sum(ci) .* lsumx) / sum(sum(ci));
    ycent(j) = sum(sum(ci,2) .* lsumy') / sum(sum(ci));
end

%% calculate radial symmetry slow
% savecorr = zeros(1, length(xcent));
% for k = 1:numel(xcent)
%     [cartx, carty] = meshgrid((1:leny)-xcent(k),(1:lenx)-ycent(k));
%     [~, rho] = cart2pol(cartx, carty);
%     tmp = cropimg(:,:,k);
%     gprMdl = fitrgp(rho(:),tmp(:));
% %     ypred = resubPredict(gprMdl);
% %     figure(1);
% %     subplot(1,2,1),
% %     imshow(tmp,[]); hold on;
% %     scatter(xcent(k),ycent(k)); hold off;
% %     subplot(1,2,2),
% %     plot(rho(:),tmp(:),'.'); hold on;
% %     plot(rho(:),ypred(:),'-'); hold off;
%     savecorr(k) = log(gprMdl.KernelInformation.KernelParameters(2)/gprMdl.KernelInformation.KernelParameters(1));
% %     pause;
% end
%%

%end
% center position relative to image boundary
xn = xcent + rect(:,1)' - 1; % -1 is to correct for matlab indexing
yn = ycent + rect(:,2)' - 1;

% indsym = savecorr>2;
% % cropimg = cropimg(:,:,indsym);
% savemass = savemass(indsym);
% % rect = rect(indsym,:);
% savecorr = savecorr(indsym);
% sigma = sigma(indsym);
%
% xn = xn(indsym);
% yn = yn(indsym);

if 0
    savecorr = savecorr(indsym,:);
    figure(1),
    imshow(log(double(imgun-min(imgun(:))+1)),[],'InitialMagnification',100);
    hold on;
    scatter(xn, yn, abs(savecorr(:,1))*50);
    hold off;
    drawnow;
end

objs = zeros(6, length(xn));
objs(1,:) = xn;
objs(2,:) = yn;
objs(3,:) = savemass;
objs(4,:) = 1:numel(xn);
objs(7,:) = sigma;
