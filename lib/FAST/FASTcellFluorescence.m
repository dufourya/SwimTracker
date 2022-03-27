function [cellFluo, bgFluo, bgSTD, spotnb, spotint, spotCoord] = FASTcellFluorescence(mask, img, varargin)

img = img-100;

maskdil = imdilate(mask,strel('sphere',5));

cellFluo = img(maskdil);
bgFluo = img(~maskdil);
a=[];
b=[];

H = swtest(bgFluo(1:ceil(numel(bgFluo)/2000):end));

if H
    cellFluo = NaN;
    bgFluo = NaN;
    spotnb = NaN;
    spotint = NaN;
    bgSTD = NaN;
    spotCoord = NaN;
else
    bgSTD = std(bgFluo);
    bgFluo = mean(bgFluo);
    cellFluo = cellFluo - bgFluo;
    cellFluo = sum(cellFluo);
    
    flr_img_hat = cwtft2(img,'wavelet','mexh','Scales',2);
    
    flr_img_filtered = squeeze(flr_img_hat.cfs(:,:,1,1));
    regional_max = imregionalmin(flr_img_filtered);

    param = evfit(flr_img_filtered(~maskdil & regional_max));

    [spot_cdf, idx] = sort(evcdf(flr_img_filtered(maskdil & regional_max),param(1),param(2)));
    ind = spot_cdf*numel(spot_cdf)./(1:numel(spot_cdf))'<= 0.05;

    [a, b] = find(maskdil & regional_max);

    a = a(idx(ind));
    b = b(idx(ind));
    
    spotnb = numel(a);
    spotCoord = [a, b];
    spotint = -flr_img_filtered(sub2ind(size(flr_img_filtered),a,b));
end

if ismember('plot',varargin)
    boundcell = bwboundaries(mask);
    bound = bwboundaries(maskdil);
    imshow(img,[],'InitialMagnification',500);
    if numel(bound)>0 && ~H
        hold on,
        plot(boundcell{1}(:,2),boundcell{1}(:,1),':g','linewidth',2);
        plot(bound{1}(:,2),bound{1}(:,1),'r','linewidth',2);
        plot(b,a,'r*');
    end
end
end
