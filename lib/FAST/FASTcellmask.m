function mask = FASTcellmask(phase,varargin)

bw = edge(phase,'log');
bw = imfill(bw,'holes');
prop = regionprops(bw,'Area','Centroid');
[~, mi] = max([prop.Area]);

mask = bwselect(bw,prop(mi).Centroid(1),prop(mi).Centroid(2),4);

if ismember('plot',varargin)
    bound = bwboundaries(mask);
    imshow(phase,[],'InitialMagnification',500);
    if numel(mask)>0
        hold on,
        plot(bound{1}(:,2),bound{1}(:,1),':g','linewidth',2);
    end
end