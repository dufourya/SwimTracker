function tracksTable = reduceTrackStruct(tracks)

for i = 1:numel(tracks)
    tracks(i).centroid_x = mean(tracks(i).x);
    tracks(i).centroid_y = mean(tracks(i).y);
    tracks(i).mean_xvec = mean(tracks(i).xvec);
    tracks(i).mean_yvec = mean(tracks(i).yvec);
end

fields = fieldnames(tracks);
ind = true(numel(fields),1);
for k = 1:numel(fields)
    if sum(cellfun(@numel,{tracks.(fields{k})}) > 1) == 0
        ind(k) = false;
    end
end

tracks = rmfield(tracks,fields(ind));
metadata = getMetadataTracks(tracks);
tracks = rmfield(tracks,'metadata');
tracks = rmfield(tracks,'trajID');
tracks(1).metadata = metadata{1};

[tracks.metadata] = metadata{:};

tracksTable = struct2table(tracks);