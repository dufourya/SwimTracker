function [metadataAll, uniqueMeta] = getMetadataTracks(tracksAll)

if istable(tracksAll)
    metadataAll = tracksAll.metadata;
else
    metadataAll = cell(numel(tracksAll),1);
    for i = 1:numel(tracksAll)
        metadataAll{i} = tracksAll(i).metadata.name_;
    end
end
uniqueMeta = unique(metadataAll);
%%