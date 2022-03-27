function updateMovieTable(varargin)

files = subdir('*.table.txt');
ind = ~contains({files.name}','/.');
files = files(ind);
fprintf('%d movies in total\n',numel(files));
%%
fprintf('Reading data...');
data = cell(1,numel(files));
metadata = cell(1,numel(files));
for k = 1:numel(data)
    data{k} = readtable(files(k).name,'ReadVariableNames',true,'Delimiter',',');
    data{k}.Properties.Description = strrep(files(k).name,'.table.txt','');
    metafile = strrep(files(k).name,'.table.txt','.meta');
    if exist(metafile,'file')
        metadata{k} = Metadata(metafile);
        metadata{k}.read();
    end
end
clear k;
%%
file_parts = cellfun(@(x) strsplit(x,'/'),{files.name},'UniformOutput',false);
file_students = cellfun(@(x) x{find(strcmp(x,'Microscope'))-1},file_parts,'UniformOutput',false);
file_names = strrep(cellfun(@(x) x{end},file_parts,'UniformOutput',false),'.table.txt','');
file_dates = cellfun(@(x) x, regexp(file_names,'20\d{6}','match','once'),'UniformOutput',false);
file_dates_folder = regexp({files.name},'20\d{2}-\d{2}-\d{2}','match','once');
ind = cellfun(@isempty,file_dates_folder);
if sum(ind)>0
    [file_dates_folder{ind}] = deal('NaT');
end
file_dates(cellfun(@isempty,file_dates)) = strrep(file_dates_folder(cellfun(@isempty,file_dates)),'-','');
%%
OD = strings(numel(file_names),1);
temperature = strings(numel(file_names),1);
catalog = strings(numel(file_names),1);

%%
notes = strings(numel(file_names),1);
buffer = strings(numel(file_names),1);
carbon = strings(numel(file_names),1);
growth = strings(numel(file_names),1);
genotype = strings(numel(file_names),1);
plasmid = strings(numel(file_names),1);
project = strings(numel(file_names),1);
folder = strings(numel(file_names),1);
inducer = strings(numel(file_names),1);
particle = strings(numel(file_names),1);
material = strings(numel(file_names),1);
%%
cumultrajtime = nan(numel(file_names),1);
medianlogdiff = nan(numel(file_names),1);
tenthlogdiff  = nan(numel(file_names),1);
ninetothlogdiff  = nan(numel(file_names),1);

exposure = nan(numel(file_names),1);
fps = nan(numel(file_names),1);
numImages = nan(numel(file_names),1);
pixel = nan(numel(file_names),1);
pixelType = strings(numel(file_names),1);
shutter = strings(numel(file_names),1);
objective = strings(numel(file_names),1);
multiplier = false(numel(file_names),1);
filterBlock = strings(numel(file_names),1);
binning = strings(numel(file_names),1);
FAST = false(numel(file_names),1);
%%
for k = 1:numel(data)
    trajtime = data{k}.trajtime;
    diffcoeff = data{k}.diffcoeff_cve_mean;
    [B,I] = sort(diffcoeff);
    cumultrajtime(k) = sum(trajtime(~isnan(B)));
    cumul = cumsum(trajtime(I))/cumultrajtime(k);
    medianlogdiff(k) = log10(min(B(cumul>=0.5)));
    tenthlogdiff(k) = log10(min(B(cumul>=0.1)));
    ninetothlogdiff(k) = log10(min(B(cumul>=0.9)));
    if ~isempty(metadata{k})
        exposure(k) = metadata{k}.getExposure;
        fps(k) = metadata{k}.getImageIntervalMs;
        pixel(k) = metadata{k}.getPixelSize;
        pixelType(k) = metadata{k}.getPixelType;
        shutter(k) = metadata{k}.getVal('Arduino-Switch','Label');
        objective(k) = metadata{k}.getVal('TINosePiece','Label');
        multiplier(k) = str2double(metadata{k}.getVal('Arduino-Input','AnalogInput0'))>900;
        filterBlock(k) = metadata{k}.getVal('TIFilterBlock1','Label');
        binning(k) = metadata{k}.getVal('Andor Zyla 4.2','Binning');
        file_dates{k} = yyyymmdd(datetime(metadata{k}.getVal('Experiment', 'Date')));
        numImages(k) = metadata{k}.getNumberOfImages;
        FAST(k) = metadata{k}.getNImagesBeforeFreeze<numImages(k);
    end
end
fprintf('done\n');
%%
varNames = {'Researcher','Folder','Date','FileName',...
    'Objective','Multiplier','FilterBlock','Shutter','Exposure','Interval',...
    'NumberImages','FASTassay','PixelSize','PixelType','Binning','CumulativeTime','MedianDiffusion',...
    'Prct10Diffusion','Prct90Diffusion',...
    'StrainNumber','Particle','Plasmid','GenotypeNumber','GrowthMedium','CarbonSource','Temperature',...
    'InducerConcentration','OpticalDensity','MotilityBuffer','Material','ProjectName','ProjectFolder','Notes'};

newT = table(string(file_students)',string(file_dates_folder)',string(file_dates)',...
    string(file_names)',objective,multiplier,filterBlock,shutter,exposure,fps,numImages,...
    FAST,pixel,pixelType,binning,cumultrajtime,...
    medianlogdiff,tenthlogdiff,ninetothlogdiff,...
    catalog,particle,plasmid,genotype,growth,carbon,temperature,inducer,OD,buffer,material,project,folder,notes,'VariableNames',varNames);

newT.Properties.RowNames = cellstr(strcat(newT.Researcher,newT.Folder,newT.FileName));

%%
if exist('Movie_table.xlsx','file')
    T = readtable('Movie_table.xlsx');
    T.Properties.RowNames = strcat(T.Researcher,cellstr(T.Folder),T.FileName);
    rowNames = intersect(T.Properties.RowNames,newT.Properties.RowNames);
    lefoverRows = setdiff(T.Properties.RowNames,newT.Properties.RowNames);
    varNames = {'StrainNumber','Particle','Plasmid','GenotypeNumber','GrowthMedium',...
        'CarbonSource','Temperature','InducerConcentration',...
        'OpticalDensity','MotilityBuffer','Material','ProjectName','ProjectFolder','Notes'};
    varNames = intersect(varNames,T.Properties.VariableNames);
    if numel(rowNames) >0
        newT(rowNames,varNames) = T(rowNames,varNames);
    end
    if ~isempty(lefoverRows)
        T = T(lefoverRows,:);
        writetable(T,'Leftover_table.xlsx');
    end
    movefile 'Movie_table.xlsx' 'Movie_table.old'
end

%%
writetable(newT,'Movie_table.xlsx');

end
