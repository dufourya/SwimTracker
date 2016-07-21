function convertExtensions(dry_run, old)
    % convert extensions used by code in Microscope repository.
    % INPUTS:
    %   dry_run  : if '1', prints the changes that would occur, but does not
    %              change any files. If '0', actually rewrites the file
    %              names.
    %   old      : a structure with field names corresponding to the output
    %              of properties('FileInfo'), and values corresponding to
    %              the old extension names.
    
    
    if nargin<1 || isempty(dry_run)
        dry_run = 1;
    end
    
    if nargin<2 || isempty(old)
        % extensions from v1.5 and earlier.
        old.picsMetaExt = '.pictures.meta.txt';
        old.picturesExt = '.pictures.mat';
        old.gradBeforeExt = '.before.tiff';
        old.gradAfterExt = '.after.tiff';
        old.metaExt = '.meta.txt';
        old.backgroundExt = '_background.tiff';
        old.gradBeforeMetaExt = '_before.meta.txt';
        old.gradAfterMetaExt = '_after.meta.txt';
        old.sourceExt = '_source.tiff';
        old.sinkExt = '_sink.tiff';
        old.tracksDat = 'tracks';
        old.tracksFinal = 'tracksFinal';
        old.aviExt = '.avi';
        old.binExt = '.bin';
    end
    
    fields = fieldnames(old);
    files = rdir(strcat('**', filesep(), '*'));

    for i=1:numel(files)
        file = files(i).name;
        %fprintf('%s\n', file);
        for j=1:numel(fields)
            field = fields{j};
            old_ext = old.(field);
            new_ext = FileInfo.(field);

            regex = strcat('(', strrep(old_ext, '.', '\.'), ')');
            %fprintf('regex: %s\n', regex);
            [tokens, matches] = regexp(file, regex, 'tokens', 'match');
            if ~isempty(matches)
                new_file = strrep(file, old_ext, new_ext);
                if ~strcmp(file, new_file)
                    fprintf('Found old extension: %s\n', file);
                    if dry_run
                        fprintf('     Would change to: %s\n', new_file);
                    else
                        fprintf('     Changing to: %s\n', new_file);
                        movefile(file, new_file);
                    end
                    break;
                end
            end
        end
    end
end