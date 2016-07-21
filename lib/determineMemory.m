function [freemem, memsize] = determineMemory()
    f_name = getFunctionName();
    os = computer();
    if (strcmp(os,'GLNXA64'))
        %http://stackoverflow.com/questions/12350598/how-to-access-memory-information-in-matlab-on-unix-equivalent-of-user-view-max
        [r,w] = unix('free -m | grep Mem');
        stats = str2double(regexp(w, '[0-9]*', 'match'));
        memsize = stats(1);
        freemem = (stats(3)+stats(end));
    elseif (strcmp(os,'PCWIN64') || strcmp(os, 'PCWIN'))
        [usr, sys] = memory();
        freemem = sys.PhysicalMemory.Available / 1024^2;
        memsize = sys.PhysicalMemory.Total / 1024^2;
    elseif (strcmp(os, 'MACI64'))
        % modified from http://apple.stackexchange.com/a/48195
        page_size_Mb = 4096/1024^2;

        [~,free]        = unix('vm_stat | grep ''free'' | awk ''{print $3}'' | sed ''s/\.//''');
        [~,inactive]    = unix('vm_stat | grep ''inactive'' | awk ''{print $3}'' | sed ''s/\.//''');
        [~,speculative] = unix('vm_stat | grep ''speculative'' | awk ''{print $3}'' | sed ''s/\.//''');

        free        = str2double(free);
        inactive    = str2double(inactive);
        speculative = str2double(speculative);

        freemem  = (free + speculative) * page_size_Mb;
        memsize  = (free + inactive) * page_size_Mb;
    else
        error(['[' f_name ']' 'don''t recognize operating system ''' os '''']);
    end
end