function [c, p, t] = startMatlabPool()
    % starts the parallel pool that works on both the cluster and 
    % local computers.
    
    if (strcmp(computer(),'GLNXA64'))
        % the os is linux-like, so we might be on the cluster.
        [~, hostname] = system('echo $HOSTNAME');
        if isempty(strfind(hostname, '.local'))
            % We are on a local computer. Let matlab decide
            % whether we need to start a new pool.
            gcp();
        else
            % We are on the cluster. Help matlab manage its workers
            % correctly.
            c=parcluster();
            t=tempname();
            mkdir(t);
            c.JobStorageLocation=t;
            delete(gcp('nocreate'));
            parpool(c);
        end
    end
end
