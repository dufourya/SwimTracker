function mmc = Initialize(scope)
    
    swimTrackerFolder = 'Microscope';
    code_path = fullfile(regexprep(userpath,';',''),swimTrackerFolder,'config');
    
    import mmcorej.*;
    mmc = CMMCore;

    if strcmp(scope, 'zenscope')
        config = fullfile(code_path, 'ZENSCOPE_MMConfig_TI.cfg');
    elseif strcmp(scope, 'emoscope')
        config = fullfile(code_path, 'EMOSCOPE_MMConfig_TI.cfg');        
    else
        error('wrong scope name!');
    end
    
    try
        mmc.loadSystemConfiguration(config);
        mmc.waitForSystem();
    catch e
        fprintf('\n\nCaught error while loading sys config:\n %s\n', getReport(e, 'extended'));
        try
            fprintf('Attempting to reset system...\n');
            mmc.reset();
            mmc.waitForSystem();
            error('System reset\n');
        catch e2
            error('\n\nCaught error while resetting:\n%s\n', getReport(e, 'extended'));
        end
    end
    fprintf('Loading sys config successful!\n');
    mmc.setConfig('System','Startup');
	mmc.waitForSystem();
    checkTemperature(mmc);
    
end