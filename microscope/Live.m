function h = Live(mmc,config)
    
    global Live
    display_size_width = 1024;
    display_size_height = 1024;
    
    if exist('Live','var') && isfield(Live,'Figure') && ishandle(Live.Figure)
        figure(Live.Figure);
    else
        Live.Figure = figure('pos',[100 10 ...
            display_size_width ...
            display_size_height+25], ...
            'Toolbar','none', ...
            'Menubar','none','Name','Live','NumberTitle', ...
            'off','IntegerHandle','off');
    end
    
    h = Live.Figure;
    
    mmc.setConfig('System','Startup');
    mmc.setConfig('Channel',config);
    mmc.waitForSystem();
    
    mmc.setAutoShutter(0);
    
    camera = char(mmc.getCameraDevice());
    camera_properties = cell(mmc.getDevicePropertyNames(camera).toArray());
    
    exposure = mmc.getExposure();
    intensity = 0;
    maxintensity = 0;
    
    if sum(strcmp('Gain',camera_properties))>0
        gain = str2double(mmc.getProperty(camera,'Gain'));
    else
        gain = -1;
    end
    
    % devices = cell(mmc.getLoadedDevices().toArray());
    
    arduino = 1;
    mmc.setShutterDevice('Arduino-Shutter');
    
    %if sum(strcmp('Arduino-Shutter',devices))>0
    %    arduino = 1;
    %    mmc.setShutterDevice('Arduino-Shutter');
    %else
    %    arduino = 0;
    %    mmc.setShutterDevice('TIDiaShutter');
    %end
    
    
    imgtmp = zeros(display_size_width,display_size_height);
    width = mmc.getImageWidth();
    height = mmc.getImageHeight();
    exitlive = false;
    zoom = false;
    normalize = false;
    crossmark = false;
    circle = false;
    remove_background = false;
    num_bg = 0;
    background_image = zeros(width,height);
    [xc,yc] = cylinder(150,500);
    xc = display_size_width/2 + round(xc(1,:));
    yc = display_size_height/2 + round(yc(1,:));
    
    epishutter = 0;
    diashutter = 0;
    shot = 1;
    
    prev_stage_pos.x = mmc.getXPosition('XYStage');
    prev_stage_pos.y = mmc.getYPosition('XYStage');
    prev_stage_pos.z = mmc.getPosition('ZStage');
    stage_tic = tic;
    
    blocks = ...
        cell(mmc.getAllowedPropertyValues('TIFilterBlock1','Label').toArray());
    currentblock = ...
        find(strcmp(char(mmc.getProperty('TIFilterBlock1','Label')),blocks));
    strblocks = blocks{1};
    for j = 2:numel(blocks)
        strblocks = strcat(strblocks,'|',blocks{j});
    end
    
    %mmc.setProperty('Arduino-Switch','State', ...
    %                num2str(diashutter + epishutter*2^5));
    
    set(Live.Figure,'CloseRequestFcn',@closeGUI);
    
    fpos = get(Live.Figure,'pos');
    lppos = [0 fpos(4)-25 fpos(3) 25];
    axpos = [0 0 fpos(3) fpos(4)-25];
    
    Live.Panel = uipanel('units','pixels','pos',lppos,'BorderType','none');
    Live.Axes = axes('Parent',Live.Figure,'units','pixels','pos',axpos, ...
        'Color',[1 1 1],'TickDir','out');
    
    Live.StartBtn = uicontrol(Live.Panel,'units','pixels','Position', ...
        [5 5 40 20],'String','Start','callback', ...
        @StartFn);
    Live.StopBtn = uicontrol(Live.Panel,'units','pixels','Position', ...
        [50 5 40 20],'String','Stop','callback',@StopFn);
    
    %Shutter
    uicontrol(Live.Panel,'units','pixels','Position',[90 2 50 20], ...
        'Style','text','String','Shutter');
    Live.shutterDia = uicontrol(Live.Panel,'units','pixels', ...
        'Position', [140 5 30 20], ...
        'Style','togglebutton', ...
        'String','Dia','enable','on', ...
        'callback',@DiaFn);
    Live.shutterEpi = uicontrol(Live.Panel,'units','pixels', ...
        'Position', [170 5 30 20], ...
        'Style','togglebutton', ...
        'String','Epi','enable','on', ...
        'callback',@EpiFn);
    
    %Zoom
    Live.norm = uicontrol(Live.Panel,'units','pixels','Position', ...
        [205 5 35 20],'Style','togglebutton','String', ...
        'Norm','enable','on','Value',false, ...
        'callback',@NormFn);
    Live.zoom = uicontrol(Live.Panel,'units','pixels','Position', ...
        [245 5 35 20],'Style','togglebutton','String', ...
        'Zoom','enable','on','callback',@ZoomFn);
    
    %Gain
    uicontrol(Live.Panel,'units','pixels','Position',[285 2 30 20], ...
        'Style','text','String','Gain','enable','on');
    Live.GainMinus = uicontrol(Live.Panel,'units','pixels','Position', ...
        [315 5 20 20],'String','-', ...
        'callback',@GainFn,'enable','off');
    Live.GainPlus = uicontrol(Live.Panel,'units','pixels','Position', ...
        [340 5 20 20],'String','+', ...
        'callback',@GainFn,'enable','off');
    Live.gaintext = uicontrol(Live.Panel,'units','pixels','Style','edit', ...
        'Position',[365 5 35 20],'String', ...
        num2str(gain),'callback',@GainFn,'enable','off');
    
    if gain >= 0
        set(Live.GainMinus, 'enable', 'on');
        set(Live.GainPlus, 'enable', 'on');
        set(Live.gaintext, 'enable', 'on');
    end
    
    %Exposure
    uicontrol(Live.Panel,'units','pixels','Position',[405 2 25 20], ...
        'Style','text','String','Expo','enable','on');
    Live.ExpoMinus = uicontrol(Live.Panel,'units','pixels','Position', ...
        [430 5 20 20],'String','-', ...
        'callback',@ExpoFn,'enable','on');
    Live.ExpoPlus = uicontrol(Live.Panel,'units','pixels','Position', ...
        [455 5 20 20],'String','+', ...
        'callback',@ExpoFn,'enable','on');
    Live.expotext = uicontrol(Live.Panel,'units','pixels','Style','edit', ...
        'Position',[480 5 40 20], ...
        'String',num2str(round(exposure)), ...
        'callback',@ExpoFn,'enable','on');
    
    %Block
    uicontrol(Live.Panel,'units','pixels','Position',[525 2 30 20], ...
        'Style','text','String','Block');
    Live.block = uicontrol(Live.Panel,'units','pixels', ...
        'Position',[555 5 70 20],'Style','popup', ...
        'String',strblocks,'Callback',@changeBlock, ...
        'enable','on','Value',currentblock);
    
    % Light
    uicontrol(Live.Panel,'units','pixels','Position',[630 2 30 20],'Style', ...
        'text','String','Light','enable','off');
    Live.intensitytext = uicontrol(Live.Panel,'units','pixels','Style', ...
        'edit','Position',[660 5 30 20], ...
        'String',num2str(intensity), ...
        'callback',@LightFn,'enable','off');
    Live.maxintensitytext = uicontrol(Live.Panel,'units','pixels','Style', ...
        'edit','Position',[695 5 30 20], ...
        'String',num2str(maxintensity), ...
        'callback',@LightFn,'enable','off');
    
    % Background subtractino
    Live.backgroung = uicontrol(Live.Panel,'units','pixels','Position', ...
        [735 5 35 20],'Style','togglebutton','String', ...
        'Bkgd','enable','on','Value',false, ...
        'callback',@backgroundSub);
    
    %Snapshot
    Live.snap = uicontrol(Live.Panel,'units','pixels', ...
        'Position', [780 5 60 20],'String','Snapshot', ...
        'Callback',@snapShot,'enable','on');
    
    %cross mark
    Live.crossmark = uicontrol(Live.Panel,'units','pixels', ...
        'Position',[840 5 40 20],'String','Cross', ...
        'Callback',@drawCross,'enable','on');
    
    %circle mark
    Live.circle = uicontrol(Live.Panel,'units','pixels', ...
        'Position',[880 5 40 20],'String','Circle', ...
        'Callback',@drawCircle,'enable','on');
    
    %focus
    Live.focus = uicontrol(Live.Panel,'units','pixels', ...
        'Position',[930 5 40 20],'String','Focus', ...
        'Callback',@fullFocus,'enable','on');
    Live.focustext = uicontrol(Live.Panel,'units','pixels','Style','edit', ...
        'Position',[970 5 40 20], ...
        'String',num2str(0),'callback',@fullFocus, ...
        'enable','off');
    
    
    function backgroundSub(hObject, eventdata);
        remove_background = ~remove_background;
        if ~remove_background
            background_image = 0 * background_image;
            num_bg = 0;
        end
    end
    
    function fullFocus(hObject, eventdata);
        mmc.setAutoShutter(0);
        mmc.setShutterOpen(1);
        mmc.fullFocus();
        set(Live.focustext, ...
            'String',num2str(round(mmc.getCurrentFocusScore())));
        mmc.setShutterOpen(1);
        mmc.waitForImageSynchro();
        mmc.snapImage();
        mmc.setShutterOpen(0);
        imgtmp = mmc.getImage();
        if mmc.getBytesPerPixel == 2
            pixelType = 'uint16';
        else
            pixelType = 'uint8';
        end
        imgtmp = typecast(imgtmp, pixelType);
        img = rotateFrame(imgtmp,mmc);
        
        if width > display_size_width
            if zoom,
                img = imcrop(img, [floor(width/2)-display_size_width/2 ...
                    floor(height/2)-display_size_height/2 ...
                    display_size_width display_size_height]);
            else
                img = imresize(img, display_size_width/width, 'bicubic');
            end
        else
            if zoom,
                img = imresize(img, 2*display_size_width/width, 'bicubic');
                img = imcrop(img, [display_size_width/2 display_size_height ...
                    display_size_width display_size_height]);
            end
        end
        img1 = double(img);
        intensity = round(std(img1(:)));
        maxintensity = sum(img1(:)== max(img1(:)));
        set(Live.intensitytext, 'String',num2str(intensity));
        set(Live.maxintensitytext, 'String',num2str(maxintensity));
        
        if normalize,
            img1 = img1 - min(img1(:));
            img1 = img1 ./ max(img1(:));
        else
            img1 = img1 ./ (2^16-1);
        end
        img1(isnan(img1)) = 0;
        
        Live.Image = imshow(img1,'parent',Live.Axes);
    end
    
    function drawCross(hObject, eventdata);
        crossmark = ~crossmark;
    end
    
    function drawCircle(hObject, eventdata);
        circle = ~circle;
    end
    
    
    function DiaFn(hObject, eventdata)
        diashutter = get(hObject,'Value');
        mmc.setShutterOpen(diashutter);
        background_image = 0 * background_image;
        num_bg = 0;
    end
    
    function EpiFn(hObject, eventdata)
        epishutter =  get(hObject,'Value');
        %         if arduino
        mmc.setProperty('Spectra','White_Enable',1);
        mmc.setProperty('Spectra','White_Level',100);
        mmc.setProperty('Spectra','State',epishutter);
        background_image = 0 * background_image;
        num_bg = 0;
        %         else
        %             mmc.setProperty('TIEpiShutter','State',epishutter);
        %         end
    end
    
    function snapShot(hObject, eventdata)
        % imwrite(uint16(imadjust(img,[])), ...
        %         strcat('snapshot',num2str(shot),'.tiff'),'tiff');
        while exist(strcat('snapshot',num2str(shot,'%0.5d'),'.tiff'),'file')
            shot = shot+1;
        end
        imwrite(rotateFrame(imgtmp,mmc),strcat('snapshot', ...
            num2str(shot,'%0.5d'),'.tiff'),'tiff');
    end
    
    function changeBlock(hObject, eventdata)
        mmc.setProperty('TIFilterBlock1', 'Label', ...
            blocks{get(hObject,'Value')})
        background_image = 0 * background_image;
        num_bg = 0;
    end
    
    function ZoomFn(hObject, eventdata)
        zoom = ~zoom;
        background_image = 0 * background_image;
        num_bg = 0;
    end
    
    function NormFn(hObject, eventdata)
        normalize = ~normalize;
        background_image = 0 * background_image;
        num_bg = 0;
    end
    
    function GainFn(hObject, eventdata)
        if hObject == Live.gaintext
            gain = int16(str2double(get(hObject,'String')));
        elseif hObject==Live.GainMinus
            gain = int16(gain / 2);
        elseif hObject==Live.GainPlus
            gain = int16(gain *2);
        end
        if gain>300,
            gain = int16(300);
        elseif gain <4,
            gain = int16(4);
        end
        
        % Must be between 4 and 300 (1000)
        mmc.setProperty(camera, 'Gain', num2str(gain))
        set(Live.gaintext, 'String',num2str(gain));
    end
    
    function ExpoFn(hObject, eventdata)
        if hObject == Live.expotext
            exposure = int16(str2double(get(hObject,'String')));
        elseif hObject==Live.ExpoMinus
            exposure = exposure - 10;
        elseif hObject==Live.ExpoPlus
            exposure = exposure + 10;
        end
        if exposure>9000, exposure = int16(9000);
        elseif exposure <1, exposure = int16(1); end
        mmc.setExposure(exposure);
        set(Live.expotext, 'String',num2str(round(exposure)));
    end
    
    function StartFn(hObject, eventdata)
        if ~mmc.isSequenceRunning()
            set(Live.GainMinus, 'enable', 'off');
            set(Live.GainPlus, 'enable', 'off');
            set(Live.ExpoMinus, 'enable', 'off');
            set(Live.ExpoPlus, 'enable', 'off');
            set(Live.gaintext, 'enable', 'off');
            set(Live.expotext, 'enable', 'off');
            set(Live.snap, 'enable', 'on');
            set(Live.focus, 'enable', 'off');
            
            camera_device = mmc.getCameraDevice();
            if strcmp('HamamatsuHam_DCAM',char(camera_device))
                mmc.setProperty(camera_device,'TRIGGER SOURCE', ...
                    'INTERNAL');
                real_bit = 16;
            elseif strcmp('Andor',char(camera_device))
                % These properties are set by 'Startup', which is loaded
                % whenever Live starts.
                %mmc.setProperty(camera_device,'Trigger','Internal');
                %mmc.setProperty(camera_device,'Shutter','Auto');
                mmc.setProperty('Arduino-Switch','Blanking Mode','Off');
                mmc.setProperty(camera_device,'FrameTransfer','On');
                real_bit = 14;
            end
            
            exitlive = false;
            mmc.startContinuousSequenceAcquisition(0);
            %             width = mmc.getImageWidth();
            %             height = mmc.getImageHeight();
            if mmc.getBytesPerPixel == 2
                pixelType = 'uint16';
            else
                pixelType = 'uint8';
            end
            
            while ~exitlive
                if (mmc.getRemainingImageCount() > 0)
                    %                     try
                    
                    imgtmp = mmc.getLastImage();
                    imgtmp = typecast(imgtmp, pixelType);
                    
                    img1 = double(rotateFrame(imgtmp,mmc));
                    img1 = img1 ./ (2^real_bit-1);
                    
                    if remove_background
                        stage_pos.x = mmc.getXPosition('XYStage');
                        stage_pos.y = mmc.getYPosition('XYStage');
                        stage_pos.z = mmc.getPosition('ZStage');
                        
                        if abs(stage_pos.z - prev_stage_pos.z) || ...
                                abs(stage_pos.x - prev_stage_pos.x) || ...
                                abs(stage_pos.y - prev_stage_pos.y)
                            background_image = 0 * background_image;
                            num_bg = 0;
                            prev_stage_pos = stage_pos;
                            stage_tic = tic;
                        elseif toc(stage_tic) < 0.5
                            background_image = 0 * background_image;
                            num_bg = 0;
                        end
                        
                        delta_background = img1 - background_image;
                        img1 = img1 - background_image + 0.5;
                        num_bg = num_bg + 1;
                        background_image = background_image + delta_background/num_bg;
                    end
                    
                    if width > display_size_width
                        if zoom,
                            img1 = imcrop(img1, ...
                                [floor(width/2)-display_size_width/2 ...
                                floor(height/2)-display_size_height/2 ...
                                display_size_width display_size_height]);
                        else
                            img1 = imresize(img1, ...
                                display_size_width/width, ...
                                'bicubic');
                        end
                    else
                        if zoom,
                            img1 = imresize(img1, ...
                                2*display_size_width/width, ...
                                'bicubic');
                            img1 = imcrop(img1, [display_size_width/2 ...
                                display_size_height/2 ...
                                display_size_width ...
                                display_size_height]);
                        end
                    end
                    
                    intensity = round(1000*std(img1(:)));
                    maxintensity = sum(round(1000*img1(:))== round(1000*max(img1(:))));
                    
                    set(Live.intensitytext, ...
                        'String',num2str(intensity));
                    set(Live.maxintensitytext, ...
                        'String',num2str(maxintensity));
                    
                    if normalize,
                        img1 = img1 - min(img1(:));
                        img1 = img1 ./ max(img1(:));
                    end
                    
                    img1(isnan(img1)) = 0;
                    
                    if crossmark,
                        img1(:,display_size_width/2) = ...
                            round(1-img1(:,display_size_width));
                        img1(display_size_height/2,:) = ...
                            round(1-img1(display_size_height,:));
                    end
                    
                    if circle,
                        for i = 1:numel(xc)
                            img1(xc(i),yc(i)) = round(1-img1(xc(i),yc(i)));
                        end
                    end
                    
                    Live.Image = imshow(img1,'parent',Live.Axes);
                    
                    %                     catch
                    %                     end
                end
                pause(0.0001);
            end
            mmc.stopSequenceAcquisition();
            mmc.waitForSystem();
        end
    end
    
    function StopFn(hObject, eventdata)
        exitlive = true;
        if gain >= 0
            set(Live.GainMinus, 'enable', 'on');
            set(Live.GainPlus, 'enable', 'on');
            set(Live.gaintext, 'enable', 'on');
        end
        set(Live.ExpoMinus, 'enable', 'on');
        set(Live.ExpoPlus, 'enable', 'on');
        set(Live.expotext, 'enable', 'on');
        %         set(Live.shutterDia,'Value', 0);
        %         set(Live.shutterEpi,'Value', 0);
        set(Live.snap, 'enable', 'off');
        set(Live.focus, 'enable', 'on');
        background_image = 0 * background_image;
        num_bg = 0;
    end
    
    function closeGUI(hObject, eventdata)
        exitlive = true;
        mmc.stopSequenceAcquisition();
        if arduino
            mmc.setProperty('Spectra','White_Enable',0);
            mmc.setProperty('Spectra','White_Level',0);
            mmc.setProperty('Spectra','State',0);
        else
            mmc.setProperty('TIEpiShutter','State',0);
        end
        mmc.setShutterOpen(0);
        mmc.waitForSystem();
        pause(0.005);
        delete(h);
    end
end
