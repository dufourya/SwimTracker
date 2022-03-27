function filterDetection(movie_name)
%%
file_info = RadFileInfo(movie_name);
metadata = Metadata(file_info.metaFile);
metadata.read();

params = [];
detected_objs_original = [];

detected_objs = [];
rad_detection = [];

if exist(strcat(file_info.detectionFile,'_original'),'file')
    load(strcat(file_info.detectionFile,'_original'),'-mat');
elseif exist(file_info.detectionFile,'file')
    load(file_info.detectionFile);
else
    error('no detection file for this movie.');
end

detected_objs_original = detected_objs;
%%
sc = get( groot, 'Screensize' );
f = figure('Position', [sc(1)+100 sc(2)+100 sc(3)-200 sc(4)-200]);
ax1 = axes('Parent',f,'position',[0 0.05 0.5 0.95]);

if exist(file_info.bin_file,'file')
    
    height = metadata.getImageHeightPixel();
    width = metadata.getImageWidthPixel();
    pixel_type = metadata.getPixelType();
    num_images = metadata.getNImagesBeforeFreeze();
    byteperpix = metadata.getBytesPerPixel();
    
    fid = fopen(file_info.bin_file);
    
    frames = zeros(width*height,5);
    ind = (1:5).*(floor(num_images/7));
    
    for i = 1:numel(ind)
        fseek(fid,ind(i)*width*height*byteperpix,'bof');
        frames(:,i) = fread(fid,width*height,pixel_type);
    end
    
    fclose('all');
    frame = double(rot90(reshape(frames(:,1),[width height])));
    scaling = [quantile(frames(:)/mean(frames(:)),0.0001) quantile(frames(:)/mean(frames(:)),0.9999)];
    img = frame./mean(frame(:));
    img = (img-scaling(1))/(scaling(2)-scaling(1));
    h1 = imshow(img, 'Parent',ax1);
    hold on
    ind_obj = detected_objs_original(5,:)==ind(1)+1;
    h0 = plot(detected_objs_original(1,ind_obj),detected_objs_original(2,ind_obj),'ro','Parent',ax1,'MarkerSize',10);
    h2 = plot(detected_objs_original(1,ind_obj),detected_objs_original(2,ind_obj),'go','Parent',ax1,'MarkerSize',10);
    
end
%%
ax2 = axes('Parent',f,'position',[0.6 0.25 0.35 0.6],'NextPlot','add');
histogram(log10(detected_objs_original(3,:)),'Normalization','pdf','Parent',ax2);

xlabel('Object mass');
ylabel('Probability density');
title('Masses of detected objects');

zeta = min(log10(detected_objs_original(3,:)));
ind_mass = log10(detected_objs_original(3,:))>=zeta;
h3 = plot(ax2,zeta,0,'d','MarkerSize',10,'MarkerFacecolor','r');

uicontrol('Parent',f,'Style','slider','Position',[10 40 (sc(3)-200)/2-20 10],...
    'value',1, 'min',1, 'max',5,'SliderStep',[1/5 1/5],'callback',{@SliderImg,h1,h2,h0});

a1 = uicontrol('Parent',f,'Style','text','Position',[(sc(3)-200)/4-10 10 40 20],'String',num2str(1));

uicontrol('Parent',f,'Style','slider','Position',[(sc(3)-200)/2+10 40 (sc(3)-200)/2-20 10],...
    'value',zeta, 'min',min(log10(detected_objs_original(3,:))), 'max',max(log10(detected_objs_original(3,:))),'callback',{@SliderFn,h2,h3});

uicontrol('Parent',f,'Style','text','Position',[(sc(3)-200)/2+10 10 40 20],'String',num2str(round(min(log10(detected_objs_original(3,:))),3,'significant')));
uicontrol('Parent',f,'Style','text','Position',[(sc(3)-200)-10 10 40 20],'String',num2str(round(max(log10(detected_objs_original(3,:))),3,'significant')));

b1 = uicontrol('Parent',f,'Style','text','Position',[3*(sc(3)-200)/4-10 10 40 20],'String',num2str(round(zeta,3,'significant')));


uicontrol('Parent',f,'Style','pushbutton','Position',[3*(sc(3)-200)/4-50 100 100 50],'String','Filter!','FontSize',20,'callback',@FilterFn);

    function SliderFn(hObject, eventdata, hplot1, hplot2)
        ind_mass = log10(detected_objs_original(3,:))>=hObject.Value;
        set(hplot1,'xdata',detected_objs_original(1,ind_obj & ind_mass));
        set(hplot1,'ydata',detected_objs_original(2,ind_obj & ind_mass));
        set(hplot2,'xdata',hObject.Value);
        set(b1,'String',num2str(round(hObject.Value,3,'significant')));
    end


    function SliderImg(hObject, eventdata, hplot1, hplot2, hplot3)
        newVal = round(hObject.Value);
        set(hObject,'Value',newVal);
        set(a1,'String',num2str(round(newVal,3,'significant')));
        frame = double(rot90(reshape(frames(:,newVal),[width height])));
        img = frame./mean(frame(:));
        img = (img-scaling(1))/(scaling(2)-scaling(1));
        set(hplot1,'CData',img);
        ind_obj = detected_objs_original(5,:)==ind(hObject.Value)+1;
        set(hplot3,'xdata',detected_objs_original(1,ind_obj));
        set(hplot3,'ydata',detected_objs_original(2,ind_obj));
        set(hplot2,'xdata',detected_objs_original(1,ind_obj & ind_mass));
        set(hplot2,'ydata',detected_objs_original(2,ind_obj & ind_mass));
    end

    function FilterFn(hObject, eventdata)
        movefile(file_info.detectionFile,strcat(file_info.detectionFile,'_original'));
        detected_objs = detected_objs_original(:,ind_mass);
        rad_detection = rad2utrack(detected_objs, metadata);
        save(file_info.detectionFile, 'detected_objs', 'rad_detection' , 'params');
    end

end