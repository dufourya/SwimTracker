function fignum = setupFigure(figname,colors)

all_fig_handles = get(0,'children');
all_fig_names   = get(all_fig_handles, 'name');

if any(strcmp(all_fig_names, figname))
    
    figs_to_close = all_fig_handles(strcmp(all_fig_names, figname));
    close(figs_to_close)
    
end

fignum = figure;
set(fignum,'name',figname)
set(fignum,'color','w')
set(fignum,'position', [2400 100 800 800])
set(fignum,'DefaultAxesLinewidth',2)
set(fignum,'DefaultAxesTickdir','out')
set(fignum,'DefaultAxesticklength',[0.03 0.02])
set(fignum,'DefaultAxesFontSize',18)
set(fignum,'DefaultLineLinewidth',3)
set(fignum,'DefaultSurfaceEdgeColor','none')

if nargin > 1
    
    set(fignum,'DefaultAxesColorOrder',colors)
    
end


end