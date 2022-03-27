function plotSwimTrackerTables(data,names,figtitle)

picformat = 'svg';

if isnumeric(names)
    names = cellstr(num2str(names));
end

unames = unique(names);
c = lines(numel(unames));

screensize = get( groot, 'Screensize' );

h = figure('Name',figtitle,'Position',[round((screensize(3)-1600)/2)...
    round((screensize(4)-1*300)/2) 1600 1*300]);

for i = 1:numel(unames)
    
    ind = contains(names,unames(i));
    temp = data(ind,:);
    
    subplot(1,4,1), hold on,
    plot(i,i,'LineWidth',2,'Color',c(i,:));
    axis off
    axis square
    
    subplot(1,4,2), hold on,
    w = 6/30;
    edges = -10:w:10;
    counts = zeros(numel(edges)-1,1);
    for k = 1:numel(counts)
        ind = log10(temp.diffcoeff_cve_mean)>=edges(k) & log10(temp.diffcoeff_cve_mean)<edges(k+1);
        counts(k) = sum(temp.trajtime(ind));
    end
    counts = counts/sum(temp.trajtime)/w;
    stairs(edges(1:end-1),counts,'LineWidth',2,'Color',c(i,:));
    xlim([-2 4]);
    ax = gca;
    ax.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8];
    ax.XTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'...
        '10^3' '10^4' '10^5' '10^6' '10^7' '10^8'};
    ax.FontSize = 14;
    xlabel('Diffusion coefficient (\mum^2/s)','FontSize',16);
    ylabel('Probability density','FontSize',16);
    axis square;
    
    subplot(1,4,3), hold on,
    w = 50/30;
    edges = 0:w:50;
    counts = zeros(numel(edges)-1,1);
    for k = 1:numel(counts)
        ind = temp.meanspeed>=edges(k) & temp.meanspeed<edges(k+1);
        counts(k) = sum(temp.trajtime(ind));
    end
    counts = counts/sum(temp.trajtime)/w;
    stairs(edges(1:end-1),counts,'LineWidth',2,'Color',c(i,:));
    xlim([0 50]);
    ax = gca;
    ax.FontSize = 14;
    xlabel('Swimming speed (\mum/s)','FontSize',16);
    ylabel('Probability density','FontSize',16);
    axis square;
    
    subplot(1,4,4), hold on,
    w = 3/30;
    edges = -10:w:10;
    counts = zeros(numel(edges)-1,1);
    for k = 1:numel(counts)
        ind = log10(temp.diffcoeff_cve_runtime)>=edges(k) & log10(temp.diffcoeff_cve_runtime)<edges(k+1);
        counts(k) = sum(temp.trajtime(ind));
    end
    counts = counts/sum(temp.trajtime)/w;
    stairs(edges(1:end-1),counts,'LineWidth',2,'Color',c(i,:));
    xlim([-2 1]);
    ax = gca;
    ax.FontSize = 14;
    ax.XTick = [-5 -4 -3 -2 -1 0 1 2 3 4 5];
    ax.XTickLabel = {'10^{-5}' '10^{-4}' '10^{-3}' '10^{-2}' '10^{-1}' '10^0' '10^1' '10^2'...
        '10^3' '10^4' '10^5'};
    xlabel('Dir. persistence time (s)','FontSize',16);
    ylabel('Probability density','FontSize',16);
    axis square;
    
end

subplot(1,4,1),
legend(unames,'Location','northwest','Interpreter','none');
% 
% set(h,'PaperPositionMode','auto');
% print(h,strcat(figtitle,sprintf('.%s',picformat)),sprintf('-d%s',picformat),'-r300','-noui');
% print(h,strcat(figtitle,sprintf('.%s','png')),sprintf('-d%s','png'),'-r300','-noui');
% crop(strcat(figtitle,sprintf('.%s','png')));

