function [rownorm columnnorm] = conf_count_plot(Q,font_size,figtitle)

%Q is the confusion matrix of raw counts
%rownormalized =1 for using row normalized, default is column
%Add axes labels

rawcounts = Q;
counts=rawcounts;
counts(size(counts,1),:)=[];
counts(:,1)=[];
rownorm=zeros(size(counts));
columnnorm=rownorm;

for row=1:size(counts,1)
    rownorm(row,:)=counts(row,:)/sum(counts(row,:));
end
rownorm(end)=NaN;

for column=1:size(counts,2)
    columnnorm(:,column)=counts(:,column)/sum(counts(:,column));
end
columnnorm(end)=NaN;

if rownorm==1
    normconfmatrix=rownorm;
    disp('row normalized');
else
    normconfmatrix=columnnorm;
    disp('column normalized');
end

xlabels=rawcounts(size(rawcounts,1),:);
xlabels(1)=[];
xlabels(end)=[];
xlabels=xlabels';
ylabels=rawcounts(:,1);
ylabels(end)=[];
ylabels(end)=[];
text_counts=counts(find(counts));
[ytext xtext]=find(counts);

[nr,nc] = size(normconfmatrix);
shiftedXtick=[1:size(counts,2)]+.5;
shiftedYtick=[1:size(counts,1)]+.5;

figure
subplot(2,1,1)
set(gca,'XTick',[],'YTick',[])
%axis square
dataAxes = axes;
set(gca,'XTick',[],'YTick',[])
invisibleAxes = axes('Visible','off','Position',get(dataAxes,'Position'),'HitTest','off','XLimMode','manual','YLimMode','manual');
linkaxes([dataAxes invisibleAxes])
pcolor([normconfmatrix nan(nr,1); nan(1,nc+1)])
shading flat
%axis square
set(gca,'XTick',shiftedXtick,'YTick',shiftedYtick, 'YTickLabel',{num2str(ylabels),' '},'XTickLabel',{num2str(xlabels),' '}, 'FontSize', font_size, 'ydir', 'reverse')


for ii=1:length(find(counts))
    icount=find(counts);
    b=normconfmatrix(icount(ii));
    color=1-((b<0.5)+((b-0.5)/0.5)*0.25);
    text('Parent',invisibleAxes,'String',num2str(text_counts(ii)),'HorizontalAlignment','Center','VerticalAlignment','Middle','Position',[xtext(ii)+.5 ytext(ii)+0.5],'FontSize',font_size,'Color',[color color color])
end

line([(size(counts,2)) (size(counts,2))],[0 size(counts,1)],'Parent',invisibleAxes,'Color','k','LineWidth',2);
line([0 (size(counts,2))],[size(counts,1) size(counts,1)],'Parent',invisibleAxes,'Color','k','LineWidth',2);
line([1 1],[0 size(counts,1)],'Parent',invisibleAxes,'Color','k','LineWidth',2);
line([0 (size(counts,2))],[1 1],'Parent',invisibleAxes,'Color','k','LineWidth',2);

colormap(flipud(gray))
%colorbar
title(figtitle)

set(gcf,'Position',[100,100,1000,1000]);