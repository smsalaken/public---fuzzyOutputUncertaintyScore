function plot_IT2MFs( center_input_MF, std_dev_UMF, std_dev_LMF )
%UNTITLED4 plots membership functions of interval typ-2 FS
%   input:
%   center_input_MF: contains center of guassian MFs on each input. Centers
%   for one inpu variable is expected in one column, For example, if there
%   is 3 input in FLS, and each input has 4 MFs on them, then the size of 
%   center_input_MF is 4x3.
%   std_dev_UMF: standard deviation of upper MF
%   std_dev_LMF: stndard deviation of lower MF

for i = 1:length(center_input_MF)
    figure;
    tmp_cntr = center_input_MF{i};
    for j = 1: length(tmp_cntr)
    x = min(tmp_cntr):1:max(tmp_cntr);
    y1 = gaussmf(x,[std_dev_UMF(i) tmp_cntr(j)]);
    plot(x,y1);
    hold on
    y2 = gaussmf(x,[std_dev_LMF(i) tmp_cntr(j)]);
    plot(x,y2);
    hold on
    X=[x,fliplr(x)];                %#create continuous x value array for plotting
    Y=[y1,fliplr(y2)];              %#create y values for out and then back
    fill(X,Y,[0.248 0.248 0.248]);                  %#plot filled area
    end
    % title(sprintf('IT2 MFs on input-%d',i))
    hold off
    set(gca,'FontSize',14)
    set(findall(gcf,'type','text'),'FontSize',14)
%     savefig(sprintf('inputMF/IT2/IT2-%d.fig',i))
%     print(sprintf('inputMF/IT2/IT2-%d',i),'-depsc','-tiff')
end
hold off

% title('MFs on each input')

end

