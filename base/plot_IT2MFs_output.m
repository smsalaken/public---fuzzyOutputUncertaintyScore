function plot_IT2MFs_output( center_output_MF, std_dev_UMF_output, std_dev_LMF_output )
%UNTITLED4 plots membership functions of interval typ-2 FS
%   input:
%   center_output_MF: contains center of guassian MFs on each output. Centers
%   for one output variable is expected in one column, For example, if there
%   is 3 output in FLS, and each output has 4 MFs on them, then the size of 
%   center_output_MF is 4x3.
%   std_dev_UMF_output: standard deviation of upper MF
%   std_dev_LMF_output: stndard deviation of lower MF
tmp = center_output_MF;
center_output_MF = [];
center_output_MF = tmp';
for i = 1:size(center_output_MF,2)
    figure;
    tmp_cntr = center_output_MF(:,i);
    for j = 1: length(tmp_cntr)
        x = min(tmp_cntr):0.01:max(tmp_cntr);
        y1 = gaussmf(x,[std_dev_UMF_output tmp_cntr(j)]);
        plot(x,y1);
        hold on
        y2 = gaussmf(x,[std_dev_LMF_output tmp_cntr(j)]);
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

