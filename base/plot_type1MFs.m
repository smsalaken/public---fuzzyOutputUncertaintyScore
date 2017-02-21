function plot_type1MFs( center_input_MF, std_dev )
%UNTITLED3 plots input membership functions of type-1 fuzzy set
%   input:
%   center_input_MF: contains center of guassian MFs on each input. Centers
%   for one inpu variable is expected in one column, For example, if there
%   is 3 input in FLS, and each input has 4 MFs on them, then the size of 
%   center_input_MF is 4x3.
%   std_dev_UMF: standard deviation of membership functions, all MFs has
%   same std_dev in this implementation
center_input_MF

for i = 1:length(center_input_MF)
    figure;
    tmp_cntr = center_input_MF{i};
    for j = 1: length(tmp_cntr)
    plot(min(tmp_cntr):0.01:max(tmp_cntr),gaussmf(min(tmp_cntr):0.01:max(tmp_cntr),[std_dev(i) tmp_cntr(j)]));
    hold on
    end
    %title(sprintf('T1 MFs on input-%d',i))
    hold off
    set(gca,'FontSize',14)
    set(findall(gcf,'type','text'),'FontSize',14)
%     savefig(sprintf('inputMF/type-1/T1-%d.fig',i))
%     print(sprintf('inputMF/type-1/T1-%d',i),'-depsc','-tiff')
end
hold off
% title('MFs on each input')
end

