function plot_type1OutputMFs( center_output_MF,  std_dev_output)
%UNTITLED5 plots output membership functions of type-1 fuzzy set
%   center_output_MF: contains center of output MFs
%   std_dev_output: standard deviation of output MFs
figure;
for i = 1:length(center_output_MF)
    plot(min(min(center_output_MF)):0.01:max(max(center_output_MF)),gaussmf(min(min(center_output_MF)):0.01:max(max(center_output_MF)),[std_dev_output center_output_MF(i)]));
    hold on
end
hold off
title('MFs on output')

end

