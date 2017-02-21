%% clearing up
clc;
clear all;
close all;
rng('shuffle')
format long

addpath('base/')

fprintf('################################################################\n')
fprintf('    This experiment uses genetic algorithm to optimize          \n')
fprintf('    the fuzzy logic system. Therefore, output values can be     \n')
fprintf('    slightly different in each run. But the observation         \n')
fprintf('    should remain similar.\n')
fprintf('################################################################\n')


%% Load data
data = csvread('data/quandl_sales_price.csv', 1, 0);

%% data split
training_parcent = 80;
testing_parcent = 20;

training_index = 1:floor(size(data,1)*training_parcent/100);
testing_index = training_index(end)+1:size(data,1);

training_data = data(training_index,:);
testing_data = data(testing_index,:);

%% FLS parameters

num_mf_input = [3 3 3 3 3 3];
num_mf_output = 3;
std_dev_UMF = {78*ones(5,1), 2*ones(3,1), 3*ones(3,1), 2*ones(3,1), 2*ones(3,1), 1.2*ones(3,1)};
std_dev_LMF = {75*ones(5,1), 1*ones(3,1), 2*ones(3,1), 1*ones(3,1), 1*ones(3,1), 0.5*ones(3,1)};
std_dev_UMF_output = {5*ones(3,1)};
std_dev_LMF_output = {2*ones(3,1)};

% generate rulebase
[ consolidated_rulebase, exhaustive_rulebase, exhaustive_rulebase_with_ruledegree,...
    center_mf_input,center_mf_output, num_input, frequency] = ...
    generate_WM_rulebase_IT2_differentInputMFs_dfrntStdDevs(training_data, training_data, ...
    num_mf_input, std_dev_UMF,std_dev_LMF,...
    std_dev_UMF_output,...
    std_dev_LMF_output, num_mf_output  );

% optimize
[ center_input_MF_tuned, center_output_MF_tuned, std_dev_LMF_tuned, std_dev_UMF_tuned ] = ...
    fitness_IT2FLS_centerOptimzd...
    ( num_mf_input, num_mf_output, consolidated_rulebase, training_data, testing_data );

plot_IT2MFs_new( center_input_MF_tuned, std_dev_UMF_tuned, std_dev_LMF_tuned )

%% with optimized FLS
defuzzified_output_tuned_IT2 = [];

weight_matrix = ones(1, size(consolidated_rulebase,1));

for k = 1: size(testing_data,1)
    new_input = testing_data(k,1:end-1);
    
    [ defuzzified_output_tuned_IT2(k), yl(k), yr(k) ] = ...
    FLS_output_WM_IT2_withWeights_differentInputMFs_dfrntStdDevs( ...
                                            new_input, consolidated_rulebase,...
                                            center_input_MF_tuned,...
                                            center_output_MF_tuned, ...
                                            std_dev_UMF_tuned, std_dev_LMF_tuned, weight_matrix);
     
end

% performance
MAPE_tuned_IT2 = mean(abs((defuzzified_output_tuned_IT2'-testing_data(:,end))./testing_data(:,end)))*100;
RMSE_tuned_IT2 = mean((defuzzified_output_tuned_IT2'-testing_data(:,end)).^2);

%% Get the maximum of 'fatness' of MF for each inputs for each test sample

for k = 1: size(testing_data,1) % for all test samples
    new_input = testing_data(k,1:end-1);
    for i = 1:length(new_input) % for all inputs in the observation
        tmp_centerMF = center_input_MF_tuned{i};
        tmp_stdDevs_UMF = std_dev_UMF_tuned{i};
        tmp_stdDevs_LMF = std_dev_LMF_tuned{i};
        for j = 1:length(tmp_centerMF) % 
            tmp_1(j) = compute_membership_grade_new(new_input(i), tmp_centerMF(j), tmp_stdDevs_UMF(j));
            tmp_2(j) = compute_membership_grade_new(new_input(i), tmp_centerMF(j), tmp_stdDevs_LMF(j));
        end
        
        fatness_allinputs(i) = sum(tmp_1-tmp_2);
        tmp_2 = []; tmp_1 = [];
       
    end
    total_fatness_each_sample(k) = sum(fatness_allinputs);
    fatness_allinputs = [];
end
%  Get the uncertainty score
SGL = min(data(:,end));
SGR = max(data(:,end));
U = (yr-yl)/(SGR-SGL);

%% Plots
figure; 
plot(defuzzified_output_tuned_IT2,'g'); 
hold on; 
plot(testing_data(:,end),'r'); 
hold on
plot(yl, 'k')
hold on
plot(yr,'k')
title('IT2FLS sales forecast')
legend('Prediction', 'Actual', 'left shoulder - centroid', 'right shoulder - centroid')

figure;
X=[1:length(yl),fliplr(1:length(yl))];                % create continuous x value array for plotting
Y=[yl,fliplr(yr)];                                    % create y values for out and then back
fill(X,Y,[0.248 0.248 0.248]);                        % plot filled area
hold on
plot(defuzzified_output_tuned_IT2,'g'); 
hold on; 
plot(testing_data(:,end),'r'); 
title('IT2FLS sales forecast')
legend('Centroid','Prediction', 'Actual')

figure;
plot(defuzzified_output_tuned_IT2)
hold on; 
plot(testing_data(:,end))
yyaxis right
plot(U*100,'c')
legend('Forecast', 'Actual', 'Uncertainy score (%)')


%  Fatness vs. uncertainty score
figure;
plot(total_fatness_each_sample,'r')
hold on
plot(U,'b')
hold off
xlabel('Test set ID')
ylabel('Sum of input uncertainty, Uncertainty score ')
legend('Uncertainty about input membership','Uncertainty score')

% Scatter plot of fatness and uncertainty score
figure; 
scatter(total_fatness_each_sample, U)

% prediction error vs uncertainty score
pe = defuzzified_output_tuned_IT2'-testing_data(:,end);
figure;
scatter(pe, U)
xlabel('Absolute Prediction Error')
ylabel('Uncertainty Score')


% proving that the max of forecast and min of U is not related, rebuttal
figure; 
plot(abs(pe)); 
hold on; 
plot(U)
xlabel('Test instances')
ylabel('Prediction error, uncertainty score')
legend('Absolute prediction error', 'Uncertainty score')




%% fatness vs uncertainty score explanation
% Extract the fatness of inputs in observations where uncertainty score is
% very high and very low

tt = sort(U);

investigation_index = [find(U==tt(1)) find(U==tt(2)) find(U==tt(3))...
    find(U==tt(end-2)) find(U==tt(end-1)) find(U==tt(end))];

absolute_error = abs(defuzzified_output_tuned_IT2'-testing_data(:,end));
fatness_investigationData = [];
fatness_investigation = [];
for k = 1: size(testing_data(investigation_index,:),1) % for all investigation samples
    new_input = testing_data(investigation_index(k),1:end-1);
    for i = 1:length(new_input) % for all inputs in the observation
        tmp_centerMF = center_input_MF_tuned{i};
        tmp_stdDevs_UMF = std_dev_UMF_tuned{i};
        tmp_stdDevs_LMF = std_dev_LMF_tuned{i};
        for j = 1:length(tmp_centerMF) % 
            tmp_1(j) = compute_membership_grade_new(new_input(i), tmp_centerMF(j), tmp_stdDevs_UMF(j));
            tmp_2(j) = compute_membership_grade_new(new_input(i), tmp_centerMF(j), tmp_stdDevs_LMF(j));
        end
        
        fatness_investigationData(i) = sum(tmp_1-tmp_2);
        tmp_2 = []; tmp_1 = [];
       
    end
  
    % format: [sum_uncertainty_per_input_all_MF(6 input in total)
    % total_uncertainty_per_observation Uncertainty_score_per_observation
    % abs_error_per_observation]
    fatness_investigation = vertcat(fatness_investigation, ...
        [investigation_index(k) fatness_investigationData sum(fatness_investigationData) U(investigation_index(k)) absolute_error(investigation_index(k))]);
    fatness_investigationData = [];
end
myheader= {'Index','I1','I2','I3','I4','I5','I6','sumU', 'UScore', 'abs_error'};
disp(array2table(fatness_investigation, 'VariableNames', myheader))
