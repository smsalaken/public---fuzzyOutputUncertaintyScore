function [ consolidated_rulebase, exhaustive_rulebase, exhaustive_rulebase_with_ruledegree,...
    center_mf,center_mf_output, num_input, frequency] = ...
    generate_WM_rulebase_IT2_differentInputMFs_dfrntStdDevs( consildated, trn_data, num_mf_array, std_dev_UMF, ...
    std_dev_LMF, std_dev_output_UMF, std_dev_output_LMF, num_mf_output )

%% generate_WM_rulebase_IT2   
% Summary: Generate Wang-Mendel rulebase for ONE
%          output type-2 FLS system


%   This function generates fuzzy rulebase for interval type-2 fuzzy sets based on
%   numerical data. It employes Gaussian MFs on input and output to
%   generate rulbase. IT USES DIFFERENT STD DEVs FOR "EACH MF" ON EVERY INPUT
%   AND OUTPUT

% Description of parameters:

% Inputs:
% consildated: matrix in the form [training_data checking_data] to determin
% maximum and minimum of inputs
% trn_data: matrix containing training samples, last column corresponds to output
% num_mf_array: number of MF on each input (defines number of FUZZY category on
% INPUT e.g. small, average, big), Assumption: all input has same number of
% MFs
% std_dev_UMF: standard deviation for input UMF (LMF is not needed for 
% rulebase construction as highest value of UMF will determine the fuzzy 
% category of input sample)
% num_mf_array_output: Number of MFs on output domain (defines number of FUZZY 
% category on OUTPUT e.g. small, average, big)

% Outputs:
% consolidated_rulebase: consoldated rulebase which contains only one rule
% from each conflict group (e.g. same antecedent with different consequent)
% based on the maximum degree of rule
% consolidated_rulebase_with_degree: consolidated rulebase with associated
% degree of rules
% rulebase: exhaustive rulebase, without any kind of modification and
% conflict removal
% center_mf = center of MFs on inputs, required to compute membership grade
% on new inputs (e.g. testing data)
% center_mf_output= center of MFs on output, required to compute wang
% mendel defuzzified output value



% Reference: Wang, L-X., and Jerry M. Mendel. "Generating fuzzy rules by 
% learning from examples." Systems, Man and Cybernetics, IEEE Transactions
% on 22.6 (1992): 1414-1427.
	

% Author: SM Salaken, March 31, 2015


%% Find input ranges

[num_samples, a] = size(trn_data);
num_input = a-1; % last columns corresponds to output samples
for i = 1:num_input
    min_inputs(i) = min(consildated(:,i));
    max_inputs(i) = max(consildated(:,i));
end

%% INPUT PORTION
% generate center of MF in the input range
for i = 1:num_input
    center_mf {i} = linspace(min_inputs(i)+0.25*(min_inputs(i)),max_inputs(i)-0.25*(min_inputs(i)),num_mf_array(i));
end


%  find MF of each sample of each input
for i = 1:num_samples
    for j = 1:num_input
        b = [];
        tmp_center = [];
        tmp_center = center_mf{j};
        tmp_std_dev_UMF = std_dev_UMF{j};
        tmp_std_dev_LMF = std_dev_LMF{j};
        for k = 1: length(tmp_center)
%             b(k) = exp(-(trn_data(i,j)-center_mf(k,j))^2/(2*std_dev^2));
            b1(k) = compute_membership_grade_new(trn_data(i,j), tmp_center(k), tmp_std_dev_UMF(k));
            b2(k) = compute_membership_grade_new(trn_data(i,j), tmp_center(k), tmp_std_dev_LMF(k));
            b(k) = (b1(k)+ b2(k))/2;
        end
        [mf(i,j), mf_cat(i,j)] = max(b);
    end
end


%% OUTPUT PORTION

center_mf_output = linspace(min(trn_data(:,end)),max(trn_data(:,end)),num_mf_output);

% find MF of each sample of output - ONE output system
z = trn_data(:,end);
t_std_dev_output_UMF = std_dev_output_UMF{1};
t_std_dev_output_LMF = std_dev_output_LMF{1};
for i = 1:length(z)
    for k = 1: length(center_mf_output)
            b_o1(k) = compute_membership_grade_new(z(i), center_mf_output(k), t_std_dev_output_UMF(k));
            b_o2(k) = compute_membership_grade_new(z(i), center_mf_output(k), t_std_dev_output_LMF(k));
            b_o(k) = (b_o1(k)+b_o2(k))/2;
    end
    [mf_o(i), mf_cat_o(i)] = max(b_o);
end

% interim rulebase - whole numbers represent fuzzy category of input
% numbers
rulebase = [];
for i = 1:num_samples
    rulebase = [rulebase; [mf_cat(i,:) trn_data(i,end)]];
end

%% exhaustive rulebase - final
% interim exhaustive rulebase, all numbers represent fuzzy category of input 
% and output
% This is EXHAUSTIVE RULEBASE
rulebase = [rulebase(:,1:end-1) mf_cat_o'];

% memerbship degree associated with each input and output of each rule
rulebase_mf = [mf mf_o']; % [MF_inputs MF_output]


%% reducing rulebase size

% find degree of each rule
degree_of_rule = prod(rulebase_mf,2);
     
% exhaustive rulebase with degree of rules
exhaustive_rulebase_with_ruledegree = [rulebase degree_of_rule];


exhaustive_rulebase = rulebase;
[consolidated_rulebase, frequency] = rulebase_pruner( exhaustive_rulebase, exhaustive_rulebase_with_ruledegree );



end

