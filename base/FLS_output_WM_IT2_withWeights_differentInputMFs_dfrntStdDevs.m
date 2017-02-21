function [ defuzzified_output, left_end, right_end ] = ...
    FLS_output_WM_IT2_withWeights_differentInputMFs_dfrntStdDevs( new_input, consolidated_rulebase, ...
    center_input_MF, center_output_MF, ...
    std_dev_UMF, std_dev_LMF, weight_matrix)
%FLS_output_WM_type1 :Return the defuzzified output of FLS
%   This function returns the output of FLS based on given rulebase and
%   related information

% Input: 
% new_input: vector containing one instance of every input variable
% cconsolidated_rulebase: rulebase of FLS
% std_dev_UMF: a cell array. Contains std dev for UMF of every MF on every
% input e.g. is the system has 3 inputs and number of MF = [3 2 5] then it
% is a 3 element cell array which contain THREE matrix of dimension [1x3],
% [1x2] and [1x5]
% std_dev_LMF : similar as above

format long


% set weights to one if not supplied
if nargin<8
    weight_matrix = ones(size(consolidated_rulebase,1),1);
end
num_input = size(new_input,2);

    combined_mg_new_input = [];
    corresponding_output_center = [];
    applicability_of_rule = [];
    for i = 1: size(consolidated_rulebase,1) % loop through every rule

        for j = 1: num_input % find membership grade for all input in every rule
            tmmmpp_centerInputMF = [];
            system_input =  consolidated_rulebase(i,j); % center of corresponding rule
            tmmmpp_centerInputMF = center_input_MF{j};
            tmp_std_dev_LMF = std_dev_LMF{j};
            tmp_std_dev_UMF = std_dev_UMF{j};
            LMF_membersip_grade_new_input(j) = compute_membership_grade_new(new_input(j), tmmmpp_centerInputMF(system_input), tmp_std_dev_LMF(system_input));
            UMF_membersip_grade_new_input(j) = compute_membership_grade_new(new_input(j), tmmmpp_centerInputMF(system_input), tmp_std_dev_UMF(system_input));
        end
%         combined_mg_new_input(i) = prod(membersip_grade_new_input);
        corresponding_output_center(i) = center_output_MF(consolidated_rulebase(i, end));
%         applicability_of_rule(i) = combined_mg_new_input(i)*corresponding_output_center(i);
        
        firing_strength_LMF(i) = prod(LMF_membersip_grade_new_input)*weight_matrix(i);
        firing_strength_UMF(i) = prod(UMF_membersip_grade_new_input)*weight_matrix(i);
            
    end
    
    rule_consequent = corresponding_output_center;
%     whos
    [defuzzified_output, left_end, right_end] = EIASC(rule_consequent,[],firing_strength_LMF, firing_strength_UMF, +1);
end

