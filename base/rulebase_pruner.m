function [ consolidated_rulebase, frequency ] = rulebase_pruner( rulebase, exhaustive_rulebase_with_ruledegree )
%UNTITLED3 find the consolidated database by removing duplicates from
%exhaustive rulebase based on rules with highest degree
%   exhaustive_rulebase_with_ruledegree: contains degree of each rule on
%   the last column
%   rulebase: exhaustive rulebase

exhaustive_rulebase = rulebase;
a = exhaustive_rulebase(:,1:end-1);
c = [exhaustive_rulebase_with_ruledegree(:,1:end-2) exhaustive_rulebase_with_ruledegree(:,end)];
d = exhaustive_rulebase_with_ruledegree(:,end-1); % consequent part only
frequency = [];
i = 1;
consolidated_rulebase_t = [];
kept_rule = [];

while size(a,1)
    counter = 0;
    temp = a(1,:);
    index_vector = zeros(size(a,1),1);
    for j = 1:size(a,1)
        
        if (temp == (a(j,:)))
            counter = counter + 1;
            index_vector (j) = 1;
            
        end
    end
    
    b = c(logical(index_vector),:);
    [~,rule_to_keep]= max(b(:,end));
    kept_rule = [kept_rule; b(rule_to_keep,:) d(rule_to_keep)];
    
    a(logical(index_vector),:) = [];
    c(logical(index_vector),:) = [];
    d(logical(index_vector),:) = [];
    consolidated_rulebase_t = [consolidated_rulebase_t; temp];
    
    frequency(i) = counter;
    i = i+1;
    
end
consolidated_rulebase = [kept_rule(:,1:end-2) kept_rule(:,end)];


end

