function [ center_input_MF_tuned, center_output_MF_tuned,...
    std_dev_LMF_tuned, std_dev_UMF_tuned ] = ...
    fitness_IT2FLS_centerOptimzd( num_mf_input, num_mf_output, rulebase,...
                                data_Trn, testing_data )

%fitness_IT2FLS_centerOptimzd Optimize IT2 FLS
%   return the optimized center and standard deviation of Membership
%   functions

% init
input_Trn = data_Trn(:, end-1);
output_Trn = data_Trn(:, end);
num_mf = [num_mf_input num_mf_output];

% how many variable to optimize?
num_variable_to_optimize = 3*sum(num_mf_input)+num_mf_output; % center of mf + std dev LMF + std dev UMF for all input and output
                                                              % No need to
                                                              % optimize
                                                              % output STD
                                                              % devs


% bounds
zz = [];
for q = 1:size(data_Trn,2)
    zz = [zz min(data_Trn(:,q))*ones(1, num_mf(q))];
end
LB = [zz 0.01*ones(1,2*sum(num_mf_input))];

zz = [];
for q = 1:size(data_Trn,2)
    zz = [zz max(data_Trn(:,q))*ones(1, num_mf(q))];
end
UB = [zz 100*ones(1,2*sum(num_mf_input))];

% initial population for genetic algorithm
Initial_pop = [2007 2011 2014 3 6 9 5 14 25 3 6 12 3 6 11 3 7 11 3 8 11 ...
    60*ones(1,num_mf(1)) 40*ones(1,num_mf(1)) 5*ones(1,5*num_mf(2)) 2*ones(1,5*num_mf(2))];

% Start with the default options
options = gaoptimset;

% Modify options setting
options = gaoptimset(options,'Display', 'iter');
options = gaoptimset(options,'PlotFcns', { @gaplotbestf });
options = gaoptimset(options, 'UseParallel', true);
options = gaoptimset(options,'PopulationSize', 250);
options = gaoptimset(options,'Generations', 150);
options = gaoptimset(options, 'InitialPopulation', Initial_pop);
% options = gaoptimset(options,'CrossoverFraction', 0.5);
% options = gaoptimset(options,'MigrationFraction', 0.7);
% options = gaoptimset(options,'TimeLimit', 15*60); % run for 5 CPU minutes

% call GA
[final_param, fval] = ga(@nested_fitness_function, ...
    num_variable_to_optimize,[],[],[],[],LB,UB,[],[],options);

% final parameters
% extract variables from final paramter list
center_MF_t = [];
r = 1;
for l = 1:length(num_mf)
    center_MF_t{l} = final_param(r:(r+num_mf(l)-1));  % center of MF for input and output
    r = num_mf(l)+1;
end
center_input_MF_tuned = center_MF_t(1:end-1);
center_output_MF_tuned = center_MF_t{end};

Temp = final_param(sum(num_mf)+1:end);
for AA = 1:length(num_mf)-1
    Temp2 = Temp(1:2*num_mf(AA));
    Temp3 = reshape(Temp2, num_mf(AA),2);
    std_dev_LMF_tuned{AA} = min(Temp3');
    std_dev_UMF_tuned{AA} = max(Temp3');
    Temp(1:2*num_mf(AA)) = [];
end




%=========================================================================
% Nested cost function
    function cost = nested_fitness_function(x)
       
       center_MF = []; 
       p = 1;
       for w = 1:length(num_mf)
            center_MF{w} = x(p:(p+num_mf(w)-1));  % center of MF for input and output
            p = num_mf(w)+1;
        end
        center_input_MF = center_MF(1:end-1);
        center_output_MF = center_MF{end};
        
        temp = x(sum(num_mf)+1:end);
        for aa = 1:length(num_mf)-1
            temp2 = temp(1:2*num_mf(aa));
            temp3 = reshape(temp2, num_mf(aa),2);
            std_dev_LMF{aa} = min(temp3');
            std_dev_UMF{aa} = max(temp3');
             temp(1:2*num_mf(aa)) = [];
        end

            
        
        weight_matrix = ones(1, size(rulebase,1));
        
        defuzzified_output_IT2 = zeros(1,size(testing_data,1));
        for k = 1: size(testing_data,1)
            new_input = testing_data(k,1:end-1);

            % call the fuzzy evaluation function for new input
            [ defuzzified_output_IT2(k), yl(k), yr(k) ] = ...
                FLS_output_WM_IT2_withWeights_differentInputMFs_dfrntStdDevs( new_input, rulebase, ...
                center_input_MF, center_output_MF, ...
                std_dev_UMF, std_dev_LMF, weight_matrix);
        end
        
        lambda = 1;
        cost = sqrt(mean((defuzzified_output_IT2'-testing_data(:,end)).^2)) + ...
            0.5*(1/numel(std_dev_LMF))*lambda*mean((cellfun(@mean, std_dev_UMF)-cellfun(@mean, std_dev_LMF)).^2);
        
        
    end % nested cost function

end % main function

