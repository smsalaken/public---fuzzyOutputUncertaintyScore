function  membership_grade  = compute_membership_grade_new( x, c, std_dev )
%compute_membership_grade Compute membership degree
%   This function compute membership degree/grade of given input based on a
%   Gaussian MF. x is the data point for which memebrship degree is
%   required, c is the center of Gaussian MF and std_dev is the standard
%   deviation of that Gaussian MF
        membership_grade = exp(-(x-c)^2/(2*std_dev^2));


end

