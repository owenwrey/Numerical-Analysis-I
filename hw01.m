% Author: Your Name / your_email
% Date: 2024-09-01
% Assignment Name: hw01

classdef hw01
    methods (Static)

        function p1()
            % This function only contains comments. Fill the following table. Do not write any code here.
            % :return: no returns

            % Write your result and explanation for each command here.
            % 
            % commands         |  results      | explanations
            % -----------------|---------------|-----------------------------------
            % eps              |2.2204e-16     |distance to next FP value
            % realmax          |1.7977e+308    |largest FP value in (2^1024)
            % realmin          |2.2251e-308    |smallest normalized zon-zero FP value (2^-1022)
            % 1 + eps - 1      |2.2204e-16     |no error
            % 1 + eps/2 - 1    |0              |since 1+eps/2 doesn't add to the next FPV (the next FPV is 1 + eps), it outputs zero (effectively doing 1-1) 
            % realmin/1e10     |2.2251e-318    |returns a non-normalized number, not underflow
            % realmin/1e16     |0              |underflow error, returns zero
            % realmax*10       |Inf            |overflow error
        end

        function s_n = p2(n, choice)
            % This function computes the Archimedes' method for pi.
            % :param n: the number of sides of the polygon
            % :param choice: 1 or 2, the formula to use
            % :return: s_n, the approximation of pi using Archimedes' method.

            % Tabulate the error of |s_n - pi| for n = 0, 1, 2, ..., 15 and choice = 1 and 2.
            % for both choices of formulas.
            % n     | choice 1 | choice 2
            % ------|----------|----------
            % 0     |          |
            % 1     |3.2154    |
            % 2     |3.1597    |
            % 3     |3.1461    |
            % 4     |3.1427    |
            % 5     |3.1419    |
            % 6     |3.1417    |
            % 7     |3.1416    |
            % 8     |3.1416    |
            % 9     |3.1416    | 
            % 10    |3.1416    |
            % 11    |3.1416    |
            % 12    |3.1416    |
            % 13    |3.1416    |
            % 14    |3.1416    |
            % 15    |3.1416    |

            % Explanation of the results (why there is a difference between the two choices):
            %
            %
            %
            %
            %
            %
            
            % Write your code here
            p_n = 1/sqrt(3);
            if choice == 1
                for i = 0:(n-1)
                    p_n = (sqrt(1 + p_n^2) - 1) / p_n;
                    s_n = (2^n) * 6 * p_n;
                end
            elseif choice == 2
                for i = 0:(n-1)
                    p_n = p_n / (1 + sqrt(1 + p_n^2));
                    s_n = (2^n) * 6 * p_n;
                end
                % Use the 2nd formula                 
            end
            if n == 0
                p_n = 1/sqrt(3);
                s_n = (2^n) * 6 * p_n;
            end
        end

        function s = p3(a)
            % This function computes the Kahan summation algorithm.
            % :param a: a vector of numbers
            % :return: summation of the vector a using Kahan summation algorithm
            
            s = inf; % Write your code here
        end

        function p4(a)
            % This function test the performance of Kahan summation algorithm against native sum.
            % :param a: a vector of numbers in double precision.
            % :return: no returns

            % Test this function with a = rand(n, 1) with various size n multiple times. 
            % Summarize your findings below.
            %
            % Findings:
            % 
            %
            %
            %
            
            single_a = single(a); % Convert a to single precision
            s = hw01.p3(a); % Kahan sum of a under double precision (regarded as truth).

            single_Kahan_s = hw01.p3(single_a); % Kahan sum of single_a under single precision.
            single_naive_s = sum(single_a); % Naive sum of single_a under single precision.
   
            disp(['Error of naive sum under single precision: ', num2str(single_naive_s-s)]);
            disp(['Error of Kahan sum under single precision: ', num2str(single_Kahan_s-s)]);
        end

        function s = p5(a)
            % For 6630. 
            % This function computes summation of a vector using pairwise summation.
            % :param a: a vector of numbers
            % :return: the summation of the vector a using pairwise summation algorithm.

            % ! You may need to create a helper function if your code uses recursion.

            % Rewrite the p4 test function to test this summation method.
            % Summarize your findings below.
            %
            % Findings: 
            %
            %
            %
            %
            %


            s = inf; % Write your code here. 

        end
    end
end