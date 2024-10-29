% Author: Owen Reynolds / owr0001
% Date: 2024-09-01
% Assignment Name: hw04


classdef hw04
    methods (Static)
        function y = p1(data, eval)
            % Implement the divided difference method to interpolate the data points, then evaluate the polynomial at the given point.
            % :param data: a matrix of size n x 2, where the first column is the x values and the second column is the y values
            % :param eval: a column vector of x values at which to evaluate the polynomial
            % :return: a vector of y values corresponding to the evaluation points

            n_data = size(data, 1);
            n_eval = size(eval, 1);
            y = zeros(n_eval, 1);

            f = zeros(n_data); % make a matrix f the size of our data points
            f(:,1) = data(:,2); % populate the first column of f with the y value data
            x = data(:, 1); % extract x data points into a vector

            
            for s = 2:n_data
                for r = s:n_data
                    f(r, s) = (f(r, s-1) - f(r-1, s-1)) / (x(r) - x(r-s+1)); % divided difference formula
                end
            end

            div_dif = zeros(n_data, 1); % vector of 0s n_data long

            for i = 1:n_data
                div_dif(i) = f(i, i); % set the div. dif. to the diagonals of f
            end

            for h = 1:n_eval
                p_list = ones(n_data, 1); % will be a list of the Newton form terms without their coefficients
                for j = 2:n_data                    
                    p_list(j) = p_list(j-1) * (eval(h) - x(j-1)); % evaluates the Newton form w/o coeffs
                end 
                y(h) = div_dif' * p_list; % multiplies the div. dif. coeffs
                %                           by their respective Newton term
                %                           and then sums these terms
                %                           (a dot product essentially)
            end

                    





            % Write your code here
            

        end


        function y = p2(data, eval)
            % for 6630 ONLY
            % Implement the divided difference method to interpolate the data 
            % points, then evaluate the polynomial at the given point.
            %
            % :param data: a cell array of size n x 1, each cell is a vector 
            %              like (x, y0, y1, y2, y3, ..., ym). 
            %              y0 = f(x), 
            %              y1 = f'(x), 
            %              y2 = f''(x),
            %              ... ,
            %              ym = f^{(m)}(x).
            % 
            %              Note, different cells may have different lengths.
            %
            % :param eval: a vector of x values at which to evaluate the polynomial
            % :return: a vector of y values corresponding to the evaluation points

            
            n = length(data);
            y = zeros(size(eval));

            % write your code here.
            
        end
    end
end
