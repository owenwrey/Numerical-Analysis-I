% Author: Owen Renyolds / owr0001@auburn.edu
% Date: 2024-09-01
% Assignment Name: hw06


classdef hw06
    methods (Static)
        
        % Problem 1

        function ret = p1(func, a, b, n, option)
            % Implement composite quadrature rules for numerical integration
            % of a function over the interval [a, b] with n subintervals.
            % The option parameter determines the type of quadrature rule to use: 
            % 1 for the midpoint rule, 2 for the trapezoidal rule, and 3 for Simpson's rule.

            %:param func: The function to integrate, provided as a function handle.
            %:param a: The lower bound of the integration interval.
            %:param b: The upper bound of the integration interval.
            %:param n: The number of subintervals to use for the integration.
            %:param option: The type of quadrature rule to use (1, 2, or 3).
            %:return: The approximate integral of the function over the interval [a, b].

            % your code here.
            if option == 1
                % your code of composite midpoint rule
                ret = 0;
                h = (b - a) / n;

                for i = 1:n
                    x_mid = a + (i - 0.5) * h;

                    ret = ret + func(x_mid) * h;
                end
            elseif option == 2
                % your code of composite trapezoidal rule
                h = (b - a) / n;
                ret = ( func(a) + func(b) ) / 2;

                for i = 1:(n-1)
                    x = a + i*h;
                    ret = ret + func(x);
                end

                ret = ret * h;

            elseif option == 3
                % your code of composite Simpson's rule
                if mod(n, 2) ~= 0
                    error('n must be even for Simpsons rule');
                end

                ret = func(a) + func(b);
                h = (b - a) / n;
                
                for i = 1:2:(n-1)
                    x = a + i * h;
                    ret = ret + 4 * func(x);
                end
                
                for i = 2:2:(n-2)
                    x = a + i * h;
                    ret = ret + 2 * func(x);
                end

                ret = ret * h / 3;


            else
                error('Invalid option: %d', option);
            end
        end


        % Problem 2

        function p2()
            % run with the following command: hw06.p2(). Do not edit this function.
            %
            % It checks the convergence of the composite quadrature rules implemented in p1.
            %
            % Here we use some examples, 
            % f_1(x) = exp(x) 
            % f_2(x) = (1 - x^2)^3, this function's first 2 derivatives at the endpoints are zero.
            % f_3(x) = (1 - x^2)^5, this function's first 4 derivatives at the endpoints are zero.
            % f_4(x) = (1 - x^2)^7, this function's first 6 derivatives at the endpoints are zero.

            % Run this function will plot the figures for the convergence of the composite quadrature rules.
            % Make comments about the results obtained from the plots. 
            %
            % > For instance, you can comment on the convergence rates of the quadrature rules, and how they compare to the theoretical rates.
            % > Here are a few example questions you can answer in your comments:
            % > Does the Runge phenomenon of f1 (Runge's function) lower the convergence rate?
            % > Does Simpson's rule have a faster convergence rate than the other two for these examples?
            % > Based on your observations, what kind of functions can have a faster convergence rate with the given composite quadrature rules?

            % Write your comments here.
            %
            % Simpson's rule converged slower than both midpoint and trapezoid for the three polynomials,
            % however Simpson's outperformed them for the exponential function.
            % None of the methods reached 8th order convergence, 6th order convergence 
            % was only reached for (1 - x^2)^7. Also, the higher the power of the polynomial,
            % the higher the convergence rate of the quadrature rule
            %

            f = {  @(x)exp(x),  @(x) (1 - x.^2 ).^3, @(x)(1 - x.^2).^5,  @(x) (1 - x.^2).^7} ;  % Define the integrand
            exact = [exp(1) - exp(-1), 32/35, 512/693 , 4096/6435];  % Define the exact integral
            n = 2.^(1:8);  % Define the number of subintervals
            for k = 1 : length(f)

                error = zeros(3, length(n));  % Initialize the error matrix with zeros

                % Calculate the approximate integral and the error for each quadrature rule and number of subintervals
                for i = 1 : length(n)
                    error(1, i) = abs(hw06.p1(f{k},-1, 1, n(i), 1) - exact(k));
                    error(2, i) = abs(hw06.p1(f{k},-1, 1, n(i), 2) - exact(k));
                    error(3, i) = abs(hw06.p1(f{k},-1, 1, n(i), 3) - exact(k));
                end

                % Plot the error against the number of subintervals using a log-log scale
                figure(k);
    
                loglog(n, error(1, :), 'r-+', 'LineWidth', 2);
                hold on;
                loglog(n, error(2, :), 'g-d', 'LineWidth', 2);
                loglog(n, error(3, :), 'b-x', 'LineWidth', 2);

                loglog(n, 1./ n.^2, 'm--', 'LineWidth', 1);
                loglog(n, 1./ n.^4, 'k-.', 'LineWidth', 1);
                loglog(n, 1./ n.^6, 'm--d', 'LineWidth', 1);
                loglog(n, 1./ n.^8, 'k--o', 'LineWidth', 1);

                xlabel('Number of subintervals');
                ylabel('Absolute error');
                title(sprintf('Convergence of composite quadrature rules for %s', functions(f{k}).function));
                legend('Midpoint rule', 'Trapezoidal rule', 'Simpson''s rule', '2nd order convergence', '4th order convergence', '6th order convergence', '8th order convergence', 'Location', 'best');
                grid on;
                hold off;
            end

        end

        
        % Problem 3

        function ret = p3(func, a, b, N, option)
            % Use your implemented Richardson extrapolation function in HW05 to implement the Romberg integration method.
            %
            % :param func: The function to integrate, provided as a function handle.
            % :param a: The lower bound of the integration interval.
            % :param b: The upper bound of the integration interval.
            % :param N: it means 2^N is the maximum number of subintervals to use for the integration. 
            %           The Romberg method will start with 2^1=2 subintervals and double the number of subintervals until 2^N
            % :param option: The type of quadrature rule to use (1, 2, or 3). See p1.
            % :return: The approximate integral of the function over the interval [a, b].

            % Note, the "powers" used in Richardson extrapolation (see hw05.m) should be [2, 4, 6, ...] for option 1 and 2. 
            % For option 3, the "powers" should be [4, 6, 8, ...].

            % your code here.
            N = 2 * N; % lol
            if option == 1 || option == 2
                powers = 2:2:2*N;
            elseif option == 3
                powers = 4:2:2*N;
            else
                error('Option must be 1, 2, or 3.')
            end
            
            R = zeros(N, N);
            %data = zeros(1, length(powers));

            for i = 1:N
                n = 2^(i-1);%powers(i);
                if option == 3 && mod(n, 2) ~= 0
                    n = n + 1; % Ensure n is even for Simpson's rule
                end
                R(i, 1) = hw06.p1(func, a, b, n, option);
            end

            
            for j = 2:N  % col
                for i = j:N  % row
                    R(i, j) = ( R(i, j-1) * 2^powers(j-1) - R(i-1, j-1) ) / (2^powers(j-1) - 1);
                end
            end
            
            ret = R(N, N);

        end



        % Problem 4
        
        function ret = p4()
            % Construct the Gauss quadrature rule using the roots of the Legendre polynomial of degree 6.
            % 
            % To evaluate Legendre polynomial of degree 6, use the helper function hw06.legendre_poly_6 defined below.
            % Its handle is @hw06.legendre_poly_6.
            % 
            % :return: A 6x2 matrix containing the roots and weights of the Gauss quadrature rule. 
            %          The first column contains the roots and the second column contains the corresponding weights.


            brackets = [-1, -3/4; -3/4, -1/4; -1/4, 0; 0, 1/4; 0.5, 3/4; 0.9, 1];
            roots = zeros(6, 1);
        
            for i = 1:6
                roots(i) = hw02.p1(@hw06.legendre_poly_6, brackets(i, 1), brackets(i, 2), 1e-14, 'regula_falsi'); 
            end
            
            weights = zeros(6, 1);
            for i = 1:6
                derivative_at_root = hw06.deriv_lengendre_poly(6, roots(i));
                weights(i) = 2 / ((1 - roots(i)^2) * derivative_at_root^2);
            end
        
            ret = [roots, weights]; 
        
        end

        % p5 is 6630 only
        function ret = p5(n)
            % For 6630 ONLY. 

            % Construct the Gauss quadrature rule using the roots of the Legendre polynomial of degree n
            %
            % :param n: The degree of the Legendre polynomial for the nodes of the Gauss quadrature rule.
            % :return: An nx2 matrix containing the roots and weights of the Gauss quadrature rule.
            % 
            % To evaluate Legendre polynomial or its derivative of a specific degree n, the handles are:
            % @(x) hw06.legendre_poly(n, x) and @(x) hw06.deriv_lengendre_poly(n, x).
            % 

            roots = zeros(n, 1);
            weights = zeros(n, 1);

            % your code here.

            ret = [roots, weights];  % Return the roots and weights of the Gauss quadrature rule

        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                                                             %
        % Helper functions below. Do not modify. You can create your own helper functions if needed.                  %
        %                                                                                                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Helper functions for p4. The following function is used to evaluate the Legendre polynomial of degree 6.
        function val = legendre_poly_6(x)
            % Compute the Legendre polynomial of degree 6 at the point x.
            %
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the Legendre polynomial of degree 6 at the point x.

            val = (231 * x^6 - 315 * x^4 + 105 * x^2 - 5) / 16;
        end

        % Helper functions for p5. The following function is used to evaluate the Legendre polynomial of degree n.
        function val = legendre_poly(n, x)
            % Compute the nth Legendre polynomial P_n at the point x.
            %
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the nth Legendre polynomial at the point x.

            if (n == 0)
                val = 1;
            elseif (n == 1)
                val = x;
            else
                val = hw06.legendre_poly(n-1, x) * x * (2 * n - 1)/n - (n - 1) * hw06.legendre_poly(n - 2, x) / n;
            end
        end

        function val = deriv_lengendre_poly(n, x)
            % Compute the derivative of the nth Legendre polynomial P_n at the point x.
            %   
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the derivative of the Legendre polynomial.
            % :return: The value of the derivative of the nth Legendre polynomial at the point x.
            val = n / (x^2 - 1) * (x * hw06.legendre_poly(n, x) - hw06.legendre_poly(n - 1, x));
        end
    end
end

