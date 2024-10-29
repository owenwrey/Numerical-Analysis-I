% Author: Owen Reynolds / owr0001
% Date: 2024-10-02
% Assignment Name: hw02

classdef hw02
    methods (Static)

        function [c, n] = p1(f, a, b, epsilon, name, f_prime)
            % p1: Implement numerical methods to find the root of a function
            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance
            % :param name: string, the name of the method
            % :param f_prime: function handle, the derivative of the function, only needed for Newton's method
            %
            % :return: 
            % c: real number, the root of the function
            % n: integer, the number of iterations
            if strcmp(name, 'bisection')
                n = 0;
                c = (a + b) / 2;
                while abs(f(c)) > epsilon
                    if f(a) * f(c) < 0
                        b = c;
                    else
                        a = c;
                    end
                    c = (a + b) / 2;
                    n = n + 1;
                end
            elseif strcmp(name, 'secant') % secant method
                % Write your code here
                c_n_minus_1 = a;
                c_n = b;
                n = 0;
                while abs(c_n - c_n_minus_1) > epsilon
                    c = c_n - f(c_n) * (c_n - c_n_minus_1) / ( f(c_n) - f(c_n_minus_1) );
                    c_n_minus_1 = c_n;
                    c_n = c;
                    n = n + 1;
                end
                
                
            elseif strcmp(name, 'newton') % Newton's method
                % Write your code here
                c = (b + a) / 2;
                n = 0;
                while abs(f(c)) > epsilon
                    c = c - f(c) / f_prime(c);
                    n = n + 1;
                end
            elseif strcmp(name, 'regula_falsi') % false position method
                % Write your code here
                x_n_minus_1 = a;
                x_n = b;
                n = 0;
                c = (b + a) / 2;
                while abs(f(c)) > epsilon
                    x_n = x_n_minus_1 - f(x_n_minus_1) * (x_n - x_n_minus_1) / (f(x_n) - f(x_n_minus_1));
                    n = n + 1;
                    c = x_n;
                end
            elseif strcmp(name, 'steffensen') % Steffensen's method
                % Write your code here
                c = (b + a) / 2;
                n = 0;
                while abs(f(c)) > epsilon
                    c = c - f(c) / ( (f(c + f(c)) - f(c)) / f(c));
                    n = n + 1;
                end
            end
        end

        function p2()
            
        %     summarize the iteration number for each method name in the table
        
        %     |name          | iter | 
        %     |--------------|------|
        %     |bisection     |18    |
        %     |secant        |10    |
        %     |newton        |4     |
        %     |regula_falsi  |9     |
        %     |steffensen    |12    |
            
        end

        function [c, n] = p3(f, a, b, epsilon) % For 6630 only
            % For 6630 only

            % Implement the Illinois method to find the root of a function

            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance

            % :return:
            % c: real number, the root of the function
            % n: integer, the number of iterations

            % Write your code here
            c = inf;
            n = inf;
        end

        function [c, n] = p4(f, a, b, epsilon) % For 6630 only
            % For 6630 only

            % Implement the Pegasus method to find the root of a function

            % :param f: function handle
            % :param a: real number, the left end of the interval
            % :param b: real number, the right end of the interval
            % :param epsilon: real number, the function tolerance

            % :return:
            % c: real number, the root of the function
            % n: integer, the number of iterations

            % Write your code here
            c = inf;
            n = inf;
        end
    end
end