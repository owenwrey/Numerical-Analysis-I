% Author: Owen Reynolds / owr0001
% Date: 2024-09-05
% Assignment Name: hw00

% The following class defines 3 functions for each problem respectively.
% Please follow the instruction inside each function. 

classdef hw00 
    methods (Static)

        function a_m = p1(m)
            % This function takes an integer m and returns the term a_m in the sequence defined by 
            % a_0 = 0, a_1 = 1, a_2 = 1, and a_n = a_{n-1} + a_{n-2} + a_{n-3} for n >= 3.
            % :param m: an integer
            % :return: the m-th term in the sequence
            n = m + 1;
            if m < 0
                error('m must be a non-negative integer')
            elseif m == 0 
                a = [0]; % Write your code here
            elseif m == 1
                a = [0, 1];
            elseif m == 2
                a = [0, 1, 1];
            else
                
                a = zeros(1, n);
                a(1) = 0;
                a(2) = 1;
                a(3) = 1;
                    for i = 4:n
                        a(i) = sum(a);
                    end
            end
            %sumA = sum(a);
            a_m = a(n);
        end

        function det_A = p2(A)
            % This function takes a matrix A of size n x n and returns the determinant of A.
            % :param A: a matrix of size n x n
            % :return: the determinant of A

            if size(A,1) ~= size(A,2)
                error('A must be a square matrix')
            else
                det_A = det(A);% Write your code here, note when you call p2 function inside your function, you need to call it like this: hw00.p2(B), where B is a matrix.
            end
        end

        function p3()
            % This function should have a run time about 1 second.
            % :return: no returns
            
            pause(1);
            
            % Write your code here
        end
    end
end
