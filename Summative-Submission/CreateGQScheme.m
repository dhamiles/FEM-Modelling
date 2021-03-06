function [gq] = CreateGQScheme(N)
%CREATEGQSCHEME Creates a GQ scheme of order N
%   This function returns a data structure populated with the Guass Weights
%   (.gsw) at each of the Gauss 'x' coordinates (.xpts). The data structure
%   also contains the number of points, initialised to the input value N.

gq.npts = N; 

% Currently this function works up to a maximum N value of 4
if (N > 0) && (N < 5)
     % Declare the size of the Gauss weight and point arrays
     gq.gsw = zeros(N,1);
     gq.xipts = zeros(N,1);
     
     % Populate the arrays with a case structure 
     switch N
         case 1
             gq.gsw(1) = 2;
             gq.xipts(1) = 0;
         case 2 
             gq.gsw(1) = 1;
             gq.xipts(1) = -sqrt(1/3);
             gq.gsw(2) = 1;
             gq.xipts(2) = sqrt(1/3);
         case 3
             gq.gsw(1) = 5/9;
             gq.xipts(1) = -sqrt(3/5);
             gq.gsw(2) = 8/9;
             gq.xipts(2) = 0;
             gq.gsw(3) = 5/9;
             gq.xipts(3) = sqrt(3/5);
         case 4 
             gq.gsw(1) = (18-sqrt(30))/36;
             gq.xipts(1) = -sqrt((3/7)+(2/7)*sqrt(6/5));
             gq.gsw(2) = (18+sqrt(30))/36;
             gq.xipts(2) = -sqrt((3/7)-(2/7)*sqrt(6/5));
             gq.gsw(3) = (18+sqrt(30))/36;
             gq.xipts(3) = sqrt((3/7)-(2/7)*sqrt(6/5));
             gq.gsw(4) = (18-sqrt(30))/36;
             gq.xipts(4) = sqrt((3/7)+(2/7)*sqrt(6/5));
     end
     
else 
    fprintf('Invalid number of Gauss Points (N) specified');
end

