function out = IntegrationGQ(func,gq)
%INTEGRATIONGQ Calculates the definite integral of a function

% Set the value of N
N = gq.npts;

% Initialise a summing variable sum to 0
sum = 0;

% Loop through each weight/point in the gq scheme
for i = 1:N
    % Add the new value to the sum variable
    sum = sum + gq.gsw(i)*func(gq.xipts(i));
end

% Output the sum variable
out = sum;

end

