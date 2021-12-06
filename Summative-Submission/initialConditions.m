function ICvector = initialConditions(ICs,mesh)
%INITIALCONDITIONS Summary of this function goes here
%   Detailed explanation goes here

% Initiate the vector to a vector of zeros
ICvector = zeros(mesh.ngn,1);

% Iterate through all the nodes in the mesh
for n = 1:mesh.ngn
    
    % Set the value of the IC vector at node n as per the IC equation
    ICvector(n) = ICs(mesh.nvec(n));

end

end

