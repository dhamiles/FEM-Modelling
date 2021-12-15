function gStiffMatrix = globalStiffnessMatrix(Dfunc,lambdafunc,mesh,order)
%GLOBALSTIFFNESSMATRIX Summary of this function goes here
%   Detailed explanation goes here

% Initiate the global stiffness matrix to a matrix of zeros 
gStiffMatrix = zeros(mesh.ngn);

% Loop through all elements in the mesh
for i = 1:mesh.ne

    e = 1+(i-1)*order; % Set the coord of the top left of matrix
    
    % Set the values of D and lambda at the current time and space coords
    D = @(x) Dfunc(x);
    lambda = @(x) lambdafunc(x);

    % Call the local diffusion matrix function
    lDiffMatrix = localDiffusionMatrix(D,i,mesh,order);
    % Add the matrix into the correct position of the global matrix
    gStiffMatrix(e:e+order,e:e+order) = gStiffMatrix(e:e+order,e:e+order)...
                                        + lDiffMatrix;

    % Call the local reaction matrix function
    lReacMatrix = localReactMatrix(lambda,i,mesh,order);
    % Add the matrix into the correct position of the global matrix
    gStiffMatrix(e:e+order,e:e+order) = gStiffMatrix(e:e+order,e:e+order)...
                                        - lReacMatrix;
    
end
end

