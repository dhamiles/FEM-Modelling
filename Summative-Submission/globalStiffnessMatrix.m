function gStiffMatrix = globalStiffnessMatrix(Dfunc,lambdafunc,mesh,t)
%GLOBALSTIFFNESSMATRIX Summary of this function goes here
%   Detailed explanation goes here

% Initiate the global stiffness matrix to a matrix of zeros 
gStiffMatrix = zeros(mesh.ngn);

% Loop through all elements in the mesh
for e = 1:mesh.ne
    
    % Set the current space coordinate (x)
    x = mesh.nvec(e);

    % Set the values of D and lambda at the current time and space coords
    D = Dfunc(x,t);
    lambda = lambdafunc(x,t);

    % Call the local diffusion matrix function
    lDiffMatrix = localDiffusionMatrix(D,e,mesh);
    % Add the matrix into the correct position of the global matrix
    gStiffMatrix(e:e+1,e:e+1) = gStiffMatrix(e:e+1,e:e+1) + lDiffMatrix;

    % Call the local reaction matrix function
    lReacMatrix = localReactMatrix(lambda,e,mesh);
    % Add the matrix into the correct position of the global matrix
    gStiffMatrix(e:e+1,e:e+1) = gStiffMatrix(e:e+1,e:e+1) - lReacMatrix;
    
end

