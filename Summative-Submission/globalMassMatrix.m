function gMassMatrix = globalMassMatrix(mesh,order)
%GLOBALMASSMATRIX Summary of this function goes here
%   Detailed explanation goes here

% Initialise the global mass matrix to a matrix of zeros
gMassMatrix = zeros(mesh.ngn);

for i = 1:mesh.ne

    e = 1+(i-1)*order; % Set the coord of the top left of matrix
    
    % Call the local mass matrix function to create the local matrix
    lMassMatrix = localMassMatrix(i,mesh,order);

    % Add the local matrix into the correct location
    gMassMatrix(e:e+order,e:e+order) = gMassMatrix(e:e+order,e:e+order)...
                                       + lMassMatrix;

end

end



