function gMassMatrix = globalMassMatrix(mesh)
%GLOBALMASSMATRIX Summary of this function goes here
%   Detailed explanation goes here

% Initialise the global mass matrix to a matrix of zeros
gMassMatrix = zeros(mesh.ngn);

for e = 1:mesh.ne

    % Call the local mass matrix function to create the local matrix
    lMassMatrix = localMassMatrix(e,mesh);

    % Add the local matrix into the correct location
    gMassMatrix(e:e+1,e:e+1) = gMassMatrix(e:e+1,e:e+1) + lMassMatrix;

end

end



