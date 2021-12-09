function gSourceVect = globalSourceVector(fFunc,mesh,t)
%GLOBALSOURCEVECTOR Summary of this function goes here
%   Detailed explanation goes here

% Initiate the global stiffness matrix to a matrix of zeros 
gStiffMatrix = zeros(mesh.ngn,1);

% Loop through all elements in the mesh
for e = 1:mesh.ne
    
    % Set the current space coordinate (x)
    x = mesh.nvec(e);

    % Set the values of D and lambda at the current time and space coords
    f = fFunc(x,t);
    
    % Call the local source vector function
    lSourceVect = localSourceVector(f,e,mesh);
    % Add the vector into the correct position of the global vector
    gSourceVect(e:e+1,1) = gStiffMatrix(e:e+1,1) + lSourceVect;
    
end

end

