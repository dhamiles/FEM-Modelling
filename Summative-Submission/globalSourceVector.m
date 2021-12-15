function gSourceVect = globalSourceVector(fFunc,mesh,order)
%GLOBALSOURCEVECTOR Summary of this function goes here
%   Detailed explanation goes here

% Initiate the global stiffness matrix to a matrix of zeros 
gSourceVect = zeros(mesh.ngn,1);

% Loop through all elements in the mesh
for i = 1:mesh.ne

    e = 1+(i-1)*order; % Set the coord of the top left of matrix
    
    % Set the values of D and lambda at the current time and space coords
    f = @(x) fFunc(x);
    
    % Call the local source vector function
    lSourceVect = localSourceVector(f,i,mesh,order);
    % Add the vector into the correct position of the global vector
    gSourceVect(e:e+order,1) = gSourceVect(e:e+order,1) + lSourceVect;
    
end

end

