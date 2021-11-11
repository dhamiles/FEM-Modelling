function gmatrix = GlobalMatrix(D,lambda,mesh)
%GLOBALMATRIX Creates a Global Matrix 
%   Function that returns a global matrix by calculating local vectors for
%   the diffusion and reaction terms, then inserting them accordingly. It
%   also takes in the mesh that the matrix is being calculated for

% Declare a matrix of zeros for the global matrix
gmatrix = zeros(mesh.ngn);

% Loop through all the elements in the mesh
for i=1:mesh.ne
    
    % Calculate the local diffusion matrix
    lmatrix = LaplaceElemMatrix(D,i,mesh);
    % Insert the local matrix into the global matrix
    gmatrix(i:i+1,i:i+1) = gmatrix(i:i+1,i:i+1) + lmatrix;
    
    % Calculate the local diffusion matrix
    lmatrix = ReactionElemMatrix(lambda,i,mesh);
    % Insert the local matrix into the global matrix
    gmatrix(i:i+1,i:i+1) = gmatrix(i:i+1,i:i+1) - lmatrix;
    
end

end

