function gvector = GlobalVector(f,mesh)
%GLOBALVECTOR Summary of this function goes here
%   Detailed explanation goes here

% Declare a vector of zeros for the global vector
gvector = zeros(mesh.ngn,1);

% Loop through all the elements in the mesh
for i=1:mesh.ne
    
    % Calculate the local vector at element i
    lvector = SourceElemVect(f,i,mesh);
    
    % Insert the vector into the global vector
    gvector(i:i+1,1) = gvector(i:i+1,1) + lvector; 
    
end

end

