function elemat = LaplaceElemMatrix(D,eID,msh)
%LAPLACEELEMMATRIX Returns 2x2 local element matrix 
%   Function returns the 2x2 local element matrix for a given mesh(msh), 
%   for the diffusion coefficient (D) at a local element (eID)

% Pull the value of J out of the mesh data structure at element eID
J = msh.elem(eID).J;

% Calculate the 2x2 local element matrix
elemat = [ D/(2*J), -D/(2*J); -D/(2*J), D/(2*J) ];
       
end

