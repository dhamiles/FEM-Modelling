function elemat = ReactionElemMatrix(lambda,eID,msh)
%REACTIONELEMMATRIX Returns 2x2 local element matrix 
%   Function returns the 2x2 local element matrix for the for a given mesh
%   (msh), for the reaction coefficient (lambda) at a local element (eID) 

% Pull the value of J out of the mesh data structure at element eID
J = msh.elem(eID).J;

% Calculate the 2x2 local element matrix
elemat = [ (2*lambda*J)/3, (lambda*J)/3; (lambda*J)/3, (2*lambda*J)/3 ];

end

