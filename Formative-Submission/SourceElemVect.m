function elevec = SourceElemVect(f,eID,msh)
%SOURCEELEMVECT Returns 2x1 local element column vector  
%   Function returns the 2x1 local element vector for a given mesh (msh), 
%   for the source term (f) at a local element (eID)
%   NOTE: this current only functions for a CONSTANT value for f!

% Pull the value of J out of the mesh data structure at element eID
J = msh.elem(eID).J;

% Calculate the 2x1 local element source term vector
elevec = [ f*J; f*J ];

end

