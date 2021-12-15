function gradient = evalBasisGrad(lnid,zeta,order)
%EVALBASISGRAD Returns the gradient of Lagrange basis function
%   evalBasisGrad returns the gradient of a first or second order Lagrange 
%   basis function at given local node id and zeta coord between [-1 1].
%   Function takes in the following inputs:
%       -lnid: the local id of a node within an element (will be [0 1] for
%              linear basis functions and [0 1 2] for quadratic functions
%       -zeta: the local position coordinate within an element to evaluate
%              the basis functions at (in range [-1 1], standard element)
%       -order: the order of the Lagrange basis functions being evaluated,
%               either 1 for linear or 2 for quadratic.

gradient=0;

if order==1
    
    sign = (-1)^(lnid+1);
    gradient = sign/2;

elseif order==2
    
    if lnid==0
        gradient = zeta - (1/2);
    elseif lnid==1
        gradient = -2*zeta;
    elseif lnid==2
        gradient = zeta + (1/2);
    end

else
   error("Invalid order provided, please input either 1 for(linear)" +...
       " or 2 for(quadratic) Lagrange basis functions");

end

end

