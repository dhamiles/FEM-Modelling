function psi = evalBasis(lnid,zeta,order)
%EVALBASIS Evaluates Lagrange basis functions for a given order 
%   evalBasis returns the value of a first or second order Lagrange basis
%   function at given local node id and zeta coord between [-1 1]. Function
%   takes in the following inputs:
%       -lnid: the local id of a node within an element (will be [0 1] for
%              linear basis functions and [0 1 2] for quadratic functions
%       -zeta: the local position coordinate within an element to evaluate
%              the basis functions at (in range [-1 1], standard element)
%       -order: the order of the Lagrange basis functions being evaluated,
%               either 1 for linear or 2 for quadratic.

if order==1
    
    sign = (-1)^(lnid+1);
    psi = (1+((sign*zeta)))/2;
    
elseif order==2
    
    if lnid==0
        psi = (zeta*(zeta-1))/2;
    elseif lnid==1
        psi = 1-(zeta^2);
    elseif lnid==2
        psi = (zeta*(1+zeta))/2;
    else
        error("Invalid local node ID, quadratic lnid values are 0,1 and 2")
    end

else
    error("Invalid order provided, please input either 1 for(linear)" +...
        " or 2 for(quadratic) Lagrange basis functions");

end

end

