function value = evalField(mesh,field,eID,zeta,order)
%EVALFIELD Evaluates the value of a field (space/material) in an element
%   evalField returns the value of a space/material field within an element
%   being modelled by either linear or quadratic Lagrangian elements. The
%   function takes the following inputs:
%       -mesh: an instance of the mesh data structure 
%       -field: the function in x (ie D = @(x)) of the field to be
%               valuated. Can be a constant or spacially varying
%       -eID: the element in the mesh at which to evaluate the field
%       -zeta: the local positional coordinate in the element
%       -order: the order of the Lagrange functions being used

% Set the length of the vectors to be filled
vecl = order + 1;
% Declare vectors of zeros for the psi and field coefficient values
psivec = zeros(1,vecl);
fcoeff = zeros(vecl,1);

value = 0;

for i = 1:vecl
    
    % Set the value of psi at local node i, at position zeta
    psivec(1,i) = evalBasis(i-1,zeta,order);
    % Set the value of the field 
    if strcmpi(class(field),'function_handle')
        % If field is function handle, call the function
        fcoeff(1,i) = field(mesh.elem(eID).x(i));
    else
        % Else, will be a vector so call the vector value at the index
        fcoeff(1,i) = field(i);
    end
    % Multiply the values and add to the value output
    value = value + psivec(1,i)*fcoeff(1,i);

end

end

