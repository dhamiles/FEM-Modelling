function [gmatrix,gvector] = BoundaryConditions(bc0val,bc0type,bc1val,bc1type,matrix,mesh,vector)
%BOUNDARYCONDITIONS Function that applies boundary conditions to a global
%matrix
%   This function takes in the following inputs:
%   - bc0val: the numerical value of the boundary condition at x=0
%   - bc0type: the type (Neumann/Dirichlet) of the BC at x=0
%   - bc1val: the numerical value of the boundary condition at x=1
%   - bc1type: the type (Neumann/Dirichlet) of the BC at x=1
%   - matrix: the current global matrix
%   - mesh: the 1D linear mesh for the problem
%   - vector: the current global vector
%
%   Currently, the types can either be 'Neumann' or 'Dirichlet'

% If the x=0 bc is a Neumann
if bc0type == "Neumann"
    
    % Subtract the value into the first position of the global vector
    vector(1,1) = vector(1,1) - bc0val;

% If the x=0 bc is a Dirichlet
elseif bc0type == "Dirichlet"
    
    % Change the first position in the global vector
    vector(1,1) = bc0val;
    % Change the first row of the global matrix
    matrix(1,1) = 1;
    matrix(1,2:mesh.ngn) = zeros(1,mesh.ngn-1);
    
end

% If the x=0 bc is a Neumann
if bc1type == "Neumann"
    
    % Add the value into the first position of the global vector
    vector(mesh.ngn,1) = vector(mesh.ngn,1) + bc1val;

% If the x=0 bc is a Dirichlet
elseif bc1type == "Dirichlet"
    
    % Change the last position in the global vector
    vector(mesh.ngn,1) = bc1val;
    % Change the last row of the global matrix
    matrix(mesh.ngn,mesh.ngn) = 1;
    matrix(mesh.ngn,1:mesh.ngn-1) = zeros(1,mesh.ngn-1);
    
end
 
% Return the new global vector and matrix with applied BC's
gmatrix = matrix;
gvector = vector;

end

