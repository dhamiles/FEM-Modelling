function lMassMatrix = localDiffusionMatrix(D,e,mesh)
%LOCALDIFFUSIONMATRIX Summary of this function goes here
%   
%

% Declare the order of the GQ scheme
N = 2;
% Set J to a simple variable to make the code more readable
J = mesh.elem(e).J;

% Create a Gaussian Quadrature scheme to integrate basis functions
gq = CreateGQScheme(N);

% Decare a functoin in x for each element
f00 = @(zeta) D/(4*J);
f01 = @(zeta) -D/(4*J);
f10 = @(zeta) -D/(4*J);
f11 = @(zeta) D/(4*J);

% Call the integrateGQ function to integrate the functions
Int00 = IntegrationGQ(f00,gq);
Int01 = IntegrationGQ(f01,gq);
Int10 = IntegrationGQ(f10,gq);
Int11 = IntegrationGQ(f11,gq);

% Insert the values into a 2x2 array 
lMassMatrix = [Int00,Int01;Int10,Int11];

end

