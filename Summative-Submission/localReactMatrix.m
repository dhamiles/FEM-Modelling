function lReactMatrix = localReactMatrix(lambda,e,mesh)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Declare the order of the GQ scheme
N = 2;
% Set J to a simple variable to make the code more readable
J = mesh.elem(e).J;

% Create a Gaussian Quadrature scheme to integrate basis functions
gq = CreateGQScheme(N);

% Create functions in zeta for psi0 and psi1
psi0 = @(zeta) (1-zeta)/2;
psi1 = @(zeta) (1+zeta)/2;

% Decare a functoin in x for each element
f00 = @(zeta) lambda*psi0(zeta)*psi0(zeta)*J;
f01 = @(zeta) lambda*psi0(zeta)*psi1(zeta)*J;
f10 = @(zeta) lambda*psi1(zeta)*psi0(zeta)*J;
f11 = @(zeta) lambda*psi1(zeta)*psi1(zeta)*J;

% Call the integrateGQ function to integrate the functions
Int00 = IntegrationGQ(f00,gq);
Int01 = IntegrationGQ(f01,gq);
Int10 = IntegrationGQ(f10,gq);
Int11 = IntegrationGQ(f11,gq);

% Insert the values into a 2x2 array 
lReactMatrix = [Int00,Int01;Int10,Int11];

end