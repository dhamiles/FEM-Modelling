function lSourceVect = localSourceVector(f,e,mesh)
%LOCALSOURCEVECTOR Summary of this function goes here
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

% Decare a function in x for each element
f0 = @(zeta) psi0(zeta)*f*J;
f1 = @(zeta) psi1(zeta)*f*J;

% Integrate the functions
Int0 = IntegrationGQ(f0,gq);
Int1 = IntegrationGQ(f1,gq);

% Return the column vector
lSourceVect = [Int0; Int1];

end

