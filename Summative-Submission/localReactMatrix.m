function lReactMatrix = localReactMatrix(lambda,e,mesh,order)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% Declare the order of the GQ scheme
if order==1
    N = 2;
elseif order==2
    N = 4;
end

% Set J to a simple variable to make the code more readable
J = mesh.elem(e).J;
% Create a Gaussian Quadrature scheme to integrate basis functions
gq = CreateGQScheme(N);

% Create placeholder matrix of zeros
lReactMatrix = zeros(order+1);

% Loop through rows of matrix
for n = 0:order
    % Loop through columns of matrix
    for m = 0:order
        % Create a 'function' in zeta to integrate with GQ
        f = @(zeta) evalField(mesh,lambda,e,zeta,order)*...
                    evalBasis(n,zeta,order)*...
                    evalBasis(m,zeta,order)*...
                    J;
        lReactMatrix(n+1,m+1) = IntegrationGQ(f,gq);
    end
end

end