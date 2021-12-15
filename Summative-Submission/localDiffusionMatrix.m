function lDiffMatrix = localDiffusionMatrix(D,e,mesh,order)
%LOCALDIFFUSIONMATRIX Summary of this function goes here
%   
%

% Declare the order of the GQ scheme
if order==1
    N = 1;
elseif order==2
    N = 3;
end

% Set J to a simple variable to make the code more readable
J = mesh.elem(e).J;
% Create a Gaussian Quadrature scheme to integrate basis functions
gq = CreateGQScheme(N);

% Create placeholder matrix of zeros
lDiffMatrix = zeros(order+1);

% Loop through rows of matrix
for n = 0:order
    % Loop through columns of matrix
    for m = 0:order
        % Create a 'function' in zeta to integrate with GQ
        f = @(zeta) evalField(mesh,D,e,zeta,order)*...
                    evalBasisGrad(n,zeta,order)*...
                    evalBasisGrad(m,zeta,order)*...
                    (1/J);
        lDiffMatrix(n+1,m+1) = IntegrationGQ(f,gq);
    end
end

end

