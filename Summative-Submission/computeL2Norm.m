function E = computeL2Norm(mesh,solution,time,dt,order)
%COMPUTEL2NORM Summary of this function goes here
%   Detailed explanation goes here

t = round(time*(1/dt))/(1/dt);
ti = round((t/dt)+1);

if order==1
    N=2;
elseif order==2
    N=3;
end

gq = CreateGQScheme(N);

E = 0;
for e = 1:mesh.ne
    
    i = 1+(e-1)*order; % Set the coord of first nodes solution in matrix
    J = mesh.elem(e).J;
    
    % Create a vector/field of solution values for the nodes in element
    cf=zeros(order+1,1);
    cf(1:order+1,1) = solution(i:i+order,ti);
    xf = mesh.elem(e).x;
    
    % Create a function in zeta to numerically integrate
    f = @(zeta) J*(TransientAnalyticSoln(evalField(mesh,xf,e,zeta,order)...
                ,t)-evalField(mesh,cf,e,zeta,order))^2;
    
    % Sum the total error on each iteration
    E = E + IntegrationGQ(f,gq);
    
end

