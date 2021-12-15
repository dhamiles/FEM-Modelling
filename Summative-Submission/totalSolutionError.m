function totalE = totalSolutionError(mesh,solution,dt,order)
%TOTALSOLUTIONERROR Summary of this function goes here
%   Detailed explanation goes here

[~,cols] = size(solution);
t = 0;
totalE = 0;

for ti = 1:cols
    
    totalE = totalE + dt*((computeL2Norm(mesh,solution,t,dt,order)));

    t = t + dt;
end

