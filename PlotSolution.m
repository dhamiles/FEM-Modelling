function output = PlotSolution(xmin,xmax,Ne,solution)
%PLOTSOLUTION Function that plots a given solution
%   Function plots the solution against x position for a given solution
%   vector (solution), with provided mesh details (xmin, xmax and Ne)

% Set the x spacing and calculate vector of x values
dx = (xmax - xmin) / Ne;
xvals = xmin:dx:xmax;

% Plot the solution against x
plot(xvals,solution,'-ok');

% Add simple axis labels
ylabel('Solution')
xlabel('x Position')

end

