function  [out] = selectivePlot(reqvals,xaxis,solution,step,titles,tag, ...
    xlab,unit)
%SELECTIVEPLOT plots a selection from a larger dataset
%   Function plots the solution at a number of variable values (ie it will
%   plot the result at a number of supplied time or x values). Takes the
%   following inputs:
%       -reqvals: vector of variable values to plot
%       -xaxis: the computed vector of x axis values
%       -solution: the 2D solution matrix 
%       -step: the step size used in calculating the solution
%       -tag: string, "a" for analytical or "n" for numerical

% Loop through all the required times to be plotted
for i = 1:length(reqvals)

    % Determine the index of the time's corresponding data
    index = round((reqvals(i)/step)+1);
    % Plot the solution
    name = strcat(unit,"=",num2str(reqvals(i))," [",tag,"]");
    % If numerical make line solid, if analytical make line dashed
    if tag=="n"
        plot(xaxis,solution(:,index),DisplayName=name,LineStyle='-');
    elseif tag=="a"
        plot(xaxis,solution(:,index),DisplayName=name,LineStyle='--');
    end
    hold on
end

legend('show',Location='best');
xlabel(xlab)
ylabel('Solution')
title(titles)

end