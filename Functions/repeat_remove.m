function [x_out,y_out,V_remove,s_out,V_out] = repeat_remove(x,y)
%
%
%function [x_out,y_out,V_out] = repeat_remove(x,y)
%
% This function serves to remove repeats in the independent variable, "x",
% and the value of the corresponding dependent variable "y" for that value
% of "x" will be the average of all "y" values that corresponded to that
% value of "x" before the repeats were removed.
%
% This function is motivated by cases in which we have lots of experimental
% data on a finite mesh.  Some of the x-values are bound to be repeated.
% This is a problem for functions like "interp1" or "spline".  So, we find
% all cases in which repeated values of "x" are found and remove the
% redundant instances.  Since this also means the y-values that correp to
% that x-value (which are in general different) will be replaced by the
% mean of those values. We also output the standard deviation of those
% y-values in the output "s_out"
%
% Inputs:
% "x": your independent variable.
% "y": your dependent variable.
%
% Outputs:
%
% "x_out": the output independent variable, which not only has the repeats
%	removed, but is also sorted ascending.
% "y_out": the output dependent variable, which has the values
%	corresponding to the repeats replaced by their average, and is also
%	sorted according to the sorting of "x_out".
% "V_remove": the indices of the original vector that were removed.
% "s_out": when the output dependent variable had the values corresponding
%	to the repeats replaced by their average, we also calculated the
%	standard deviation of the repeated values.
% "V_out": the indices of the input vector (after sorting) that are
%	repeated.

[x,isort] = sort(x);
if isrowvec(x), x = x'; wasrowvec = true; else wasrowvec = false; end
if ~exist('y','var')
	y = zeros(size(x));
else
	if (isrowvec(y) && ~wasrowvec) || (~isrowvec(y) && wasrowvec)
		error('"x" and "y" inputs must be the same size')
	elseif isrowvec(y)
		y = y';
	end
	y = y(isort,:);
end

s = NaN(size(y));
V_out = repeatcheck(x);
for j = 1:length(V_out)
	x(V_out{j}(2:end)) = NaN;
	s(V_out{j}(1),:) = std(y(V_out{j},:));
	y(V_out{j}(1),:) = mean(y(V_out{j},:));
	y(V_out{j}(2:end),:) = NaN;
	
end
V_remove = isnan(x);
x_out = x(~V_remove);
y_out = y(~V_remove,:);
s_out = s(~V_remove,:);

if wasrowvec
	x_out = x_out';
	y_out = y_out';
end