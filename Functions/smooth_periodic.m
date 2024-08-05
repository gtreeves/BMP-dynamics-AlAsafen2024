function [ysmooth,varargout] = smooth_periodic(varargin)
%Smooths periodic vectors
%
%function ysmooth = smooth_periodic(varargin)
%
% This function smooths a vector that is periodic.  That is, element(1) is
% equal to element("end+1"), element(2) equal to element("end+2"), etc.
% This is typically the case for data that come from a circle.
%
% Optional argument varargin can consist of these cases:
%	* [y] (If 1 argument)
%	* [x,y] (If 2 arguments and if arg(1) and arg(2) are the same size).
%		Note that if you do this, it does not make sense unless you are in
%		"duplicate_ends" (see below), so don't use this option for varargin
%		unless you are in that case.
%	* [y,nsmooth] (If 2 arguments and if arg(1) and arg(2) are not the same
%		size and if arg(2) is a double scalar)  
%	* [y,duplicate_ends] (If 2 arguments and if arg(1) and arg(2) are not
%		the same size and if arg(2) is a logical)  
%	* [x,y,nsmooth] (If 3 arguments and if arg(1) and arg(2) are same size
%		and if arg(3) is a double scalar)
%	* [x,y,duplicate_ends] (If 3 arguments and if arg(1) and arg(2) are
%		same size and if arg(3) is a logical). 
%	* [x,y,endpoints] (If 3 arguments and if arg(1) and arg(2) are
%		same size and if arg(3) is a 2-elt double). 
%	* [y,nsmooth,duplicate_ends] (If 3 arguments and if arg(1) and arg(2)
%		are not the same size) 
%	* [x,y,nsmooth,duplicate_ends] (If 4 arguments and if the fourth
%		argument is a logical)
%	* [x,y,nsmooth,endpoints] (If 4 arguments and if the fourth
%		argument is a 2-elt double)
%
% Which case we are in will be determined by the number of arguments and
% their sizes and classes.
%
% NOTE: using the "x" variable doesn't make sense if you don't have true
%	duplicate ends or if the natural endpoints of the x-coordinate are not
%	passed.
%
% "x": the x-axis vector.  If not passed, this is (1:length(y))' by
%	default.  This must be a 1D vector and not a 2D array.
% "y": the periodic vector to be smoothed (required).  This must be a 1D
%	vector and not a 2D array.
% "nsmooth": the number of points to be included in the sliding window.  If
%	not passed, this is 5 by default.  This cannot be greater than
%	length(y)/2.
% "duplicate_ends": Pass this as logical true if y(1) == y(end).
%	Otherwise, the program will assume that is not the case (i.e., default
%	false).

%	
%
% Output: "ysmooth" is the smoothened vector.

%% ========================================================================
% Checking input arguments
% =========================================================================
nArg = size(varargin,2); 
switch nArg
	case 1
		y = varargin{1};
% 		nsmooth = varargin{2};
	
	case 2
		if length(varargin{1}) == length(varargin{2})
			x = varargin{1};
			y = varargin{2};
			duplicate_ends = true;
			if y(1) ~= y(end)
				error('If passing only x and y, you must have duplicate ends')
			end
% 			nsmooth = varargin{2};
		else
			y = varargin{1};
			if isnumeric(varargin{2})
				nsmooth = varargin{2};
			elseif islogical(varargin{2})
				duplicate_ends = varargin{2};
			else
				error('The second of two inputs (after y) needs to be either a scalar or logical')
			end
		end
		
	case 3
		if length(varargin{1}) == length(varargin{2})
			x = varargin{1};
			y = varargin{2};
			if isnumeric(varargin{3})
				if isscalar(varargin{3})
					nsmooth = varargin{3};
					if ~exist('duplicate_ends','var')
						duplicate_ends = true;
					end
				elseif numel(varargin{3}) == 2
					endpoints = varargin{3};
				else
					error('The third input after x,y needs to be either a scalar, logical, or 2-elt vector')
				end
			elseif islogical(varargin{3})
				duplicate_ends = varargin{3};
			else
				error('The third input after x,y needs to be either a scalar, logical, or 2-elt vector')
			end
		else
			y = varargin{1};
			if isnumeric(varargin{2})
				nsmooth = varargin{2};
			else
				error('The second of three inputs (after y) needs to be a scalar')
			end
			if islogical(varargin{3})
				duplicate_ends = varargin{3};
			else
				error('The third of three inputs (after y,nsmooth) needs to be a logical')
			end
		end
		
	case 4
		x = varargin{1};
		y = varargin{2};
		nsmooth = varargin{3};
		
		if islogical(varargin{4})
			duplicate_ends = varargin{4};
		elseif isnumeric(varargin{4}) 
			if numel(varargin{4}) == 2
				endpoints = varargin{4};
			else
				error('The fourth input needs to be either a 2-elt vector or a logical')
			end
		else
			error('The fourth input needs to be either a 2-elt vector or a logical')
		end
		
	otherwise
		error('Too many inputs')
end

%
% Check "duplicate_ends" input
%
if ~exist('duplicate_ends','var')
	duplicate_ends = false;
end
if duplicate_ends && y(1) ~= y(end)
	error('If you pass "duplicate_ends" as "true", you must really have duplicate ends')
end

%
% Check y input
%
if ~isrowvec(y) && ~isrowvec(y')
	error('Your input to be smoothed must be a vector, not a 2D array')
end
if isrowvec(y)
	wasrow = true;
	y = y';
else
	wasrow = false;
end
n = length(y);

%
% Check nsmooth input
%
if ~exist('nsmooth','var')
	nsmooth = 5;
end
if nsmooth > n/2
	error('You are trying to smooth with too many points.')
end


%% ========================================================================
% Running the algorithm
% =========================================================================

if exist('x','var')
	
	%
	% Check "x" input
	%
	if ~isrowvec(x) && ~isrowvec(x')
		error('Your input vectors must be 1D, not 2D arrays')
	end
	x = x(:);
	dx = diff(x);
	
	if duplicate_ends
		x_begin = x(1);
		x_end = x(end);
		x(end) = [];
		y(end) = [];
		
% 		x1 = [x(end-nsmooth+1:end);x;x(1:nsmooth)];
	elseif exist('endpoints','var')
		x_begin = endpoints(1);
		x_end = endpoints(2);
	else
		error('If you pass "x", you must either have duplicate ends, or pass "endpoints"')
	end
	
	xleft = x_begin - flipud(cumsum(flipud(dx(end-nsmooth+1:end))));
	xright = x_end + [0;cumsum(dx(1:nsmooth-1))];
	x1 = [xleft;x;xright];
	
	
	y1 = [y(end-nsmooth+1:end);y;y(1:nsmooth)];
	% x1 = [x;x;x];
	% y1 = [y;y;y];
	ysmooth1 = smooth(x1,y1,nsmooth);

else
	if duplicate_ends
		y(end) = [];
	end
	y1 = [y(end-nsmooth+1:end);y;y(1:nsmooth)];
	ysmooth1 = smooth(y1,nsmooth);
	
end
ysmooth = ysmooth1(nsmooth+1:end-nsmooth);



%% ========================================================================
% Tidying up
% =========================================================================

if duplicate_ends
	ysmooth(end+1) = ysmooth(1);
end
if wasrow
	ysmooth = ysmooth';
end
if exist('endpoints','var') && exist('x1','var')
	[x_out,y_out] = repeat_remove(x1,ysmooth1);
	y_endpoints = interp1(x_out,y_out,[x_begin x_end]);
	varargout = {y_endpoints(1),y_endpoints(2)};
else
	varargout = {};
end














