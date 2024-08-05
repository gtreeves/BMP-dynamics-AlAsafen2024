function [data,F] = ftn_gradfit(data,varargin)
%Fits Dorsal/Nuc data (from lish_border) to a time/z-series of Gaussians
%
%function [data,F] = lish_gradfit(data,varargin)
%
% This function reads in the data output from "lish_border.m" and fits
% them to determine a,b,sig, and mu for each timepoint and z-slice (see
% Liberman et al, 2009).  The best fit mu for the whole time/z-series
% becomes "s_mid", which is considered the ventral midline for each
% timepoint. Then these four fields (A,B,Sig,s_mid) are appended to the end
% of "data", which is saved back to the mat file called:
%	[data.pth,'data.mat']
%
% "data": structure containing the data and metadata.
%
% Optional argument varargin can consist of these things, in this order:
%	(1) "yesplot": whether you want to plot the outcome.  Default, "false".
%		If this is not specified, but you still want to specify other
%		arguments, put empty brackets -- [] -- in place of this argument.
%
% "data": the output structure is the same as the input, plus the fit
%	parameters
% "F": the movie of the plots, only made if asked for by "yesplot"

%
% Unpacking varargin.
%
nArg = size(varargin,2); iArg = 1;
if nArg >= iArg && ~isempty(varargin{iArg})
	yesplot = varargin{iArg}; else
	yesplot = false;
end%, iArg = iArg + 1;

% tstart = data.tstart;
% zstart = data.zstart;
% tend = data.tend;
% zend = data.zend;
% T = tend - tstart + 1;
% if data.thstart ~= data.thend, error('not ready for theta stack'), end
% D = zend - zstart + 1;
T = data.D;
D = 1;
dz = 0;

S = data.S;
RatioS = data.R2;

dt = data.dt; % timestep in seconds.
z0 = data.z0; % how deep the first slice is
A = zeros(T,D); B = A; Mu = A; Sig = A; Rsquare = A;
dA = A; dB = A; dSig = A;
st = A; mo = A;
Rsmooth = cell(T,D);
x = linspace(-1,1);
t = linspace(0,(T-1)*dt,T)'/60; % time vector in minutes
z = linspace(z0,z0+(D-1)*dz,D)'; % depth in microns

%
% First we must find an estimate for the midline.
%
for i = 1:T
	for j = 1:D
% 		S{i,j} = S{i,j}/pi/data.R(i,j);
		s = S{i,j};
		r = RatioS{i,j};
		
		%
		% Smoothing by "p" points
		%
		p = 5;
		s1 = [s(end-p+1:end)-2;s;s(1:p)+2];
		r1 = [r(end-p+1:end);r;r(1:p)];
		rsmooth = smooth(s1,r1,p);
		r = rsmooth(p+1:end-p);
		
		Rsmooth{i,j} = r;
		st(i,j) = std(r);
		mo(i,j) = mode(r);
	end
end

%
% We will eventually calculate the symmetry score, but we only want to do
% that on time points that have clear peaks. This means their distributions
% have a large standard deviation compared to the mode. (Shouldn't I have 
% used median?)
%
stm = st./mo;
m = max(stm(:));
% idx = find(stm > 0.25*m);
idx = find(stm > 0.5*m);

r = []; s = [];
for i = idx'
	r = [r;Rsmooth{i}]; %#ok<AGROW>
	s = [s;S{i}]; %#ok<AGROW>
end
[s,i] = sort(s);
r = r(i);

%
% Calculating the symmetry score
%
th = linspace(-pi,pi,301)';

p = 10; % smooth (even more) by 10 points
s1 = [s(end-p+1:end)-2;s;s(1:p)+2];
r1 = [r(end-p+1:end);r;r(1:p)];
rsmooth = smooth(s1,r1,p);
r2 = rsmooth(p+1:end-p);

% f = circnormpdf(pi*s,th',0.1);
f = circnormpdf(pi*s,th',10);
symscore = f.*repmat(r2,1,length(th));
symscore = sum(symscore); [~,imax] = max(symscore,[],2);


%
% Fitting the timecourse of dl gradients to a timecourse of gaussians.  We
% first do a rough pass, and then examine the timepoints with the best-fit
% r-square value.  The average of the presumptive midlines for these will
% become the set-in-stone midline value for the whole movie, and then a
% second-pass is run using that value fixed.
%
s_mid1 = th(imax)/pi;
% s_mid1 = -0.4175;
s_mid2 = 0;
mu = 0; muL = -.3; muU = .3;
sig = 0.15; sigL = 0.05; sigU = 0.4; sigL1 = sigL; sigU1 = sigU;

for g = 1:3
	for i = 1:T
		for j = 1:D
			s = S{i,j};
			s = mod(s-s_mid1-s_mid2+1,2) - 1;
			r = Rsmooth{i,j};
			rmax = max(r); rmin = min(r);
			
			if g < 3
				a = rmax - rmin; aL = 0; aU = 10*a;
				b = rmin; bL = 0; bU = a/2 + rmin;
			else
				a = A(i,j); aL = a - 1e-4; aU = a + 1e-4;
				b = B(i,j); bL = b - 1e-4; bU = b + 1e-4;
			end
			
			f = fittype('a*exp(-(x-mu)^2/2/sigma^2)+b');
			opts = fitoptions('Method','NonlinearLeastSquares',...
				'Startpoint',[a b mu sig],...
				'Lower',[aL bL muL sigL],...
				'Upper',[aU bU muU sigU]);
			try
				[cfun,gof] = fit(s,r,f,opts);
				coeffvals = coeffvalues(cfun);
				cint68 = confint(cfun,(1-2*(1-normcdf(1,0,1))));
			catch
				coeffvals = NaN(1,4);
				gof.rsquare = NaN;
				cint68 = NaN(2,4);
			end
			
			A(i,j) = coeffvals(1);
			B(i,j) = coeffvals(2);
			Mu(i,j) = coeffvals(3);
			Sig(i,j) = coeffvals(4);
			Rsquare(i,j) = gof.rsquare;
			
			if g < 3
				dA(i,j) = coeffvals(1) - cint68(1,1);
				dB(i,j) = coeffvals(2) - cint68(1,2);
			end
			dSig(i,j) = coeffvals(4) - cint68(1,4);
			
		end
	end
	if g == 1
		idx = find(Rsquare > 0.95);
		if isempty(idx)
			[~,kk] = sort(Rsquare);
			idx = kk(end-9:end);
		end
		s_mid2 = meanDU(Mu(idx));
		sig = meanDU(Sig(idx));
		sigL = sig - 1e-4; sigU = sig + 1e-4;
		muL = -1e-4; muU = 1e-4;
	elseif g == 2
		sigL = sigL1; sigU = sigU1;
	end
end
s_mid = s_mid1 + s_mid2; % our ventral midline, in s.

%
% Plotting and making movie, if asked for.
%
if yesplot
	figure
	set(gcf,'Position',[100 100 400 300],...
		'Paperpositionmode','auto','Color',[1 1 1])
	YL = 1.2*max(A(:)+B(:));
	for i = 1:T
		for j = 1:D
			
			%
			% Extract and plot points
			%
			s = mod(S{i,j}-s_mid+1,2) - 1;
			[s,isort] = sort(s);
			r = RatioS{i,j};
			r = r(isort);
			plot(s,r,'.')
			hold on
			
			%
			% Plot Gaussian fit
			%
			y = A(i,j)*exp(-x.^2/2/Sig(i,j)^2)+B(i,j);
			if isfield(data,'M') && ~isnan(data.M)
				y = y + M(i,j)*abs(x);
			end
			plot(x,y)
			
			%
			% Plot smoothed curve
			%
			p = 5;
			rsmooth = Rsmooth{i,j};
			rsmooth = rsmooth(isort);
			plot(s,rsmooth)
			
			%
			% Add time annotation
			%
% 			timestr = sprintf('%6.2f min, %6.2f microns',t(i),z(j));
			timestr = sprintf('%6.2f min',t(i));
% 			htb = annotation('textbox',[0.5450 0.8105 0.3500 0.1000]);
			htb = annotation('textbox',[0.3050 0.8105 0.6000 0.1000]);
			set(htb,'String',timestr,'Fontsize',12,...
				'HorizontalAlignment','right','LineStyle','none')
			
			%
			% Other annotations
			%
			set(gca,'Fontsize',12)
			xlabel('s [fractional half circumference]')
			ylabel('r [relative intensity]')
			
			hold off			
			xlim([-1,1])		
			if ~exist('ymin','var')
				ymin = 0;
			end
			if ~exist('ymax','var')
				ymax = YL;
			end
			ylim([ymin ymax])
			
			%
			% Cleaning up
			%
			F(i,j) = getframe(gcf); %#ok<AGROW>
			delete(htb)
		end
	end
else
	F = 'no plot asked for';	
end

%
% Concluding remarks
%
data.t = t;
data.z = z;
data.A = A;
data.B = B;
data.Sig = Sig;
data.s_mid = s_mid;
data.gof = Rsquare;
data.dA = dA;
data.dB = dB;
data.dSig = dSig;
try
	save([data.pth,data.prefix,'data.mat'],'data');
end

