% script_medea_analysis2
%
% The old function was "1" because we skipped some lsm files and this one
% is a lot better organized. The files processed will be the same as can be
% found in an excel table "Bleaching_stats.xlsx".
%
% This function is also possible because of the re-structuring of the
% folders in which the files are stored. Bad images are sent to the "Bad"
% folder, and zstacks are found in the "zstack" folder.
%
% This processes the directories in pths, below, which should be
% photobleaching time series. 

clear
close all
yesplot = false;
load medea_filenames % after running "script_gathermedeafiles"


%
% Loop through each directory
%
soln = [];
for ii = 13:length(filenames)
	
	filename = filenames{ii};
	
	data = run_analyze_medea(filename,genotype,2,1);
	if ~isstruct(data)
		disp([data,' in:'])
		disp(filename)
		disp(' ')
	else
		soln = [soln;data];
	end
	disp(['ii = ',num2str(ii),' out of ',num2str(length(filenames))])
end



clk1 = clock;
clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2),'_',...
	num2strDU(clk1(4),2),'-',num2strDU(clk1(5),2),'-',num2strDU(round(clk1(6)),2)];
% clk = [num2str(clk1(1)),'-',num2strDU(clk1(2),2),'-',num2strDU(clk1(3),2)];

save(['Mat/',clk,'_Soln'],'soln')



