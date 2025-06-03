function [lam_um,lam_sfg,lam_dfg,sfgids,dfgids] = pcpm_lambda(simWin)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
arguments
	simWin
end
	lam_um = (fliplr(simWin.LambdanmPlot.*1e-3));	% [um]  Ascending
	[lam_sfg, sfgids] = sfg_lambda(lam_um',lam_um,lam_um);
	[lam_dfg, dfgids] = dfg_lambda(lam_um',lam_um,lam_um);
end