function [nu_central,bandwidth] = line2band(lambda_central,linewidth)
%
%	Sebastian C. Robarts 2025 - sebrobarts@gmail.com

	nu_central = c ./ lambda_central;
	lambda_edges = lambda_central + [-linewidth linewidth]./2;
	nu_edges = c ./ lambda_edges;
	bandwidth = range(nu_edges);

end