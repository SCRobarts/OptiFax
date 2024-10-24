function [x_binned,bin_ids] = xbin(xs,bins)
%
%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com

	n_bins = length(bins);
	flipped = 0;
	if length(bins)>1
		if bins(1) > bins(2)
			bins = fliplr(bins);
			flipped = 1;
		end
		bin_ids = discretize(xs,bins);
	else
		bin_ids = 1;
	end
	% bin_ids(isnan(bin_ids)) = n_bins;
	% bins(n_bins) = 0;
	bin_ids(isnan(bin_ids)) = 1;
	bins(1) = 0;
	bin_ids = uint16(bin_ids);
	x_binned = bins(bin_ids);
	% x_binned = sparse(x_binned);
	x_binned = single(x_binned);
	if flipped
		bin_ids = fliplr(bin_ids);
	end
end