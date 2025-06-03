function Mbinned = bin3d(M,ids)
	%
	%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
	
sflag = 0;

sz = size(M);
if sz(2) == 1
	M = repmat(M,1,sz(1));
	sz = size(M);
end
nz = size(M,3);

if issparse(M)
	sflag = 1;
	M = full(M);
end

if gpuDeviceCount
	M = gpuArray(M);
end

indices = arrayfun(@(l) uint16(1:l), sz, 'UniformOutput', false); %vector of indices in each dimension
[indices{:}] = ndgrid(indices{:});
indices{2} = repmat(ids,[1 1 nz]);
indices = cell2mat(cellfun(@(v) v(:), indices, 'UniformOutput', false));

Mbinned = accumarray(indices,M(:),sz,@max);
Mbinned = gather(Mbinned);

if sflag
	Mbinned = ndSparse(Mbinned);
end


end