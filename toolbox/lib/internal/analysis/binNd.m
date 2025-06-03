function Mbinned = binNd(M,ids)
	%
	%	Sebastian C. Robarts 2024 - sebrobarts@gmail.com
	
if iscell(M)
	cflag = 1;
	nz = size(M,2);
	sz = size(M{1});
else
	cflag = 0;	
	nz = size(M,3);
	sz = size(M);
end

indices = arrayfun(@(l) uint16(1:l), sz, 'UniformOutput', false); %Vector of indices in each dimension
[indices{:}] = ndgrid(indices{:});
indices{2} = ids;
indices = cell2mat(cellfun(@(v) v(:), indices, 'UniformOutput', false));

if gpuDeviceCount 
	if cflag
		M = cellfun(@gpuArray,M,UniformOutput=false);
	else
		M = gpuArray(M);
	end
end

Mbinned = cell(1,nz);
for zi = 1:nz
	if cflag
		Mzi = full(M{zi});
	else
		Mzi = full(M(zi));
	end
	% Mbinned{zi} = sparse(accumarray(indices,Mzi(:),sz,@max));
	Mbinned{zi} = sparse(accumarray(indices,Mzi(:),sz));	% sum function
	Mbinned{zi} = gather(Mbinned{zi});
end

end