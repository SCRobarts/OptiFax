function mat = struct2mat(struc)
	cells = struct2cell(struc);
	mat   = cell2mat(cells);
end