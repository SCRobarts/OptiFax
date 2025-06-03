function savesimfig(hfig,pathstr,folderstr,filestr,extstr)
arguments
	hfig
	pathstr
	folderstr
	filestr = 'simplot';
	extstr = '.png';
end

	r_now=char(datetime('now'), 'yyyy_MM_dd__HH_mm');
	folderpath = [pathstr,'\',folderstr];
	filename = [folderpath,'\',r_now,'_',filestr,extstr];
	
	while exist(folderpath,"dir") ~=7
		try
			mkdir(pathstr,folderstr);
		catch
			warning('Failed to create folder, potentially due to dropbox sync timings, trying again...')
			pause(1)
			delete(folderpath)
		end
	end

	wid = 'MATLAB:print:ExportExcludesUI';
	warning('off',wid);

	f = getframe(hfig);
	while exist(filename,"file") ~=2
		try
			imwrite(f.cdata,filename);
		catch
			warning('Failed to save fig on first try, potentially due to dropbox sync timings, trying again...')
			pause(1)
			delete(filename)
		end
	end
end