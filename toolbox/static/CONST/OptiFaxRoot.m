function pathstr = OptiFaxRoot(devFlag)
arguments
	devFlag = 0;
end
if devFlag
	pathstr = string(pwd);
	pathstr = pathstr + filesep + 'toolbox';
else
	addonPath = matlab.internal.addons.util.retrieveAddOnsInstallationFolder;
	pathstr = addonPath + filesep + 'Toolboxes' + filesep + 'OptiFax';
end
end