function pathstr = OptiFaxRoot()
addonPath = matlab.internal.addons.util.retrieveAddOnsInstallationFolder;
pathstr = addonPath + filesep + 'Toolboxes' + filesep + 'OptiFax';
end