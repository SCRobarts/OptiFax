classdef SimPlotter < matlab.mixin.Copyable 
%SIMPLOTTER An object to create and manage plots in optical simulations
%  
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

	properties
		Parent			OpticalSim
		Screen								= 2;
		FontSize							= 14;
		CMap								= jet;
		SpecLabel							= "Wavelength / (nm)";
		TimeLabel							= "Delay / (fs)";
		YLabel
		SpecLimits							= [350 2000];
		TempLimits							= [];
		ExtraLimits							= [];
		YData
		Type								= "Evolution";
	end
	properties (Transient)
		ProgressFigure			
		EvoTiles			
		SpectralEvoAxes	
		SpectralEvoPlot
		TemporalEvoAxes	
		TemporalEvoPlot
		OutTiles
		InTiles
		SpectralAxes
		SpectralMagPlot
		SpectralPhiPlot
		TemporalAxes
		TemporalMagPlot
		TemporalPhiPlot
		TemporalText
		ScreenPositions
	end
	properties (Dependent)
		TitleString
		TimeLimits
		YLimits
	end

	methods(Static)
		function tl = createTiles(panH)
			tl = tiledlayout(panH,"horizontal","TileSpacing","loose","Padding","compact");
		end
	end

	methods
		%% Construction
		function obj = SimPlotter(optSim,ydat,ylab,kplotlims,plottype)
			arguments
				optSim	OpticalSim
				ydat 
				ylab = "Distance / (mm)";
				kplotlims = [350 2000];
				plottype = "Evolution";
			end
		%SIMPLOTTER Construct an instance of this class
			set(groot(), 'DefaultFigureUnits','pixels');
			obj.Parent = optSim;
			obj.YData = ydat;
			obj.YLabel = ylab;
			obj.SpecLimits = kplotlims(1,:);
			if length(kplotlims(:,1))>1
				obj.TempLimits = kplotlims(2,:);
			end
			obj.ExtraLimits = obj.Parent.SimWin.PlotLimits;
			obj.Type = plottype;
			obj.ScreenPositions = groot().MonitorPositions + 50*[1 1 -2 -3];
			if length(obj.ScreenPositions(:,1)) < obj.Screen
				obj.Screen = length(obj.ScreenPositions(:,1));
			end
			if strcmp(obj.Type,"Evolution")
				obj.evofig;
			elseif strcmp(obj.Type,"Spectrogram")
				obj.specfig;
			elseif strcmp(obj.Type,"Scan")
				obj.scanfig;
			end
		end
		
		%% Core Figure Initialisation
		function evofig(obj)
			obj.ProgressFigure = figure;
			figPos = obj.ScreenPositions(obj.Screen,:);
			set(obj.ProgressFigure,'Color',[1 1 1],'Position',figPos, 'Visible', 'off')
			fontsize(obj.ProgressFigure,obj.FontSize,'points');

			posEvo = [0 0.2 1 0.6];
			ph1 = uipanel(obj.ProgressFigure,"Position",posEvo);
			obj.EvoTiles = obj.createTiles(ph1);
			title(obj.EvoTiles,'Round Trip Number: ');

			posOut = posEvo;
			posOut(2) = mod(posEvo(2)+posEvo(4),1);
			posOut(4) = 1 - posOut(2);
			ph2 = uipanel(obj.ProgressFigure,"Position",posOut);
			obj.OutTiles = obj.createTiles(ph2);

			posIn = posOut;
			posIn(1:2) = [0 0];
			ph3 = uipanel(obj.ProgressFigure,"Position",posIn);
			obj.InTiles = obj.createTiles(ph3);

			obj.SpectralEvoAxes = obj.createEvoAxes;
			xlabel(obj.SpectralEvoAxes,obj.SpecLabel)
			xlim(obj.SpectralEvoAxes, obj.SpecLimits);
			obj.SpectralEvoPlot = obj.cplot(obj.SpectralEvoAxes, obj.Parent.SimWin.LambdanmPlot);
			
			obj.TemporalEvoAxes = obj.createEvoAxes;
			xlabel(obj.TemporalEvoAxes,obj.TimeLabel)
			xlim(obj.TemporalEvoAxes,obj.TimeLimits)
			obj.TemporalEvoPlot = obj.cplot(obj.TemporalEvoAxes, obj.Parent.SimWin.TimesfsPlot);
			
			obj.SpectralAxes = obj.createInOutAxes(obj.OutTiles,"l");
			[obj.SpectralMagPlot, obj.SpectralPhiPlot] = obj.Parent.XOutPulse.lplot(obj.SpecLimits);

			obj.TemporalAxes = obj.createInOutAxes(obj.OutTiles,"t");
			[obj.TemporalMagPlot,obj.TemporalPhiPlot,obj.TemporalText] = obj.Parent.XOutPulse.tplot(1,obj.TimeLimits);

			obj.SpectralAxes(2) = obj.createInOutAxes(obj.InTiles,"l");
			[obj.SpectralMagPlot(2), obj.SpectralPhiPlot(2)] = obj.Parent.XInPulse.lplot(obj.SpecLimits);

			obj.TemporalAxes(2) = obj.createInOutAxes(obj.InTiles,"t");
			[obj.TemporalMagPlot(2),obj.TemporalPhiPlot(2),obj.TemporalText(2)] = obj.Parent.XInPulse.tplot(1,obj.TimeLimits);

			drawnow('limitrate'); 
		end

		function specfig(obj)
			% if isempty(obj.TempLimits)
			% 	obj.TempLimits = 
			% end
			obj.ProgressFigure = figure;
			figPos = obj.ScreenPositions(1,:);
			set(obj.ProgressFigure,'Color',[1 1 1],'Position',figPos, 'Visible', 'off')
			fontsize(obj.ProgressFigure,obj.FontSize,'points');
			
			posIn = [0 0.5 0.5 0.5];
			phIn = uipanel(obj.ProgressFigure,"Position",posIn);
			obj.InTiles = obj.createTiles(phIn);
			title(obj.InTiles,"Crystal In");
			
			obj.SpectralAxes = obj.createInOutAxes(obj.InTiles,"s");
			[obj.SpectralMagPlot] = obj.Parent.XInPulse.spectrogram(obj.SpectralAxes,obj.SpecLimits,tlims=obj.TimeLimits);
			title(obj.SpectralAxes,"Crystal In");

			posOut = posIn;
			posOut(1) = mod(posIn(1)+posIn(3),1);
			phOut = uipanel(obj.ProgressFigure,"Position",posOut);
			obj.OutTiles = obj.createTiles(phOut);
			delaystr = ['Delay = ', int2str(obj.Parent.Delay.*1e15), 'fs'];
			title(obj.OutTiles,delaystr);
		
			obj.SpectralAxes(2) = obj.createInOutAxes(obj.OutTiles,"s");
			obj.SpectralMagPlot(2) = obj.Parent.XOutPulse.spectrogram(obj.SpectralAxes(2),obj.SpecLimits,tlims=obj.TimeLimits);
			title(obj.SpectralAxes(2),"Crystal Out");

			% posEvo = posIn;
			% posEvo(1:2) = [0 0];
			posEvo = [0 0 1 0.5];
			phEvo = uipanel(obj.ProgressFigure,"Position",posEvo);
			obj.EvoTiles = obj.createTiles(phEvo);
			title(obj.EvoTiles,"Difference (Xtal GD offset)");

			obj.SpectralEvoAxes = obj.createInOutAxes(obj.EvoTiles,"s");
			obj.SpectralEvoPlot = obj.Parent.XDiffPulse.spectrogram(obj.SpectralEvoAxes,obj.SpecLimits,tlims=obj.TimeLimits);
			title(obj.SpectralEvoAxes,"Signal");

			obj.TemporalEvoAxes = obj.createInOutAxes(obj.EvoTiles,"s");
			pumpLims = obj.Parent.Source.Wavelength*1e9 + obj.Parent.Source.LineWidth.*[-2e9 2e9];
			obj.TemporalEvoAxes.YLim = pumpLims;
			obj.TemporalEvoPlot = obj.Parent.XDiffPulse.spectrogram(obj.TemporalEvoAxes,pumpLims,tlims=obj.TimeLimits);
			title(obj.TemporalEvoAxes,"Pump");

			% posDep = posEvo;
			% posDep(1) = mod(posEvo(1)+posEvo(3),1);
			% phDep = uipanel(obj.ProgressFigure,"Position",posDep);
		end

		function scanfig(obj)
			lnm = obj.Parent.SimWin.LambdanmPlot;
			% idlerLims = [dfg_lambda() dfg_lambda()]

			obj.ProgressFigure = figure;
			figPos = obj.ScreenPositions(1,:);
			set(obj.ProgressFigure,'Color',[1 1 1],'Position',figPos, 'Visible', 'off')
			fontsize(obj.ProgressFigure,obj.FontSize,'points');

			posEvo = [0 0 1 1];
			ph1 = uipanel(obj.ProgressFigure,"Position",posEvo);
			obj.EvoTiles = obj.createTiles(ph1);

			obj.SpectralAxes = obj.createEvoAxes;
			xlabel(obj.SpectralAxes,'Power / (W)')
			obj.SpectralMagPlot = obj.lplot(obj.SpectralAxes);
			if ~isempty(obj.ExtraLimits)
				hold(obj.SpectralAxes,"on")
				obj.SpectralMagPlot(2) = obj.lplot(obj.SpectralAxes);
				obj.SpectralMagPlot(2).Color = 'red';
				hold(obj.SpectralAxes,"off")
			end

			obj.SpectralEvoAxes = obj.createEvoAxes;
			xlabel(obj.SpectralEvoAxes(1),obj.SpecLabel)
			xlim(obj.SpectralEvoAxes(1), obj.Parent.PumpLimsnm);
			obj.SpectralEvoPlot = obj.cplot(obj.SpectralEvoAxes, within(lnm,obj.Parent.PumpLimsnm));
			clim(obj.SpectralEvoAxes(1),[-15 0]);

			obj.SpectralEvoAxes(2) = obj.createEvoAxes([1 2]);
			xlabel(obj.SpectralEvoAxes(2),obj.SpecLabel)
			xlim(obj.SpectralEvoAxes(2), obj.SpecLimits);	 % Allows for panning in signal range?
			obj.SpectralEvoPlot(2) = obj.cplot(obj.SpectralEvoAxes(2), within(lnm,obj.Parent.SignalLimsnm));
			clim(obj.SpectralEvoAxes(2),[-20 0]);

			if ~isempty(obj.ExtraLimits)
				obj.SpectralEvoAxes(3) = obj.createEvoAxes;
				xlabel(obj.SpectralEvoAxes(3),obj.SpecLabel)
				xlim(obj.SpectralEvoAxes(3), obj.ExtraLimits);
				obj.SpectralEvoPlot(3) = obj.cplot(obj.SpectralEvoAxes(3), within(lnm,obj.ExtraLimits));
				clim(obj.SpectralEvoAxes(3),[-20 0]);
			end

			cb = colorbar(obj.SpectralEvoAxes(1));
			cb.Layout.Tile = 'east';
			drawnow('limitrate'); 
		end

		%% Axes Creation
		function ax = createInOutAxes(obj,tH,type)
			ax = nexttile(tH);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			if strcmp(type,"t")
				xlabel(ax,obj.TimeLabel)
				xlim(ax, obj.TimeLimits);
			elseif strcmp(type,"l")
				xlabel(ax,obj.SpecLabel)
				xlim(ax, obj.SpecLimits);
			elseif strcmp(type,"s")
				xlabel(ax,obj.TimeLabel)
				xlim(ax, obj.TimeLimits);
				ylabel(ax,obj.SpecLabel)
				ylim(ax,obj.YLimits);
				hold on
			end
		end

		function ax = createEvoAxes(obj,tsize)
			arguments
				obj
				tsize = [1 1];
			end
			t = obj.EvoTiles;
			ax = nexttile(t,tsize);
         	ylabel(obj.YLabel);
			disableDefaultInteractivity(ax);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			hold(ax,"on");
		end

		%% Plotting
		function ph = lplot(obj,ax)
			x = ones(1,length(obj.YData));
			ph = plot(ax,x,obj.YData,'y-x');
			ylim(ax,obj.YLimits);
			ax.Color = 'k';
			ph.LineWidth = 1.5;
		end

		function ph = cplot(obj,ax,x)
			z = ones(length(obj.YData),length(x));
			ph = surf(ax,x,obj.YData,z);
			ylim(ax,obj.YLimits)
			shading(ax,'interp');
			colormap(ax,obj.CMap);
			drawnow('limitrate');
		end

		function ioplots(obj,n,pulse,pn)
			arguments
				obj		SimPlotter
				n
				pulse	OpticalPulse
				% pn = pulse.NumberOfPulses;
				pn = 1;
			end
			% obj.SpectralMagPlot(n).YData = pulse.ESD_pJ_THz(pn,:);
			obj.SpectralMagPlot(n).YData = pulse.CombinedESD_pJ_THz(pn,:);
			obj.SpectralPhiPlot(n).YData = pulse.SpectralPhase(pn,:);
			updatepeaks(obj.SpectralMagPlot(n));
			% obj.TemporalMagPlot(n).YData = gather(pulse.TemporalIntensity(pn,:));
			obj.TemporalMagPlot(n).YData = pulse.CombinedTemporalIntensity(pn,:);
			obj.TemporalPhiPlot(n).YData = pulse.TemporalPhase(pn,:);

			% tstr = {['Energy = ', num2str(pulse.Energy(pn,:)*1e9,3), ' nJ'],...
			tstr = {['Power = ', num2str(pulse.CombinedPower(pn,:),3), ' W'],...
					['FWHM = ', num2str(pulse.DurationCheck(pn)*1e15,3), ' fs']};

			obj.TemporalText(n).String = tstr;
			obj.updatetitles(n,pulse);
		end
		
		function roundtripplots(obj)
			optSim = obj.Parent;
			if optSim.ProgressPlotting
				obj.scalefigs;
			end
			trip = optSim.SimTripNumber;
			roundstr = ['Evolution over ', int2str(trip), ' round trips'];
			obj.ProgressFigure.Name = roundstr;
			obj.EvoTiles.Title.String = roundstr;

			obj.SpectralEvoPlot.ZData = optSim.StoredPulses.CombinedESD_pJ_THz;
			obj.SpectralEvoAxes.ColorScale = 'log';
			obj.TemporalEvoPlot.ZData = optSim.StoredPulses.CombinedTemporalIntensity(:,1:optSim.SimWin.Granularity:end);

			outPulse = optSim.OutputPulse;
			obj.ioplots(1,outPulse);

			inPulse = optSim.InputPulse;
			obj.ioplots(2,inPulse);
			
			if ~obj.ProgressFigure.Visible
				obj.ProgressFigure.Visible = "on";
			end
			drawnow

		end
		
		%% Internal Update Functions
		function updateplots(obj)
			optSim = obj.Parent;
			trip = optSim.TripNumber;
			tripstr = ['Round Trip Number: ', int2str(trip)];
			obj.ProgressFigure.Name = tripstr;

			switch obj.Type
				case "Evolution"
					obj.EvoTiles.Title.String = tripstr;
					obj.SpectralEvoPlot.ZData = optSim.IkEvoData;
					obj.TemporalEvoPlot.ZData = optSim.ItEvoData;
					obj.ioplots(1,optSim.XOutPulse);
					obj.ioplots(2,optSim.XInPulse);
				case "Spectrogram"
					obj.InTiles.Title.String = tripstr;
					obj.Parent.XInPulse.spectrogram(obj.SpectralAxes(1),obj.SpecLimits,tlims=obj.TimeLimits);
					obj.Parent.XOutPulse.spectrogram(obj.SpectralAxes(2),obj.SpecLimits,tlims=obj.TimeLimits);
					obj.Parent.XDiffPulse.spectrogram(obj.SpectralEvoAxes,obj.SpecLimits,tlims=obj.TimeLimits);
					pumpLims = obj.Parent.Source.Wavelength*1e9 + obj.Parent.Source.LineWidth.*[-2e9 2e9];
					obj.Parent.XDiffPulse.spectrogram(obj.TemporalEvoAxes,pumpLims,tlims=obj.TimeLimits);
				case "Scan"
					% esdpumpdep = optSim.OutputPulse.CombinedESD_pJ_THz;
					% esdout = optSim.OutputPulse.CombinedESD_pJ_THz;
					esdpumpdep = optSim.CombinedESDPumpDep;
					esdout = optSim.CombinedESDOut;
					pumpdep = within(esdpumpdep,obj.Parent.PumpLimsnm,obj.Parent.SimWin.LambdanmPlot);
					signal = within(esdout,obj.Parent.SignalLimsnm,obj.Parent.SimWin.LambdanmPlot);
					% obj.SpectralEvoPlot.ZData = optSim.CombinedESDOut;
					obj.SpectralEvoPlot(1).ZData = esd2irel(pumpdep,-100);
					obj.SpectralEvoPlot(2).ZData = esd2irel(signal,-100);
					obj.SpectralMagPlot(1).XData = obj.Parent.CombinedPowerOut;
					if ~isempty(obj.ExtraLimits)
						esdxout = optSim.CombinedESDXOut;
						spec = within(esdxout,obj.ExtraLimits,obj.Parent.SimWin.LambdanmPlot);
						obj.SpectralEvoPlot(3).ZData = esd2irel(spec,-100);
						pow = sum(spec.*1e-24,2).*obj.Parent.SimWin.DeltaNu.*obj.Parent.Source.RepetitionRate;
						obj.SpectralMagPlot(2).XData = pow;
					end

					obj.EvoTiles.Title.String = obj.TitleString;
					obj.EvoTiles.Title.Interpreter = 'none';
					obj.EvoTiles.Subtitle.String = obj.Parent.InfoString;
			end

			if ~obj.ProgressFigure.Visible
				obj.ProgressFigure.Visible = "on";
			end
			drawnow('limitrate','nocallbacks');
		end

		function updateYData(obj,ydat)
			obj.YData = ydat;
			obj.SpectralEvoPlot.YData = obj.YData;
			obj.TemporalEvoPlot.YData = obj.YData;
			ylim(obj.SpectralEvoAxes,obj.YLimits);
			ylim(obj.TemporalEvoAxes,obj.YLimits);
		end

		function updatetitles(obj,n,pulse)
			specstr = pulse.Name + ' Spectral in ' + pulse.Medium.Bulk.Material;
			tempstr = pulse.Name + ' Temporal in ' + pulse.Medium.Bulk.Material;
			obj.SpectralAxes(n).Title.String = specstr;
			obj.TemporalAxes(n).Title.String = tempstr;
		end

		function scalefigs(obj)
			optSim = obj.Parent;
			progPos = optSim.ProgressPlotter.ProgressFigure.Position;
			if progPos(3) > (obj.ScreenPositions(obj.Screen,3)/2 + 1)
				progPos(3) = progPos(3)./2;
				optSim.ProgressPlotter.ProgressFigure.Position =  progPos;
			end

			% if obj.ProgressFigure.Position(3) > obj.ScreenPositions(obj.Screen,3)
				progPos(1) = progPos(1) + progPos(3);
				obj.ProgressFigure.Position = progPos;
			% end
		end

		function updatelimits(obj)
			if ~isempty(obj.SpectralAxes) && ~strcmp(obj.Type,"Scan")
				[obj.SpectralAxes.XLim] = deal(obj.SpecLimits);
			end
			if~isempty(obj.SpectralEvoAxes)
				if length(obj.SpectralEvoAxes) > 1
					obj.SpectralEvoAxes(2).XLim = obj.SpecLimits;
					if~isempty(obj.ExtraLimits)
						obj.SpectralEvoAxes(3).XLim = obj.ExtraLimits;
					end
				end
			end
		end
		%% Update Property Functions
		function set.SpecLimits(obj,speclims)
			obj.SpecLimits = speclims;
			obj.updatelimits;
		end

		function set.ExtraLimits(obj,speclims)
			obj.ExtraLimits = speclims;
			obj.updatelimits;
		end

		%% Dependent Get Functions
		function tstr = get.TitleString(obj)
			cavname = [char(obj.Parent.System.Name)];
			xtalinf = obj.Parent.System.Xtal.InfoString;
			lasname = [char(obj.Parent.Source.Name)];
			powin	= [num2str(obj.Parent.Source.AveragePower','% .2f'), 'W'];
			lasinf	= [obj.Parent.Source.InfoString];

			tstr = [cavname,', ',xtalinf,', ',lasname,' ',powin,' ',lasinf];
		end

		function tls = get.TimeLimits(obj)
			if isempty(obj.TempLimits)
				tls = 0.95*[min(obj.Parent.SimWin.Timesfs) max(obj.Parent.SimWin.Timesfs)];
			else
				tls = obj.TempLimits;
			end
		end

		function yls = get.YLimits(obj)
			yls = [min(obj.YData) max(obj.YData)];
		end

	end

end