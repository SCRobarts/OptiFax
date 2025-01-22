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
		YData

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
		TimeLimits
		YLimits
	end

	methods(Static)
		function tl = createTiles(panH)
			tl = tiledlayout(panH,"horizontal","TileSpacing","loose","Padding","compact");
		end
	end

	methods
		function obj = SimPlotter(optSim,ydat,ylab,kplotlims)
			arguments
				optSim	OpticalSim
				ydat 
				ylab = "Distance / (mm)";
				kplotlims = [350 2000];
			end
		%SIMPLOTTER Construct an instance of this class
			set(groot(), 'DefaultFigureUnits','pixels');
			obj.Parent = optSim;
			obj.YData = ydat;
			obj.YLabel = ylab;
			obj.SpecLimits = kplotlims;
			obj.ScreenPositions = groot().MonitorPositions + 50*[1 1 -2 -3];
			if length(obj.ScreenPositions(:,1)) < obj.Screen
				obj.Screen = length(obj.ScreenPositions(:,1));
			end
			obj.evofig;
		end
		
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

		function ax = createInOutAxes(obj,tH,type)
			ax = nexttile(tH);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			if strcmp(type,"t")
				xlabel(ax,obj.TimeLabel)
				xlim(ax, obj.TimeLimits);
			else
				xlabel(ax,obj.SpecLabel)
				xlim(ax, obj.SpecLimits);
			end
		end

		function ax = createEvoAxes(obj)
			t = obj.EvoTiles;
			ax = nexttile(t);
         	ylabel(obj.YLabel);
			disableDefaultInteractivity(ax);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			hold(ax,"on");
		end

		function ph = cplot(obj,ax,x)
			z = ones(length(obj.YData),length(x));
			ph = surf(ax,x,obj.YData,z);
			ylim(ax,obj.YLimits)
			shading(ax,'interp');
			colormap(ax,obj.CMap);
			drawnow('limitrate');
		end

		function updateplots(obj)
			optSim = obj.Parent;
			trip = optSim.TripNumber;
			tripstr = ['Round Trip Number: ', int2str(trip)];
			obj.ProgressFigure.Name = tripstr;
			obj.EvoTiles.Title.String = tripstr;
			obj.SpectralEvoPlot.ZData = optSim.IkEvoData;
			obj.TemporalEvoPlot.ZData = optSim.ItEvoData;

			obj.ioplots(1,optSim.XOutPulse);
		
			obj.ioplots(2,optSim.XInPulse);

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

		function updatetitles(obj,n,pulse)
			specstr = pulse.Name + ' Spectral in ' + pulse.Medium.Bulk.Material;
			tempstr = pulse.Name + ' Temporal in ' + pulse.Medium.Bulk.Material;
			obj.SpectralAxes(n).Title.String = specstr;
			obj.TemporalAxes(n).Title.String = tempstr;
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
			% obj.SpectralEvoPlot.ZData = optSim.StoredPulses.ESD_pJ_THz;
			% obj.TemporalEvoPlot.ZData = optSim.StoredPulses.TemporalIntensity(:,1:optSim.SimWin.Granularity:end);
			obj.SpectralEvoPlot.ZData = optSim.StoredPulses.CombinedESD_pJ_THz;
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

		function tls = get.TimeLimits(obj)
			tls = 0.95*[min(obj.Parent.SimWin.Timesfs) max(obj.Parent.SimWin.Timesfs)];
		end

		function yls = get.YLimits(obj)
			yls = [min(obj.YData) max(obj.YData)];
		end

	end

end