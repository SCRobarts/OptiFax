classdef SimPlotter < matlab.mixin.Copyable 
%SIMPLOTTER An object to create and manage plots in optical simulations
%  
%
%	Sebastian C. Robarts 2023 - sebrobarts@gmail.com

	properties
		Parent			OpticalSim
		Screen								= 0;
		FontSize							= 14;
		CMap								= jet;
		SpecLabel							= "Wavelength / (nm)";
		TimeLabel							= "Delay / (fs)";
		SpecLims							= [350 2000];
		TimeLims
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
		YData
	end

	methods(Static)
		function tl = createTiles(panH)
			tl = tiledlayout(panH,"horizontal","TileSpacing","compact","Padding","tight");
		end
	end

	methods
		function obj = SimPlotter(optSim,ydat)
		%SIMPLOTTER Construct an instance of this class
			obj.Parent = optSim;
			obj.YData = ydat;
			obj.TimeLims = [min(obj.Parent.SimWin.Timesfs) max(obj.Parent.SimWin.Timesfs)];
			obj.createfigure;
		end
		
		function createfigure(obj)
			obj.ProgressFigure = figure;
			set(obj.ProgressFigure,'Color',[1 1 1],'Position',[(200+(obj.Screen*1920)) 50 1200 900], 'Visible', 'on')
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
			xlim(obj.SpectralEvoAxes, obj.SpecLims);
			obj.SpectralEvoPlot = obj.cplot(obj.SpectralEvoAxes, obj.Parent.SimWin.Lambdanm, obj.Parent.IkEvoData);
			
			obj.TemporalEvoAxes = obj.createEvoAxes;
			xlabel(obj.TemporalEvoAxes,obj.TimeLabel)
			% xlim(obj.TemporalAxes,[-2 2].*(obj.Parent.Pulse.DurationCheck - obj.Parent.Delay))
			obj.TemporalEvoPlot = obj.cplot(obj.TemporalEvoAxes, obj.Parent.SimWin.Timesfs, obj.Parent.ItEvoData);
			
			obj.SpectralAxes = obj.createInOutAxes(obj.OutTiles,"l");
			[obj.SpectralMagPlot, obj.SpectralPhiPlot] = obj.Parent.XOutPulse.lplot(obj.SpecLims);

			obj.TemporalAxes = obj.createInOutAxes(obj.OutTiles,"t");
			[obj.TemporalMagPlot,obj.TemporalPhiPlot,obj.TemporalText] = obj.Parent.XOutPulse.tplot(obj.TimeLims);

			obj.SpectralAxes(2) = obj.createInOutAxes(obj.InTiles,"l");
			[obj.SpectralMagPlot(2), obj.SpectralPhiPlot(2)] = obj.Parent.XInPulse.lplot(obj.SpecLims);

			obj.TemporalAxes(2) = obj.createInOutAxes(obj.InTiles,"t");
			[obj.TemporalMagPlot(2),obj.TemporalPhiPlot(2),obj.TemporalText(2)] = obj.Parent.XInPulse.tplot(obj.TimeLims);

			drawnow('limitrate');
		end

		function ax = createInOutAxes(obj,tH,type)
			% t = obj.OutTiles;
			ax = nexttile(tH);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			% hold(ax,"on");
			if strcmp(type,"t")
				xlabel(ax,obj.TimeLabel)
				xlim(ax, obj.TimeLims);
			else
				xlabel(ax,obj.SpecLabel)
				xlim(ax, obj.SpecLims);
			end
		end

		function ax = createEvoAxes(obj)
			t = obj.EvoTiles;
			ax = nexttile(t);
         	ylabel('Distance (mm)')
			ylim([0 obj.Parent.System.Xtal.Length*1e3])
			disableDefaultInteractivity(ax);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			hold(ax,"on");
		end

		function ph = cplot(obj,ax,x,z)
			ph = surf(ax,x,obj.YData,z);
			shading(ax,'interp');
			colormap(ax,obj.CMap);
			drawnow('limitrate');
		end

		function updateplots(obj,optSim)
			trip = optSim.TripNumber;
			obj.EvoTiles.Title.String = ['Round Trip Number: ', int2str(trip)];
			obj.SpectralEvoPlot.ZData = optSim.IkEvoData;
			obj.TemporalEvoPlot.ZData = optSim.ItEvoData;

			specplot(1,optSim.XOutPulse);
			timeplot(1,optSim.XOutPulse);
		
			specplot(2,optSim.XInPulse);
			timeplot(2,optSim.XInPulse);

			drawnow

			function specplot(n,pulse)
				obj.SpectralMagPlot(n).YData = pulse.ESD_pJ_THz;
				obj.SpectralPhiPlot(n).YData = pulse.SpectralPhase;
				updatepeaks(obj.SpectralMagPlot(n));
			end

			function timeplot(n,pulse)
				obj.TemporalMagPlot(n).YData = gather(pulse.TemporalIntensity);
				obj.TemporalPhiPlot(n).YData = pulse.TemporalPhase;

				tstr = {['Energy = ', num2str(pulse.Energy*1e9,3), ' nJ'],...
						['FWHM = ', num2str(pulse.DurationCheck*1e15,3), ' fs']};

				obj.TemporalText(n).String = tstr;
			end

		end
	end

end