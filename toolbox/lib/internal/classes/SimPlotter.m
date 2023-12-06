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
		Figure			
		EvoTiles			
		SpectralEvoAxes	
		SpectralEvoPlot
		TemporalEvoAxes	
		TemporalEvoPlot
		InOutTiles
		SpectralInOutAxes
		SpectralMagPlot
		SpectralPhiPlot
		TemporalInOutAxes
		TemporalMagPlot
		TemporalPhiPlot
		TemporalText
		YData
	end

	methods(Static)
		function tl = createTiles(pH)
			tl = tiledlayout(pH,"horizontal","TileSpacing","compact","Padding","tight");
		end
	end

	methods
		function obj = SimPlotter(optSim,ydat)
		%SIMPLOTTER Construct an instance of this class
			obj.Parent = optSim;
			obj.YData = ydat;
			obj.createfigure;
		end
		
		function createfigure(obj)
			obj.Figure = figure;
			set(obj.Figure,'Color',[1 1 1],'Position',[(200+(obj.Screen*1920)) 50 1200 900], 'Visible', 'on')
			fontsize(obj.Figure,obj.FontSize,'points');

			ph1 = uipanel(obj.Figure,"Position",[0 0.2 1 0.8]);
			obj.EvoTiles = obj.createTiles(ph1);
			title(obj.EvoTiles,'Round Trip Number: ');

			ph2 = uipanel(obj.Figure,"Position",[0 0 ph1.Position(3) ph1.Position(2)]);
			obj.InOutTiles = obj.createTiles(ph2);

			obj.SpectralEvoAxes = obj.createEvoAxes;
			xlabel(obj.SpectralEvoAxes,obj.SpecLabel)
			xlim(obj.SpectralEvoAxes, obj.SpecLims);
			obj.SpectralEvoPlot = obj.cplot(obj.SpectralEvoAxes, obj.Parent.SimWin.Lambdanm, obj.Parent.IkEvoData);
			
			obj.TemporalEvoAxes = obj.createEvoAxes;
			xlabel(obj.TemporalEvoAxes,obj.TimeLabel)
			% xlim(obj.TemporalAxes,[-2 2].*(obj.Parent.Pulse.DurationCheck - obj.Parent.Delay))
			obj.TemporalEvoPlot = obj.cplot(obj.TemporalEvoAxes, obj.Parent.SimWin.Timesfs, obj.Parent.ItEvoData);
			
			obj.SpectralInOutAxes = obj.createInOutAxes;
			xlabel(obj.SpectralInOutAxes,obj.SpecLabel)
			xlim(obj.SpectralInOutAxes, obj.SpecLims);
			[obj.SpectralMagPlot, obj.SpectralPhiPlot] = obj.Parent.PumpPulse.lplot(obj.SpecLims);
			% [obj.SpectralMagPlot, obj.SpectralPhiPlot] = obj.Parent.PumpPulse.lplot(obj.SpectralInOutAxes, obj.SpecLims);

			obj.TemporalInOutAxes = obj.createInOutAxes;
			xlabel(obj.TemporalInOutAxes,obj.TimeLabel)
			obj.TimeLims = [min(obj.Parent.SimWin.Timesfs) max(obj.Parent.SimWin.Timesfs)];
			% xlim(obj.TemporalInOutAxes, obj.TimeLims);
			[obj.TemporalMagPlot,obj.TemporalPhiPlot,obj.TemporalText] = obj.Parent.PumpPulse.tplot(obj.TimeLims);


			drawnow('limitrate');
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

		function ax = createInOutAxes(obj)
			t = obj.InOutTiles;
			ax = nexttile(t);
			ax.Interactions = [];
			ax.Toolbar.Visible = 'off';
			% hold(ax,"on");
		end

		function ph = cplot(obj,ax,x,z)
			ph = surf(ax,x,obj.YData,z);
			shading(ax,'interp');
			colormap(ax,obj.CMap);
			drawnow('limitrate');
		end

		function updateplots(obj,optSim)
			pulse = optSim.Pulse;
			obj.EvoTiles.Title.String = ['Round Trip Number: ', int2str(optSim.TripNumber)];
			obj.SpectralEvoPlot.ZData = optSim.IkEvoData;
			obj.TemporalEvoPlot.ZData = optSim.ItEvoData;
			obj.SpectralMagPlot.YData = pulse.ESD_pJ_THz;
			% obj.SpectralMagPlot = peaksplot(optSim.SimWin.LambdanmPlot,pulse.ESD_pJ_THz,50,obj.SpectralInOutAxes);
			obj.SpectralPhiPlot.YData = pulse.SpectralPhase;
			obj.TemporalMagPlot.YData = gather(pulse.TemporalIntensity);
			obj.TemporalPhiPlot.YData = pulse.TemporalPhase;

			tstr = {['Energy = ', num2str(pulse.Energy*1e9,3), ' nJ'],...
					['FWHM = ', num2str(pulse.DurationCheck*1e15,3), ' fs']};

			obj.TemporalText.String = tstr;
			% obj.TemporalMagPlot.YData = optSim.ItEvoData(end,:);
			% obj.TemporalPhiPlot.YData = optSim.Pulse.SpectralPhase;
			% obj.SpectralInOutAxes;
			% obj.Parent.Pulse.lplot(obj.SpecLims);
			% obj.TemporalInOutAxes;
			% obj.Parent.Pulse.tplot;
			drawnow
		end
	end

end