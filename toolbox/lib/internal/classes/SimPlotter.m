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
		SpecLabel							= "Wavelength (nm)";
		TimeLabel							= "Delay (s)";
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
		SpectralInOutPlot
		TemporalInOutAxes
		TemporalInOutPlot
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
			obj.TemporalEvoPlot = obj.cplot(obj.TemporalEvoAxes, obj.Parent.SimWin.Times, obj.Parent.ItEvoData);
			
			obj.SpectralInOutAxes = obj.createInOutAxes;
			xlabel(obj.SpectralInOutAxes,obj.SpecLabel)
			xlim(obj.SpectralInOutAxes, obj.SpecLims);
			% obj.SpectralInOutPlot = obj.Parent.PumpPulse.lplot(obj.SpecLims);
			obj.SpectralInOutAxes = obj.Parent.PumpPulse.lplot(obj.SpecLims);

			obj.TemporalInOutAxes = obj.createInOutAxes;
			xlabel(obj.TemporalInOutAxes,obj.TimeLabel)
			% obj.TemporalInOutPlot = obj.Parent.PumpPulse.tplot;
			obj.TemporalInOutAxes = obj.Parent.PumpPulse.tplot;

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
			obj.EvoTiles.Title.String = ['Round Trip Number: ', int2str(optSim.TripNumber)];
			obj.SpectralEvoPlot.ZData = optSim.IkEvoData;
			obj.TemporalEvoPlot.ZData = optSim.ItEvoData;
			% obj.SpectralInOutAxes;
			% obj.Parent.Pulse.lplot(obj.SpecLims);
			% obj.TemporalInOutAxes;
			obj.Parent.Pulse.tplot;
			drawnow
		end
	end

end