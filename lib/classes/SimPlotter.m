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
		Figure			
		Tiles			
		SpectralAxes	
		SpectralPlot
		TemporalAxes	
		TemporalPlot
		YData
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

			ph1 = uipanel(obj.Figure,"Position",[0 0.1 1 0.8]);
			obj.Tiles = tiledlayout(ph1,"horizontal","TileSpacing","compact","Padding","tight");
			
			obj.SpectralAxes = obj.createaxes;
			xlabel(obj.SpectralAxes,"Wavelength (nm)")
			xlim(obj.SpectralAxes, [350 2000]);
			obj.SpectralPlot = obj.cplot(obj.SpectralAxes, obj.Parent.SimWin.Lambdanm, obj.Parent.IkEvoData);
			
			obj.TemporalAxes = obj.createaxes;
			xlabel(obj.TemporalAxes,"Delay (s)")
			% xlim(obj.TemporalAxes,[-2 2].*(obj.Parent.Pulse.DurationCheck - obj.Parent.Delay))
			obj.TemporalPlot = obj.cplot(obj.TemporalAxes, obj.Parent.SimWin.Times, obj.Parent.ItEvoData);
			
			drawnow('limitrate');
		end
		
		function ax = createaxes(obj)
			t = obj.Tiles;
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
	end

end