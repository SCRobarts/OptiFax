%% gnlsefibsim
function [Z, AT, AW, W, B] = gnlsefibsim(laser,simWin,fibre)

tdisp = 1e12;
lambdanm = simWin.Lambdanm;
taxis = simWin.TemporalRange;

laser.simulate(simWin);
laser.Pulse.Radius = fibre.ModeFieldDiameter./2;
laser.Pulse.plot;
drawnow
fibre.simulate(simWin);

P0 = laser.Pulse.PeakPower;
nr = laser.Pulse.Medium.Bulk.RefractiveIndex; 
Esq2I = nr.*c.*eps0 / 2;
ITmag = Esq2I .* (abs(laser.Pulse.TemporalField)).^2;
ITphi = laser.Pulse.TemporalPhase;

ITmag = ITmag .* laser.Pulse.Area;
IT = ITmag .* exp(1i*ITphi);

AT = sqrt(IT);
% AT = sqrt(laser.Pulse.TemporalIntensity .* laser.Area);

t = simWin.Times;
w0 = simWin.ReferenceOmega;

gamma0 = fibre.Gamma0;
betas = fibre.Betas;
fR = fibre.RamanFraction;
hR = fibre.RamanResponse;
L = fibre.Length;
loss = 0;
nplots = 21;

[Z, AT, AW, W, B] = gnlsefib(t, AT, w0, gamma0, betas, ...
                             loss, fR, hR, L, nplots);

laser.Pulse.TemporalRootPower = (AT(end,:));

% === plot output
figure();

WL = 2*pi*c./W; iis = (WL>450e-9 & WL<1600e-9); % wavelength grid
lIW = 10*log10(abs(AW).^2 .* 2*pi*c./WL'.^2);	% log scale spectral intensity
mlIW = max(max(lIW));							% max value, for scaling plot
subplot(1,2,1);             
pcolor(WL(iis).*1e9, Z, lIW(:,iis));			% plot as pseudocolor map
clim([mlIW-40.0, mlIW]);  xlim([450,1600]); 
colormap jet; shading interp; 
xlabel('Wavelength / nm'); ylabel('Distance / m');

lIT = 10*log10(abs(AT).^2); % log scale temporal intensity
mlIT = max(max(lIT));       % max value, for scaling plot
subplot(1,2,2);
pcolor(t.*tdisp, Z, lIT);   % plot as pseudocolor map
clim([mlIT-40.0, mlIT]);  xlim([-5,5]); 
colormap jet; shading interp;
xlabel('Delay / ps'); ylabel('Distance / m');

figure(7),subplot(311),plot(t*tdisp,abs(AT(end,:)).^2 ./ laser.Pulse.Area)
		hold on
		plot(t*tdisp,abs(AT(1,:)).^2 ./ laser.Pulse.Area)
        xlabel('Time (ps)')
        ylabel('Power (W)')
        title(['z = ' num2str(L) 'm'])
        axis([-taxis*tdisp/2 taxis*tdisp/2 0 P0*1.2/laser.Pulse.Area])
		hold off

figure(7),subplot(312),plot((W/(2*pi))*1e-12,lIW(end,:))
		hold on
		plot((W/(2*pi))*1e-12,lIW(1,:))
		xlabel('frequency (THz)'); ylabel('Spectrum (a.u.)')
		set(gca,'xlim',[230 500])
		set(gca,'ylim',[mlIW-30.0, mlIW])
		hold off

figure(7),subplot(313),plot(lambdanm,abs(AW(end,:)).^2)
		hold on
		plot(lambdanm,abs(AW(1,:)).^2)
		xlabel('\lambda (nm)'); ylabel('Spectrum (a.u.)')
		set(gca,'xlim',[350,1600])
		hold off

end