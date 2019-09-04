%% Lorentzian fitting
%  This function fits data 'sca' of any scattering spectrum to a Lorentian
%  lineshape. The inital conditions are fully calculated, therefore no more
%  user input is needed for initial conditions of fitting.
%  input:       sca ---- the scattering spectrum
%               enei---- the wavelength vector (same length as sca)
%  output:      fit ---- fitted curve, the same length as enei
%               fit_Lorentz ---- fitted parameters ( Background linewidth/2pi peak linewidth  )
%               ndir ---- the subfolder or naming given by the user used
%               for saving results
%               type ----- choose to plot the sca of abs spectrum
%               PlotResults ----- if 'Yes', do plotting; if 'No', plotting
%               is held.
function [ fit, fit_Lorentz, lorentz_wl ] = LorentzFitYW( sca, enei, ndir, type )


for j = 1:length(sca)
    if sca(j) <= 0
        sca(j) = 0;
        enei(j) = 0;
    end
end

lorentz_eV_plot_only =  linspace( 1248/min(enei), 1248/max(enei), length( enei )*3 );
lorentz_wl_plot_only = 1248./lorentz_eV_plot_only;

sca(sca==0) = []; sca = reshape(sca,1,[]);
enei(enei==0) = []; enei = reshape(enei,1,[]);

%  photon setup for Lorentzian fit
enei_ev = 1248./enei ;
lorentz_eV = linspace( 1248/min(enei), 1248/max(enei), length( enei )*3 );
lorentz_wl = 1248./lorentz_eV;

[ max_sca, idx_max ] = max( real( sca(:) ) ) ; % return the index of wavelength where scattering is the biggest
[ min_sca, ~ ] = min( real( sca(:) ) ) ; % return the index of wavelength where scattering is the smallest
init_lw = 2/(pi*max_sca)*trapz( 1248./enei, real(sca(:)) ); % estimated value for Lorentzian fit, the linewidth. With this value, other parameters of lorentzian fit can be estimated with the maximum value of sca.

initial_guess = [min_sca init_lw/(2*pi) 1248/enei(idx_max) init_lw];  % 3rd is the SP - put expected value according to aspect ratio

options = optimoptions('lsqcurvefit','Display','off'); % do not show output of lsqcurvefit
[ fit_Lorentz, ~, ~, ~ ] = lsqcurvefit( @Lorentzfunction, initial_guess, enei_ev, real(sca(:)') , [], [], options );

fit_Lorentz(4) = abs(fit_Lorentz(4));  % sometimes lsqcurvefit finds a negative linewidth. This is the correction.

fit = fit_Lorentz(1) + fit_Lorentz(2)./( (lorentz_eV_plot_only - fit_Lorentz(3) ).^2 + ( 0.5*fit_Lorentz(4)).^2 ); % generate fitted curve on our own energy scale

%% FINAL PLOT

figure
plot( enei , real (sca) , 'o--' , lorentz_wl_plot_only , fit , 'r-' );
xlabel( 'Wavelength (nm)' );
switch type
    case 'abs'
        ylabel( 'Absorption cross section (nm^2)' );
        ndir = [ndir,'_abs'];
    otherwise
        ylabel( 'Scattering cross section (nm^2)' );
end
legend( ndir );
title( ['Lorentzian fitting - ',ndir] );
hold on
dim = [.6 .3 .5 .3];
str = {['\Gamma_{0}=',num2str(round(fit_Lorentz(4),3)), ' eV'],
    ['\Omega_{0}=',num2str(round(fit_Lorentz(3),3)), ' eV'],
    ['\Lambda_{0}=', num2str(round(1248./fit_Lorentz(3),2)), ' nm']};
t = annotation('textbox',dim,'String',str,'FitBoxToText','on');
t.FontSize=9;
legend('Data points','Lorentzian fit','Location','NorthEast')
set(gcf,'Visible','off')
set(gcf,'PaperPositionMode','auto','Color','white'); % maintain aspect ratio, background white
saveas (gcf, [ndir,'_Spectrum.fig'], 'fig')
saveas (gcf, [ndir,'_Spectrum.png'], 'png')


% save results
switch type
    case 'abs'
        fileID = fopen('abs.txt','w');
        fprintf(fileID,'%6s %6s \r\n','enei','abs');
        for j = 1 : length ( enei )
            fprintf(fileID,'%4.2f  %4.2f \r\n', enei (j), sca(j));
        end
        fclose(fileID);
        save ([ ndir,'_abs.mat'] , 'sca');
        %save ([ ndir,'_abs_fit.mat'] , 'fit');
    case 'sca'
        fileID = fopen('sca.txt','w');
        fprintf(fileID,'%6s %6s \r\n','enei','sca');
        for j = 1 : length ( enei )
            fprintf(fileID,'%4.2f  %4.2f \r\n', enei (j), sca(j));
        end
        fclose(fileID);
        save ([ ndir,'_sca.mat'] , 'sca');
        %save ([ ndir,'_sca_fit.mat'] , 'fit');
    case 'ext'
        fileID = fopen('ext.txt','w');
        fprintf(fileID,'%6s %6s \r\n','enei','txt');
        for j = 1 : length ( enei )
            fprintf(fileID,'%4.2f  %4.2f \r\n', enei (j), sca(j));
        end
        fclose(fileID);
        save ([ ndir,'_sca.mat'] , 'sca');
end

fileID = fopen([ndir,'_Lorentz.txt'],'w');
fprintf(fileID,'%10s %10s \r\n','lorentz_wl','fit');
for j = 1 : length (lorentz_wl)
    fprintf (fileID, '%6.2f % 6.2f \r\n', lorentz_wl (j), fit (j) );
end
fclose(fileID);

save ([ ndir,'_enei.mat'], 'enei');
save ([ ndir,'_params.mat'], 'fit_Lorentz');

end