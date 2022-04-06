% plot_fig = [0 0 0 0 1];
plot_fig = [0 0 0 0 1];

load Ie_top.mat
if 1%plot_fig(1)
  load IeFlickering-01.mat
%   theta_lims = acos(mu_lims)*180/pi;

end

E_max = size(Ie_top,3);
%%
load msis20051008_3.dat

h_atm = (74:1:610)'*1e3;h_atm = (100:1:610)'*1e3;
dz = @(n) 150 + 150/200*(0:(n-1))' +1.2*exp(((0:(n-1))-150)/17)';
z_atm = 100e3 + cumsum(dz(321)) - dz(1);
%%
load J.mat
figure
plot(t,[J_up;J_down])
%% Some CFL plots

nbr_subplot = 3;

figure
for iIe = 1:4
  
  filename = {['IeFlickering-0',num2str(iIe),'.mat']};
  load(filename{:})
  E_max = size(Ie_ztE,3);
   
  for i1 = 1:nbr_subplot  
    i11 = round(i1*size(Ie_ztE,2)/nbr_subplot); 

    I_up = 0;
    I_down = 0;
    E_up = 0;
    E_down = 0;
    for iZ = 0:320
      Ie_at_one_z = Ie_ztE(321-iZ:321:10*321,:,:);
      E_at_one_z_up = 0;
      E_at_one_z_down = 0;
      I_at_one_z_up = 0;
      I_at_one_z_down = 0;
      for i2 = 6:10   % upward fluxes
        E_at_one_z_up = E_at_one_z_up + ...
                        trapz(E(1:E_max),E(1:E_max).*squeeze(Ie_at_one_z(i2,i11,:)).');
        I_at_one_z_up = I_at_one_z_up + ...
                        trapz(E(1:E_max),squeeze(Ie_at_one_z(i2,i11,:)).');
      end
      for i2 = 1:5    % downward fluxes
        E_at_one_z_down = E_at_one_z_down + ...
                        trapz(E(1:E_max),E(1:E_max).*squeeze(Ie_at_one_z(i2,i11,:)).');
        I_at_one_z_down = I_at_one_z_down + ...
                        trapz(E(1:E_max),squeeze(Ie_at_one_z(i2,i11,:)).');
      end
      
      I_up(iZ+1) = I_at_one_z_up;
      I_down(iZ+1) = I_at_one_z_down;
      E_up(iZ+1) = E_at_one_z_up;
      E_down(iZ+1) = E_at_one_z_down;
    end
    
  subplot(2,6,(iIe-1)*nbr_subplot+i1)
%   plot(E_up,z_atm([end:-1:1])/1e3)
  semilogx(I_up,z_atm([end:-1:1])/1e3)
  hold on
%   plot(E_down,z_atm([end:-1:1])/1e3)
  semilogx(I_down,z_atm([end:-1:1])/1e3)
%   xlim([1e15 1e18]);
  xlim([1e12 1e15]);
  title(['t =', num2str(t_run(i11)),' s'])
  hold off
  
  end
end


%% Correct the error in Ie_top
% Ie_top has been converted in electron/ster/eV/m^2/s in make_all_Ie_top.m
% script. But the conversion in per eV has been done with the wrong indices
% (only last 10 energy bins have been scaled)
for i3 = size(Ie_top,1):-1:1,
  % converting back the last 10 energy bins to initial values (i.e. not per
  % eV)
  Ie_top(:,:,i3) = Ie_top(:,:,i3)*dE(i3);
end
for i3 = size(Ie_top,3):-1:1,
  % and now doing it properly!
  Ie_top(:,:,i3) = Ie_top(:,:,i3)/dE(i3);
end
%% Fig 1
% Make a subplot of the e- fluxes for the different upward streams. Plotted
% as a heatmap over time and energy
Ie_plot = Ie_top_raw;

if plot_fig(1)
  figure(1)
  tiledlayout(2,5)
  for i_UP = size(Ie_top,1):-1:1  %size(Ie_top,1)/2+1
    nexttile
    pcolor(t,E(1:E_max),log10(max(0,real(squeeze(Ie_plot(i_UP,:,:))))).')
    shading flat
    xlabel('time (s)')
    ylabel('energy (eV)')
    title([num2str(theta_lims(i_UP+1)),' - ',num2str(theta_lims(i_UP)), '°'])
    caxis([3 14])
    colormap(jet)
    set(gca,'YScale','log');
  end
  cb = colorbar;
  cb.Layout.Tile = 'east';
  cb.Label.String = 'log10(#/m2/s/eV/ster)';
end
%% Fig 2
% Plot the total upward e- flux (sum over the different upward streams)
% at a specific time (time can be chosen below as t_to_plot)
t_to_plot = 1; % (s)
Ie_plot = Ie_top_raw;

Ie_plot_UP = 0;
%loop over upward streams
for i_UP = size(Ie_plot,1):-1:size(Ie_plot,1)/2+1
  Ie_plot_UP = Ie_plot_UP + Ie_plot(i_UP,:,:);
end

Ie_plot_UP = squeeze(Ie_plot_UP);
[~,index_t] = min(abs(t-t_to_plot));

if plot_fig(2)
  figure(2)
  semilogy(E(1:E_max),squeeze(Ie_plot_UP(index_t,:)));
  xlabel('E (eV)')
  ylabel('Flux (#/s/m²/eV/ster)')
  title(['Upward electron flux @ t = ',num2str(t(index_t)),'s'])
end
%% Fig 3
% Plot the total upward e- flux (sum over the different upward streams)
% integrated over time, as a function of E.
Ie_plot = Ie_top_raw;

Ie_plot_UP = 0;
%loop over upward streams
for i_UP = size(Ie_plot,1):-1:size(Ie_plot,1)/2+1
  Ie_plot_UP = Ie_plot_UP + Ie_plot(i_UP,:,:);
end
Ie_plot_UP = squeeze(Ie_plot_UP);
dt = t(2) - t(1);
%integrate over time
Ie_plot_UP = sum(Ie_plot_UP*dt,1);

if plot_fig(3)
  figure(3)
  semilogy(E(1:E_max),squeeze(Ie_plot_UP));
  xlabel('E (eV)')
  ylabel('Flux (#/m²/eV/ster)')
  title('total Upward electron flux integrated over time')
end
%% Fig 4
% Compute the total upward e- flux and total downward e- flux, integrate
% them over time, and plot the ratio of the former (up) over the latter
% (down)
% c_o_mu = mu_avg(mu_lims);
c_o_mu = ones(1,10);
% Ie_plot = Ie_ztE(321:321:3210,:,:) .* abs(c_o_mu).';
Ie_plot = Ie_top_raw .* abs(c_o_mu).';

% dt = diff(t); dt = [dt dt(end)];
dt = ones(1,11);
% dt = dt(1:6) 

Ie_plot_UP = 0;
Ie_plot_DOWN = 0;
%loop over upward streams
for i_UP = size(Ie_plot,1):-1:size(Ie_plot,1)/2+1
  Ie_plot_UP = Ie_plot_UP + Ie_plot(i_UP,:,:);
end
Ie_plot_UP = squeeze(Ie_plot_UP);
% Ie_plot_UP = Ie_plot_UP(1:101,:);
Ie_plot_UP = dt * Ie_plot_UP;

%loop over downward streams
for i_DOWN = 1:size(Ie_top,1)/2
  Ie_plot_DOWN = Ie_plot_DOWN + Ie_plot(i_DOWN,:,:);
end
Ie_plot_DOWN = squeeze(Ie_plot_DOWN);
% Ie_plot_DOWN = Ie_plot_DOWN(1:101,:);
Ie_plot_DOWN = dt*Ie_plot_DOWN;

if plot_fig(4)
  figure(4)
  plot(E(1:E_max),Ie_plot_UP./Ie_plot_DOWN)
  xlabel('E (eV)')
  title('Ratio of upward over downward flux of e^-')
end
% AND DISP TOTAL ENERGY FLUX RATIO

Ie_at_one_z = Ie_ztE(321:321:3210,:,:) .* abs(c_o_mu).';
E_flux_up = 0
E_flux_down = 0
  for i1 = size(Ie_at_one_z,2)
    for i2 = 6:10
      E_flux_up = E_flux_up + ...
                  trapz(E(1:E_max),E(1:E_max).*squeeze(Ie_at_one_z(i2,i1,:)).');
    end
    for i2 = 1:5
      E_flux_down = E_flux_down + ...
                    trapz(E(1:E_max),E(1:E_max).*squeeze(Ie_at_one_z(i2,i1,:)).');
    end
  end
  E_flux_up
  E_flux_down
disp(['Energy back = ', num2str(E_flux_up/E_flux_down*100), '%']);

Energy_flux_UP = Ie_plot_UP * E(1:E_max).'
Energy_flux_DOWN = Ie_plot_DOWN * E(1:E_max).'

disp(['Energy back = ', num2str(Energy_flux_UP/Energy_flux_DOWN*100), '%']);
%% Fig 5 & 6
% Plot the total upward flux integrated over energies, as a function of t
% if ~exist('Ie_plot_UP','var') 

% c_o_mu = mu_avg(mu_lims);
c_o_mu = ones(1,10);
Ie_plot = Ie_top_raw .* abs(c_o_mu).';


Ie_plot_UP = 0;
Ie_plot_DOWN = 0;
%loop over upward streams
for i_UP = size(Ie_plot,1):-1:size(Ie_plot,1)/2+1
  Ie_plot_UP = Ie_plot_UP + Ie_plot(i_UP,:,:);
end
Ie_plot_UP = squeeze(Ie_plot_UP);
%loop over downward streams
for i_DOWN = 1:size(Ie_plot,1)/2
  Ie_plot_DOWN = Ie_plot_DOWN + Ie_plot(i_DOWN,:,:);
end
Ie_plot_DOWN = squeeze(Ie_plot_DOWN);

if plot_fig(5)
  figure(5)

  subplot(3,2,1)
  plot(t,sum(Ie_plot_UP,2))
  ylabel('Ie^{UP} (#e/m2/s)')
  title('Flux e^- UP')

  subplot(3,2,2)
  plot(t,Ie_plot_UP*E(1:E_max).')
  ylabel('Ie^{UP} (eV/m2/s)')
  title('Flux energy UP')

  subplot(3,2,3)
  plot(t,sum(Ie_plot_DOWN,2),'r')
  ylabel('Ie^{DOWN} (#e/m2/s)')
  title('Flux e^- DOWN')

  subplot(3,2,4)
  plot(t,Ie_plot_DOWN*E(1:E_max).','r')
  ylabel('Ie^{DOWN} (eV/m2/s)')
  title('Flux energy DOWN')

  subplot(3,2,5)
  plot(t,sum(Ie_plot_UP,2))
  hold on
  plot(t,sum(Ie_plot_DOWN,2))
  hold off
  xlabel('t (s)')
  ylabel('Ie (#e/m2/s)')
  title('Flux e^-')
  legend('UP', 'DOWN','location','northwest')

  subplot(3,2,6)
  plot(t,Ie_plot_UP*E(1:E_max).')
  hold on
  plot(t,Ie_plot_DOWN*E(1:E_max).')
  hold off
  xlabel('t (s)')
  ylabel('Ie (eV/m2/s)')
  title('Flux energy')
  legend('UP', 'DOWN','location','northwest')
end




figure(6)
plot(t,cumsum(Ie_plot_UP * E(1:E_max).' .* dt.'),'linewidth',2)
hold on
plot(t,cumsum(Ie_plot_DOWN*E(1:E_max).' .* dt.'),'linewidth',2)
hold off
xlabel('t (s)')
ylabel('Energy (eV)')
title('Total of the energy passing through the top of the ionosphere')
legend('UP', 'DOWN','location','northwest')
%% Fig 7
% some test
if plot_fig(7)

  Ie_plot = Ie_top_raw;
  Ie_plot_DOWN = 0;
  %loop over downward streams
  for i_DOWN = 1:size(Ie_plot,1)/2
    Ie_plot_DOWN = Ie_plot_DOWN + Ie_plot(i_DOWN,:,:);
  end
  Ie_plot_DOWN = squeeze(Ie_plot_DOWN);
  
  Ie_plot = I_TEST_TOP_SAVE_TOT;
  Ie_plot_DOWN2 = 0;
  %loop over downward streams
  for i_DOWN = 1:size(Ie_plot,1)/2
    Ie_plot_DOWN2 = Ie_plot_DOWN2 + Ie_plot(i_DOWN,:,:);
  end
  Ie_plot_DOWN2 = squeeze(Ie_plot_DOWN2);
  
  
  figure(7)

  subplot(3,2,1)
  plot(t,sum(Ie_plot_DOWN2,2))
  ylabel('Ie^{UP} (#e/m2/s)')
  title('Flux e^- DOWN raw')

  subplot(3,2,2)
  plot(t,Ie_plot_DOWN2*E(1:E_max).')
  ylabel('Ie^{UP} (eV/m2/s)')
  title('Flux energy DOWN raw')

  subplot(3,2,3)
  plot(t,sum(Ie_plot_DOWN,2))
  ylabel('Ie^{DOWN} (#e/m2/s)')
  title('Flux e^- DOWN')

  subplot(3,2,4)
  plot(t,Ie_plot_DOWN*E(1:E_max).')
  ylabel('Ie^{DOWN} (eV/m2/s)')
  title('Flux energy DOWN')

  subplot(3,2,5)
  plot(t,sum(Ie_plot_DOWN2,2))
  hold on
  plot(t,sum(Ie_plot_DOWN,2))
  hold off
  xlabel('t (s)')
  ylabel('Ie (#e/m2/s)')
  title('Flux e^-')
  legend('DOWN raw', 'DOWN','location','northwest')

  subplot(3,2,6)
  plot(t,Ie_plot_DOWN2*E(1:E_max).')
  hold on
  plot(t,Ie_plot_DOWN*E(1:E_max).')
  hold off
  xlabel('t (s)')
  ylabel('Ie (eV/m2/s)')
  title('Flux energy')
  legend('DAWN raw', 'DOWN','location','northwest')
end
%% Comparison of energies up/down

AWA_10_top = [24 28 34]; AWA_20_top = [56 71 96];
AWA_10_top_raw = [104 141 214]; AWA_20_top_raw = [96 138 220];

figure(8)
plot([1 2 4],AWA_10_top,'b+-')
hold on
plot([1 2 4],AWA_20_top,'bo-')
plot([1 2 4],AWA_10_top_raw,'r+-')
plot([1 2 4],AWA_20_top_raw,'ro-')
hold off
xlim([0 5])
ylim([0 inf])
legend('AWA-10-top','AWA-20-top','AWA-10-top-raw','AWA-20-top-raw', ...
  'location','northwest')
xlabel('x E')
ylabel('E back in %')



%
PnLh_3keV_top = [37 67 97];
PnLh_3keV_topraw = [440 460 480];
PnLh_5keV_top = 61;
PnLh_5keV_topraw = 753;
PnLh_7keV_top = 80;
PnLh_7keV_topraw = 996;

figure(9)
plot([10 30 60],PnLh_3keV_top,'b+-')
hold on
plot([10 30 60],PnLh_3keV_topraw,'r+-')
plot(10,PnLh_5keV_top,'bo')
plot(10,PnLh_5keV_topraw,'ro')
plot(10,PnLh_7keV_top,'b^')
plot(10,PnLh_7keV_topraw,'r^')
% set(gca,'YScale','log');
hold off
xlim([0 90])
legend('PnLh-3keV-top','PnLh-3keV-topraw', ...
  'PnLh-5keV-top','PnLh-5keV-topraw', ...
  'PnLh-7keV-top','PnLh-7keV-topraw',...
  'location','northeast')
xlabel('pitch angle width of precipitation')
ylabel('E back in %')



%
G_10_top = [30 65 94 108];
G_10_topraw = [201 749 1157 1359];
G_90_top = [311 717 836];
G_90_topraw = [176 691 906];

figure(10)
plot([1 3 5 7],G_10_top,'b+-')
hold on
plot([1 3 5 7],G_10_topraw,'r+-')
plot([1 3 5],G_90_top,'bo-')
plot([1 3 5],G_90_topraw,'ro-')
hold off
xlim([0 8])
xlabel('energy Gaussian (keV)')
ylabel('E back in %')
legend('G-10-top','G-10-topraw', ...
  'G-90-top','G-90-topraw', ...
  'location','northwest')




