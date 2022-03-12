function [Ie_zE,mu_pars_out] = Ie_M_stream_4_auroraB(h,mag_ze,B,E,mu_lims,mu_pars,Ie_primary,p_e_q,ne,Te,varargin)
% Ie_M_stream_4_aurora - steady-state multi-stream electron transport
% 
% Calling:  [Ie_zE,mu_pars_out] = Ie_M_stream_4_aurora(h,mag_ze,E,mu_lims,mu_pars,Ie_primary,p_e_q,ne,Te,...
%                                     nO,O_levels,XsO,@O_e_2nd_dist,@phase_fcn_O,...
%                                     nN2,N2_levels,XsN2,@N2_e_2nd_dist,@phase_fcn_N2,...
%                                     % Additional rows for each
%                                     % sufficiently abundant species
%                                     ni,E_Li,Xsi,@i_e_2nd_dist,@phase_fcn_i);
% Input:
%  h       - altitude (m), double array [nZ x 1]
%  mag_ze  - magnetic zenith angle (radians), double scalar
%  B       - Magnetic field or magnetic field-strength array. 
%  E       - energy grid (eV), double array [1 x nE]
%  mu_lims - cosine-of-pitch-angle limits of beams, double array
%            [1 x (n_beams+1)], should start and stop with -1 and
%            +1 for field-aligned down and up respectively
%  mu_pars - cell-array with the output from e_scattering_beamdistribution,
%            that is the three-dimensional array with probabilities
%            for scattering into each beam and the mateix with the
%            solid angles for each pitch-angle in each stream:
%            {Pmu2mup,theta2beamW}, if not
%            e_scattering_beamdistribution will be calculated with
%            181 steps in angles (which seems to be a bit on the
%            thin side)
%  Ie_primary - cell-array with function-handles to functions for
%            the electron fluxes in each stream, for example:
%            I_primary{1} = @(t,E) Ie_smooth4(t,E). The function
%            should return the electron flux in an energy bin for
%            the correspoding stream for the full array of times
%            energy-by-energy.
%  p_e_q   - internal source of energetic electrons
%            (m^-3s^-1 dE^-1), double array, either [nZ x nE] (in
%            which case the production is distributed isotropically
%            over the pitch-angle-streams) or [(nZ*n_beams) x nE]
%            (in which the production is take to be per
%            pitch-angle-stream, sorted like the electron fluxes)
%  ne      - electron concentration (m^-3), double array [nZ x 1]
%  Te      - electron temperature (K), double array [nZ x 1]
%
%  ni      - number density (m^-3) of species with collision cross
%            sections, double array [nZ x 1]
%  E_Li    - excitation thresholds (eV) of states in ni
%            (dE_el,dE_i,E_ion), with number of secondary electrons
%            produced for each excitation in the second column,
%            double array [nLevels x 2]
%  Xsi     - collision cross sections (m^2) (Xs_el;Xs_i;Xs_ion)
%            double array [nLevels x nE]
%  second_e_E_fcn_i - function for spectrum of secondary electrons
%  phase_fcn_i - phase-functions for angular scattering as a
%                function of energy and angles, function-handle
%                that return the elastic inelastic angular
%                scattering probability
%  
%  The number of species that the electrons collides with is
%  arbitrary. It is only required that N_i, E_Lij, Xs_ij,
%  second_e_E_fcn_i(E) and phase_fcn are provided for all. 
% 
% SEE also: e_scattering_beamdistribution, Etrp_*.m

%  Copyright © Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
%  This is free software, licensed under GNU GPL version 2 or later

% wbh = waitbar(1,'Working in the coal mine...');

second_e_isotropic = 1;
nZ = numel(h);

% The functions calculates the electron flux in energy-ranges from
% E_{i} to E_{i+1}, here we calculate dE
gradE = diff(E);
gradE = gradE([1:end end]);

if iscell(mu_pars)
  n_dirs = size(mu_pars{1},2);
else
  n_dirs = 181;
end
theta = linspace(0,pi,n_dirs);

% Array for degrading electrons
deg_Ie = zeros(numel(h)*numel(mu_lims(1:end-1)),numel(E));



Ie_zE = zeros(numel(h)*numel(mu_lims(1:end-1)),numel(E));

%% 2 Calculate the altitude-time-variation of electron fluxes down
%  and up through the thermosphere, starting at the highest
%  energy and iterate down in energy. For each energy calculate the
%  excitation and ionisation of all states in O, N2 and O2 and the
%  corresponding energy-degradation and production of secondary
%  electrons.
if iscell(mu_pars) && numel(size(mu_pars{1})) == 3
  Pmu2mup     = mu_pars{1};
  theta2beamW = mu_pars{2};
  mu_pars_out = mu_pars;
else
  [Pmu2mup,theta2beamW] = e_scattering_beamdistribution(mu_lims,n_dirs);
  mu_pars_out = {Pmu2mup,theta2beamW};
end
c_o_mu = mu_avg(max(-1,min(1,mu_lims)));
BeamWeights = sum(theta2beamW,2);
BeamWeights(c_o_mu<0) = BeamWeights(c_o_mu<0)/sum(BeamWeights(c_o_mu<0));
BeamWeights(c_o_mu>0) = BeamWeights(c_o_mu>0)/sum(BeamWeights(c_o_mu>0));

[Mi2ip1,Mi2im1] = dI_mirror_dthetadB(acos(mu_lims),B);

elsc_b2b = zeros(numel(h),size(Pmu2mup,3),size(Pmu2mup,3));
for iE = length(E):-1:1,
  
  elsc_b2b = 0*elsc_b2b;
  try
    waitbar(iE/length(E),wbh);
  catch
    % Whatever, it becomes tedious...
  end
  %% 1: calculate the total and elastic scattering cross-section with
  %  accounting for the energy-dependent back-scattering. This does
  %  not depend on time, only on altitude, through the neutral
  %  density variation, and energy through the cross-sections.
  
  % Elastic
  % and inelastic collision cross sections
  tot_sc = zeros([length(h) 1]);
  % Loop over all species to calculate total and elastic cross
  % sections, and the corresponding back-scattering ratio:
  for i1 = 1:5:length(varargin),
    
    n_i = varargin{i1};   % Neutral density
    xs_i = varargin{i1+2};% Array with collision cross sections
    dE_i = varargin{i1+1};% Energy levels
                          % second_e_fcn = varargin{i1+3}; % Energy spectra for secondary electrons
    ang_scat_fcn = varargin{i1+4}; % Back scattering ratio (energy
                                   % dependent)
                                   % add the elastic back scattering cross sections together
    [curr_phfc_e,curr_phfc_i] = ang_scat_fcn(E(iE),theta(:));
    curr_phfc_e = curr_phfc_e.*sin(theta');
    curr_phfc_e = curr_phfc_e/sum(curr_phfc_e);
    B2B = beams2beams(curr_phfc_e,Pmu2mup,theta2beamW);
    curr_phfc_i = curr_phfc_i.*sin(theta');
    curr_phfc_i = curr_phfc_i/sum(curr_phfc_i);
    A2B = beams2beams(curr_phfc_i,Pmu2mup,theta2beamW);
    for i_1 = size(B2B,1):-1:1,
      for i_2 = size(B2B,2):-1:1
        elsc_b2b(:,i_1,i_2) = elsc_b2b(:,i_1,i_2) + n_i*(xs_i(1,iE)*B2B(i_1,i_2));
      end
    end
    % add the elastic back scattering to the total
    tot_sc = tot_sc + n_i*xs_i(1,iE);
    % add the inelastic collisions to the total
    for j2 = 2:size(xs_i,1)-1,
      
      % The second factor corrects for the case where the energy loss
      % dE_i(j2,1) is smaller than the width in energy in the energy bin
      tot_sc = tot_sc + n_i*(xs_i(j2,iE)); %.*min(1,dE_i(j2,1)./gradE(iE)));
      for i_1 = size(B2B,1):-1:1,
        for i_2 = size(B2B,2):-1:1
          elsc_b2b(:,i_1,i_2) = ( elsc_b2b(:,i_1,i_2) + ...
                                  n_i*(xs_i(j2,iE)*A2B(i_1,i_2)).*max(0,1-dE_i(j2,1)./gradE(iE)));
        end
      end
      % That is when gradE(iE) > dE_i(j2,1) only the fraction
      % dE_i(j2,1)/gradE is lost from the energy bin [E(iE) E(iE)+gradE(iE)]
      
    end
    % Add the ionisation to the total scattering. This requires
    % special treatment since in addition to the ionisation potential
    % primary electrons lose energy to the secondary electrons
    j2 = size(xs_i,1);
    % So we calculate the secondary electron spectra
    % second_e_spectra = feval(second_e_fcn,E,E(iE));
    % second_e_spectra = (second_e_spectra + second_e_spectra([2:end end]))/2.*gradE;
    % And the average energy lost to them
    % dE_s = cumsum(second_e_spectra.*E)./cumsum(second_e_spectra);
    % This makes the fraction of electrons lost from the energy bin
    % [E(iE) E(iE)+gradE(iE)]
    % tot_sc = tot_sc + n_i*(xs_i(j2,iE).*min(1,(dE_i(j2,1)+dE_s(iE))./gradE(iE)));
    %% This is a QD-fix since this is not the correct
    %  phase-function for the ionising collisions, but this will
    %  have to do for now.
    tot_sc = tot_sc + n_i*(xs_i(j2,iE)); %.*min(1,dE_i(j2,1)./gradE(iE)));
    for i_1 = size(B2B,1):-1:1,
      for i_2 = size(B2B,2):-1:1
        elsc_b2b(:,i_1,i_2) = ( elsc_b2b(:,i_1,i_2) + ...
                                n_i*(xs_i(j2,iE)*A2B(i_1,i_2)).*max(0,1-dE_i(j2,1)./gradE(iE)));
      end
    end
  end
  
  % Add the loss due to electron-electron collisions:
  A = tot_sc + dEds_ee(E(iE),ne,Te)/gradE(iE);
  
  if size(p_e_q,1) == numel(h)
    % redistribute isotropically
    for i_mu = numel(c_o_mu):-1:1
      if c_o_mu < 0
        W = BeamWeights(i_mu)/sum(BeamWeights(c_o_mu<0));
      else
        W = BeamWeights(i_mu)/sum(BeamWeights(c_o_mu>0));
      end
      p_e_Q((1:numel(h))+(i_mu-1)*numel(h)) = p_e_q(:,iE)*W;
    end
  else % p_e_q should be per pitch-angle-stream
    p_e_Q = p_e_q(:,iE);
  end
  for i_mu = 1:numel(c_o_mu),
    Ie_p{i_mu} = Ie_primary{i_mu}(E(iE))*gradE(iE);
  end
  DE = gradE([iE,min(numel(gradE),iE+1)]);
  
  % This is the call to the solver of the coupled difference
  % equations, with up-stream differences for the spatial gradients
  [Ie_ziE] = de_M_streamB_us(h/abs(cos(mag_ze)),c_o_mu,...
                             Ie_p,...
                             deg_Ie(:,iE)+(p_e_Q(:)),...
                             A,elsc_b2b,...
                             Mi2im1,Mi2ip1, ...
                             0);
  % TODO: Check if this is often applied outside of at
  % too-low-altitudes:
  Ie_ziE(Ie_ziE(:)<0) = 0;
  
  Ie_zE(:,iE) = Ie_ziE;
  
  %% Calculation of the energy degradation and production of
  %  secondary electrons, starting off with adding the energy
  %  degradation due to electron-electron collisions:
  if iE > 1
    deg_Ie(:,iE-1) = deg_Ie(:,iE-1) + ...
        repmat(dEds_ee(E(iE),ne,Te)/gradE(iE),size(c_o_mu(:),1),1).*Ie_zE(:,iE);
  end 
  % then continuing doing species by species
  for i1 = 1:5:length(varargin),
    
    n_i = varargin{i1};   % Neutral density
    xs_i = varargin{i1+2};% Array with collision cross sections
    dE_i = varargin{i1+1};% Energy levels
    second_e_fcn = varargin{i1+3}; % Energy spectra for secondary electrons
    ang_scat_fcn = varargin{i1+4}; % Back scattering ratio (energy
                                   % dependent)
                                   % add the elastic back scattering cross sections together
    [~,curr_phfc_i] = ang_scat_fcn(E(iE),theta(:));
    curr_phfc_i = curr_phfc_i.*sin(theta');
    curr_phfc_i = curr_phfc_i/sum(curr_phfc_i);
    B2B = beams2beams(curr_phfc_i,Pmu2mup,theta2beamW);
    idx1 = [];
    idx2 = [];
    Ab2b = [];
    for i_1 = size(B2B,1):-1:1,
      for i_2 = size(B2B,2):-1:1,
        idx1 = [idx1,(1:numel(h))+(i_1-1)*numel(h)];
        idx2 = [idx2,(1:numel(h))+(i_2-1)*numel(h)];
        Ab2b = [Ab2b, n_i*B2B(i_1,i_2)];
      end
    end
    A_b2b = sparse(idx1,idx2,Ab2b);
    
    % Now loop over the inelastic collision cross sections (except
    % ionisation):
    for idE = 2:size(dE_i,1),
      if dE_i(idE,2) > 0
        break;
        % then we're off to ionisation's
      end
      % The flux of electrons degraded from energy bin [E(iE) E(iE)+gradE(iE)]
      Ie_degraded = (xs_i(idE,iE).*min(1,dE_i(idE,1)./gradE(iE)))*A_b2b*Ie_ziE;
      % to any lower energy bin by excitation of the idE-th state
      % of the current species.
      
      % Find the energy bins that the electrons in the current
      % energy bin will degrade to when losing dE_i(idE,1) eV
      i_upper = find( E + gradE > E(iE)-dE_i(idE,1) & ...
                      E < E(iE)+gradE(iE)-dE_i(idE,1) );
      partition_fract = zeros(size(i_upper));
      if ( ~isempty(i_upper) && i_upper(1) < iE ) 
        % Distribute the degrading electrons between those bins
        % (don't ask, figure it out for yourself if you really want
        % to know... ...I have sadly forgotten)
        partition_fract(1) = min(1,...
                                 (E(i_upper(1))+gradE(i_upper(1)) - ...
                                  (E(iE)-dE_i(idE,1)))/gradE(iE));
        if length(i_upper)>2
          partition_fract(2:end-1) = min(1,gradE(i_upper(2:end-1))/gradE(iE));
        end
        partition_fract(end) =  min(1,...
                                    (E(iE)+gradE(iE)-dE_i(idE,1) - ...
                                     E(i_upper(end)))/gradE(iE));
        if i_upper(end) == iE
          partition_fract(end) = 0;
        end
        partition_fract = partition_fract/sum(partition_fract);
        
        for i_u = find(partition_fract~=0),%1:numel(i_upper),
          deg_Ie(:,i_upper(i_u)) = ( deg_Ie(:,i_upper(i_u)) + ...
                                     max(0,Ie_degraded) * ...
                                     partition_fract(i_u) );
        end
      end
      
    end
    % What remains is the ionisation
    for idE2 = idE:size(dE_i,1),
      
      E_p_d = -( E - ( E(iE)+gradE(iE) - dE_i(idE2,1) ) );
      % Then take care of the ionisation
      i_upper = find( E + gradE > E(iE)-dE_i(idE2,1) & ...
                      E < E(iE)+gradE(iE)-dE_i(idE2,1) );% & ...
      if ( ~isempty(i_upper) && i_upper(1) < iE )
        % Calculate the ionisation rate:
        if second_e_isotropic % TODO: Insert proper 2nd-ary e
                              % production rates, direction for
                              % direction just as for elastic and
                              % inelastic collisions - except here
                              % we have to do it energy-by-energy/BG-20180921
          Ionization = 0*Ie_ziE;
          % Altitude-profiles of total ionizations
          Ionizing = repmat(n_i,numel(c_o_mu),1).*(xs_i(idE2,iE) * Ie_ziE);
          for i_mu = numel(c_o_mu):-1:1,
            % Altitude-profiles of total ionizations in each
            % pitch-angle-stream
            Ionization((1:nZ)+nZ*(i_mu-1),:) = max(0,...
                                                   (n_i).*(xs_i(idE2,iE) * Ie_ziE((1:nZ)+nZ*(i_mu-1)))*...
                                                   BeamWeights(i_mu)/sum(BeamWeights));
          end
        else
          Ionizing   = max(0,...
                           repmat(n_i,numel(c_o_mu),1).*(xs_i(idE2,iE) * Ie_ziE));
          Ionization = max(0,...
                           repmat(n_i,numel(c_o_mu),1).*(xs_i(idE2,iE) * Ie_ziE));
        end
        
        % Calculate the spectra of the secondary electrons, using the
        % average energy of the electrons in the current energy bin:
        second_e_spectra = second_e_fcn(E,E(iE),dE_i(idE2,1),'s');
        second_e_spectra = (second_e_spectra + second_e_spectra([2:end end]))/2.*gradE;
        % And the distribution of the ionising electrons (that have
        % to lose the corresponding amount of energy)
        deg_sp = second_e_fcn(E,E(iE),dE_i(idE2,1),'c');
        
        if sum(second_e_spectra) > 0
          % Multiply with the number of electrons created - double
          % ionization, double dissociative ionization gives two...
          second_e_spectra = dE_i(idE2,2)*second_e_spectra/(sum(second_e_spectra));
        end
        if sum(deg_sp) > 0
          deg_sp = deg_sp/sum(deg_sp);
        end
        % TODO: Check that the sum of s_e_s and deg_Sp are the same
        e_ionized_distribution = max(second_e_spectra,deg_sp); % or just sum
        e_ionized_distribution(~isfinite(e_ionized_distribution)) = 0;
        
        
        % Add the ionisation energy and the energy of the secondary
        % electrons to get the degradation energy of the ionising
        % electrons 
        
        
        % and add this to the degrading electron fluxes, here there
        % should be "no" angular-scattering:
        % This part should be possible to move out of the loop over
        % the species
        i_l_s_e_s = find(second_e_spectra>0,1,'last'); 
        deg_Ie(:,1:i_l_s_e_s) = deg_Ie(:,1:i_l_s_e_s) + Ionization*second_e_spectra(1:i_l_s_e_s);
        i_f_ds = find(deg_sp>0,1,'first'); 
        deg_Ie(:,i_f_ds:(iE-1)) = deg_Ie(:,i_f_ds:(iE-1)) + Ionizing*deg_sp(i_f_ds:(iE-1));
        
      end
    end
    
  end
  try
    load('stopnow.dat')
    if stopnow == 1 || stopnow==iE
      keyboard
    end
  catch
  end
  
end

try
  close(wbh)
catch
end

