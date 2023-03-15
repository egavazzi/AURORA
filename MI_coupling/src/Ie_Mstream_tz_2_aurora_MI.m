function [Ie_ztE,mu_pars_out] = Ie_Mstream_tz_2_aurora_MI(AURORA_root_directory,n_loop,h,mag_ze,E,mu_lims,mu_pars,t,I0,p_e_q,ne,Te,OPS,varargin)
    % Ie_Mstream_tz_2_aurora - time-dependent multi-stream electron transport
    % 
    % Calling: 
    %  [Ie_ztEG,mu_pars_out] = Ie_Mstream_tz_2_aurora(h,mag_ze,E,mu_lims,mu_scatterings,t,...
    %                                     I0,Ie_primary,p_e_q,ne,Te,...
    %                                     OPS
    %                                     ni,E_Li,Xsi,@i_e_2nd_dist,@phase_fcn_i,...
    %                                     nO,O_levels,XsO,@O_e_2nd_dist,@phase_fcn_O,...
    %                                     nN2,N2_levels,XsN2,@N2_e_2nd_dist,@phase_fcn_N2,...
    %                                     % Additional rows for each
    %                                     % sufficiently abundant species
    %                                     nO2,O2_levels,XsO2,@O2_e_2nd_dist,@phase_fcn_O2);
    % Input:
    %  h       - altitude (m), double array [nZ x 1]
    %  mag_ze  - magnetic zenith angle (radians), double scalar
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
    %  t       - time array (s), double array [1 x n_t]
    %  I0      - electron flux profile at the start, i.e. initial
    %            condition (#e/m^2/s/dE), double array [(nZ x n_beams) x 1]
    %            with the streams stacked in order from most parallel
    %            to B downward to most paralell to B upwards.
    %  p_e_q   - internal source of energetic electrons
    %            (m^-3s^-1 dE^-1), double array [nZ x nE]
    %  ne      - electron concentration (m^-3), double array [nZ x 1]
    %  Te      - electron temperature (K), double array [nZ x 1]
    %  OPS     - options struct, optional argument, with parameters
    %            controlling the electron transport. Fields used:
    %            OPS.second_e_isotropic - 1 for isotropic production of
    %            secodary electrons, 0 for forward production.
    %            OPS.CscD_e - scaling factor for time-of-flight
    %            spreading/diffusion, ought to be between 0 and 1, zero
    %            disables time-of-flight smearing due to spread of
    %            field-aligned velocities within the energy-pitch-angle
    %            cells.
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
    
    %  Copyright � Bjorn Gustavsson 20180603, bjorn.gustavsson@uit.no
    %  This is free software, licensed under GNU GPL version 2 or later
    
    %wbh = waitbar(1,'Working in the coal mine...');
    
    dOPS.second_e_isotropic = 1;
    dOPS.CscD_e = 1;
    
    if ~isempty(OPS) && isstruct(OPS)
      dOPS = merge_structs(dOPS,OPS);
    end
    
    second_e_isotropic = dOPS.second_e_isotropic;
    CscD_e = dOPS.CscD_e
    nZ = numel(h);
    
    % The functions calculates the electron flux in energy-ranges from
    % E_{i} to E_{i+1}, here we calculate dE
    gradE = diff(E);
    gradE = gradE([1:end end]);
    
    % average pitch-angle cosine
    %beta_m = 1/3^.5; 
    
    if iscell(mu_pars)
      n_dirs = size(mu_pars{1},2);
    else
      n_dirs = 181;
    end
    theta = linspace(0,pi,n_dirs);
    
    % Array for degrading electrons
    deg_Ie = zeros(numel(h)*numel(mu_lims(1:end-1)),numel(t),numel(E));
    
    
    % electron electron energy loss function...
    %L = esup_e_loss(E,ne,Te);
    %[dLdE,~] = gradient(L,E,h);
    
    % ...multiplied with electron density
    %neL = repmat(ne,size(E)).*L;
    
    Ie_ztE = zeros(numel(h)*numel(mu_lims(1:end-1)),numel(t),numel(E));
    Ie_ztEpdE = zeros(numel(h)*numel(mu_lims(1:end-1)),numel(t));
    
    % AX2 = -inf;
    
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
    
    
    % Here the normalization factor is calculated anew with theta2beamW
    % instead of using the already calculated BeamW. This is to ensure
    % conservation of the number of electrons during scattering as the
    % calculations are discretized in pitch angle, whereas BeamW was 
    % calculated in a more continuous way (integral of sin(theta)). 
    % This leads to two slightly different normalization factors, the one
    % calculated with theta2beamW being the one to use here if we don't want 
    % to lose some electrons.
    BeamWeights = sum(theta2beamW,2);
    BeamWeights(c_o_mu<0) = BeamWeights(c_o_mu<0)/sum(BeamWeights(c_o_mu<0));
    BeamWeights(c_o_mu>0) = BeamWeights(c_o_mu>0)/sum(BeamWeights(c_o_mu>0));
    
    elsc_b2b = zeros(numel(h),size(Pmu2mup,3),size(Pmu2mup,3));
    theta_lims = acos(mu_lims);
    % D_e = time_of_flight_DOMS_diffusion(E(10:20:end),gradE(10:20:end),theta_lims,0);
    D_e = time_of_flight_DOMS_diffusion(E,gradE,theta_lims,0); % seems
                                                               % to
                                                               % take ~5s
    D_e = D_e.*CscD_e;
    
    
    for iE = length(E):-1:1,
      
      elsc_b2b = 0*elsc_b2b;
    %   try
    %     waitbar(iE/length(E),wbh);
    %   catch
    %     % Whatever, it becomes tedious...
    %   end
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
        % TODO: In with curr_phfc_2nde and  curr_phfc_primary_e
        [curr_phfc_e,curr_phfc_i] = ang_scat_fcn(E(iE),theta(:));
        curr_phfc_e = curr_phfc_e.*sin(theta');
        curr_phfc_e = curr_phfc_e/sum(curr_phfc_e);
        B2B = beams2beams(curr_phfc_e,Pmu2mup,theta2beamW);
        curr_phfc_i = curr_phfc_i.*sin(theta');
        curr_phfc_i = curr_phfc_i/sum(curr_phfc_i);
        A2B = beams2beams(curr_phfc_i,Pmu2mup,theta2beamW);
        
        % add the elastic back scattering to the total
        tot_sc = tot_sc + n_i*xs_i(1,iE);
        for i_1 = size(B2B,1):-1:1,
          for i_2 = size(B2B,2):-1:1
            elsc_b2b(:,i_1,i_2) = elsc_b2b(:,i_1,i_2) + n_i*(xs_i(1,iE)*B2B(i_1,i_2));
          end
        end
        
        % add the inelastic collisions and ionisations to the total
        for j2 = 2:size(xs_i,1),
          tot_sc = tot_sc + n_i*(xs_i(j2,iE));
          for i_1 = size(B2B,1):-1:1,
            for i_2 = size(B2B,2):-1:1
              % The second factor corrects for the case where the energy loss
              % dE_i(j2,1) is smaller than the width in energy in the energy
              % bin.
              % That is when gradE(iE) > dE_i(j2,1) only the fraction
              % dE_i(j2,1)/gradE is lost from the energy bin [E(iE) E(iE)+gradE(iE)]
              elsc_b2b(:,i_1,i_2) = ( elsc_b2b(:,i_1,i_2) + ...
                                      n_i*(xs_i(j2,iE)*A2B(i_1,i_2)).*max(0,1-dE_i(j2,1)./gradE(iE)));
            end
          end
        end
        
      end
      
      % Add the loss due to electron-electron collisions:
      A = tot_sc + dEds_ee(E(iE),ne,Te)/gradE(iE);
      I0_beams = I0(:,min(iE,size(I0,2)));
      for i_mu = numel(c_o_mu):-1:1
        if c_o_mu < 0
          W = BeamWeights(i_mu)/sum(BeamWeights(c_o_mu<0));
        else
          W = BeamWeights(i_mu)/sum(BeamWeights(c_o_mu>0));
        end
        p_e_Q((1:numel(h))+(i_mu-1)*numel(h)) = p_e_q(:,iE)*W;
      end
      
      
      
      load Ie_incoming.mat Ie_total % load the incoming flux
      
      % for the first run
      if size(Ie_total{1},2) == 1
        for i_mu = 1:(length(mu_lims)-1)
        Ie_p{i_mu} = Ie_total{i_mu}(iE) * ones(1,length(t)); %(#e/s/m^2)
        end
      % for the following runs
      else
        for i_mu = 1:(length(mu_lims)-1)
          Ie_p{i_mu} = Ie_total{i_mu}(iE,(1:51)+(n_loop-1)*50);
        end
      end
      
      
      
      DE = gradE([iE,min(numel(gradE),iE+1)]);
      % This is the call to the standard Crank-Nicholson
      % PDE-integrating function, with  central differences for the
      % spatial gradient
      %    [Ie_zt] = de_M_stream_CNzt(h/cos(mag_ze),t,c_o_mu,...
      %                               Ie_p,...
      %                               I0_beams,...
      %                               v_of_E(E(iE)),...
      %                               deg_Ie(:,:,iE)+repmat(p_e_Q(:),size(t)),...
      %                               A,elsc_b2b,...
      %                               ne.*L(:,iE),ne.*dLdE(:,iE),Ie_ztEpdE,DE,...
      %                               0);
      % This is the call to the modified Crank-Nicolson PDE-integrating
      % function, with up-stream differences for the spatial gradients
      % D_e = min(1e-3,DaFcn(E(iE)));
      
      [Ie_zt] = de_M_stream_CNztusD_CFL64(h/cos(mag_ze),t,c_o_mu,...
                                    Ie_p,...
                                    I0_beams,...
                                    v_of_E(E(iE)),...
                                    deg_Ie(:,:,iE)+repmat(p_e_Q(:),size(t)),...
                                    A,elsc_b2b,...
                                    D_e(iE,:));
      % TODO: Check if this is often applied outside of at
      % too-low-altitudes:
      Ie_zt(Ie_zt(:)<0) = 0;
      
      Ie_ztE(:,:,iE) = Ie_zt;
      Ie_ztEpdE = Ie_zt;
      
      %% Calculation of the energy degradation and production of
      %  secondary electrons, starting off with adding the energy
      %  degradation due to electron-electron collisions:
      if iE > 1
        % This was sure an error, but will it still be?
        deg_Ie(:,:,iE-1) = deg_Ie(:,:,iE-1) + ... % dEds_ee(E(iE),ne,Te)/gradE(iE).*Ie_ztE(:,:,iE).*gradient(h);
            repmat(dEds_ee(E(iE),ne,Te)/gradE(iE),size(c_o_mu(:),1),size(t,2)).*Ie_ztE(:,:,iE);
        %repmat(dEds_ee(E(iE),ne,Te)/gradE(iE).*gradient(h),size(c_o_mu(:),1),size(t,2)).*Ie_ztE(:,:,iE);
      end 
      % then continuing done species by species
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
        % for idE = 2:length(dE_i)-1,
        for idE = 2:size(dE_i,1),
          if dE_i(idE,2) > 0
            break;
            % then we're off to ionisation's
          end
            
          % The flux of electrons degraded from energy bin [E(iE) E(iE)+gradE(iE)]
          Ie_degraded = (xs_i(idE,iE).*min(1,dE_i(idE,1)./gradE(iE)))*A_b2b*Ie_zt;
          % to any lower energy bin by excitation of the idE-th state
          % of the current species.
          
          % Find the energy bins that the electrons in the current
          % energy bin will degrade to when losing dE_i(idE,1) eV
          i_upper = find( E + gradE > E(iE)-dE_i(idE,1) & ...
                          E < E(iE)+gradE(iE)-dE_i(idE,1) );
          partition_fract = zeros(size(i_upper));
          if ( ~isempty(i_upper) && i_upper(1) < iE ) %Â¤
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
            % Calculate the flux of degrading electrons, sum of those
            % back-scattered, and those that scatter in the forward
            % direction:
            for i_u = find(partition_fract~=0),%1:numel(i_upper),
              deg_Ie(:,:,i_upper(i_u)) = ( deg_Ie(:,:,i_upper(i_u)) + ...
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
              Ionization = 0*Ie_zt;
              Ionizing = repmat(n_i,numel(c_o_mu),numel(t)).*(xs_i(idE2,iE) * Ie_zt);
              for i_mu = numel(c_o_mu):-1:1,
                Ionization((1:nZ)+nZ*(i_mu-1),:) = max(0,...
                                                       repmat(n_i,1,numel(t)).*(xs_i(idE2,iE) * Ie_zt((1:nZ)+nZ*(i_mu-1),:))*...
                                                       BeamWeights(i_mu)/sum(BeamWeights));
              end
            else
              Ionizing   = max(0,...
                               repmat(n_i,numel(c_o_mu),numel(t)).*(xs_i(idE2,iE) * Ie_zt));
              Ionization = max(0,...
                               repmat(n_i,numel(c_o_mu),numel(t)).*(xs_i(idE2,iE) * Ie_zt));
            end
            % Calculate the spectra of the secondary electrons, using the
            % average energy of the electrons in the current energy bin:
            second_e_spectra = feval(second_e_fcn,E,E(iE),dE_i(idE2,1),'s',AURORA_root_directory);
            second_e_spectra = (second_e_spectra + second_e_spectra([2:end end]))/2.*gradE;
            % And the distribution of the ionising electrons (that have
            % to lose the corresponding amount of energy)
            % deg_sp = abs(feval(second_e_fcn,max(1,min(( E(iE) - dE_i(idE2,1) )/2,max(0,E_p_d))),E(iE),3));
            deg_sp = second_e_fcn(E,E(iE),dE_i(idE2,1),'c',AURORA_root_directory);
            
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
            
            % and add this to the degrading electron fluxes, here there
            % should be "no" angular-scattering:
            % This part should be possible to move out of the loop over
            % the species
            % TODO: insert the corresponding angular scattering for
            % ionizing primary electrons here!/BG-20180921
            for iI = 1:iE-1,
              deg_Ie(:,:,iI) = ( deg_Ie(:,:,iI) + ...
                                 Ionization*second_e_spectra(iI) +... % Isotropic secondaries
                                 Ionizing*deg_sp(iI) ); % Forward-scattering primary-e
              if e_ionized_distribution(iI) == 0
                break
              end
            end
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
%       waitbar((length(E)-iE)/length(E),wbh,sprintf('%d',iE));
    end
% close(wbh)