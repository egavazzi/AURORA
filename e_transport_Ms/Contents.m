% E_TRANSPORT_MS, 
% (ver 20180705)
% 
% Multi-stream electron transport equations, steady-state and
% time-dependent electron-transport functions
% 
% Time-dependent electron transport functions
%  Ie_Mstream_tz_2_aurora - time-dependent multi-stream electron transport
%  de_M_stream_CNzt       - multistream Crank-Nicholson central difference electron transport
%  de_M_stream_CNztus     - multistream Crank-Nicholson upstream difference electron transport
%
% Steady-state electron transport functions
%  Ie_M_stream_4_aurora  - steady-state multi-stream electron transport
%  Ie_M_stream_2_photo_e - 2-hemispheres steady-state multi-stream electron transport
%  de_M_stream_CNztus    - multistream electron transport solver
%  de_M_stream2HS_us     - 2-hemispheres multistream electron transport solver
%
% Multi-stream functions functions for pitch-angle scattering
%  e_scattering_beamdistribution - angular redistribution PDF for multi-stream electron transport
%  beams2beams                   - multi-beam scattering redistribution coefficients 
%  mu_avg                        - pitch-angle averages between limits of cosine-of-pitch-angles
%
% Correction for finite energy-pitch-angle cells:
%  time_of_flight_DOMS_diffusion - Energy-pitch-angle time-of-flight widening
%  de_M_stream_CNztusD           - multistream Crank-Nicholson electron
%                                  transport with diffusion for correcting for
%                                  finite energy-pitch-angle cells:

%   Copyright © 2019-2020 Bjorn Gustavsson, <bjorn.gustavsson@irf.se>
%   This is free software, licensed under GNU GPL version 2 or later
