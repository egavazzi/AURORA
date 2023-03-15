% EXTRACTS THE DISTRIBUTION FUNCTION OF IONOSPHERIC
% ELECTRONS AT THE IONOSPHERIC BOUNDARY (BOUNDARY CONDITION).
speciesnumber = 3;


%% -------------------------------------------------------------------------- %%
ccc=pwd;
[a0,a1]=system('uname -n'); [a0,a2]=system('ps |grep -i matlab');
dummy=[a1 a2 datestr(now)];
dummy=dummy(double(dummy)~=10); dummy=dummy(double(dummy)~=32);

if ~exist('absolutelynomessages')
  absolutelynomessages=logical(0);
end

% ------- input data ------- %
% The easiest way to get the general parameters
inputb6;

% If the fields_per_file parameter wasn't defined in the input file,
% default to 1 field per file.
if ~exist('fields_per_file')
  fields_per_file = 1;
end

% Now read the species specific data
particle=struct();
fid=fopen('inputb6.m');
for ii=1:Nspecies
  theline=fgetl(fid);
  if length(theline)<5,
    theline=[theline '     '];
  end
  while ~(strcmp(theline(1:5),'%SPEC') | strcmp(theline(1:5),'%spec'))
    theline=fgetl(fid);
    if length(theline)<5,
      theline=[theline '     '];
    end
  end
  while ~(strcmp(theline(1:4),'%END') | strcmp(theline(1:4),'%end'))
    eval(theline)
    theline=fgetl(fid);
    if length(theline)<4,
      theline=[theline '    '];
    end
  end
  particle(ii).Nvz=Nvz;
  particle(ii).vzmin=vzmin;
  particle(ii).vzmax=vzmax;
  particle(ii).Nmu=Nmu;
  particle(ii).mumin=mumin;
  particle(ii).mumax=mumax;
  particle(ii).muexp=muexp;
  particle(ii).mass=mass;
  particle(ii).charge=charge;
  particle(ii).n0=n0;
  particle(ii).vz0=vz0;
  particle(ii).kTz=kTz;
  particle(ii).kTp=kTp;
  particle(ii).n0L=n0L;
  particle(ii).vz0L=vz0L;
  particle(ii).kTzL=kTzL;
  particle(ii).kTpL=kTpL;
  particle(ii).n0R=n0R;
  particle(ii).vz0R=vz0R;
  particle(ii).kTzR=kTzR;
  particle(ii).kTpR=kTpR;
end
fclose(fid);

for ii=1:Nspecies
  particle(ii).dvz=(particle(ii).vzmax-particle(ii).vzmin)/particle(ii).Nvz;
  particle(ii).vzcorn=particle(ii).vz0 + ...
      particle(ii).vzmin + particle(ii).dvz*[0:particle(ii).Nvz];
  particle(ii).vz = ...
      0.5*(particle(ii).vzcorn(1:end-1)+particle(ii).vzcorn(2:end));

  vmu=[1:particle(ii).Nmu];
  particle(ii).mu = particle(ii).mumin + ...
      0.5*((vmu.^particle(ii).muexp+(vmu-1).^particle(ii).muexp) / ...
           particle(ii).Nmu^particle(ii).muexp) * ...
      (particle(ii).mumax-particle(ii).mumin);
  particle(ii).dmu = ((vmu.^particle(ii).muexp - ...
                       (vmu-1).^particle(ii).muexp) / ...
                      particle(ii).Nmu^particle(ii).muexp) * ...
      (particle(ii).mumax-particle(ii).mumin);
  particle(ii).mucorn = particle(ii).mu-0.5*particle(ii).dmu; 
  particle(ii).mucorn = [particle(ii).mucorn particle(ii).mu(end) + ...
                      0.5*particle(ii).dmu(end)];
end


% --- distribution function fBC(z,vz,mu) --- %
% one file per timestep, containing a structured array of all species,
% with a three-dimensional array for f(z,vz,mu) in each.
cd outp

% Prevent two processes from performing simultaneous conversions
if exist([pwd 'lock.fBC'])
  if ~absolutelynomessages
    disp('Another process is already working on this directory.')
    disp('If this is not the case, remove the file')
    disp([pwd 'lock.fBC'])
  end
else
  all_is_fine=logical(0);
  try
    dlmwrite('lock.fBC',dummy,'')
    fid=fopen('lock.fBC','r');
    lock=textscan(fid,'%s');
    fclose(fid);
    if strcmp(dummy,lock{1})
      all_is_fine=logical(1);
    end
  catch
    all_is_fine=logical(0);
  end

  if all_is_fine
    try
      dd=dir('datfiles/fBC/');
      fBC=struct();
      
      for speciesnumber=Nspecies:-1:1
        if ~isnan(particle(speciesnumber).mass) & ...
              ~isinf(particle(speciesnumber).mass)
          break
        end
      end
      for ii=3:length(dd)
        if length(dd(ii).name)>=10
          if strcmp(dd(ii).name(1:2),'fR')
            if strcmp(dd(ii).name(11:12),num2str(speciesnumber,'%0.2d')) & ...
                  strcmp(dd(ii).name(13:end),'.ketchup.dat')
              jj = ii;
            end
          end
        end
      end
      
      for speciesnumber = 1:Nspecies
        if ~isnan(particle(speciesnumber).mass) & ...
                        ~isinf(particle(speciesnumber).mass)

          fBC(speciesnumber).f = ...
                zeros(particle(speciesnumber).Nvz,particle(speciesnumber).Nmu);
          infile = ['datfiles/fBC/' dd(jj).name(1:10) ...
                          num2str(speciesnumber,'%0.2d') ...
                          '.ketchup.dat'];
          fid=fopen(infile,'r');
          if fid<0
            error(['Error reading file ' infile])
          end

          for kk = 1
            instruct = textscan(fid,'%f');
            fBC(speciesnumber).f(:,:) = ...
                        reshape(instruct{1}, [particle(speciesnumber).Nmu ...
                        particle(speciesnumber).Nvz]).';
            fclose(fid);
          end
        end
      end
      
%       outfile = ['fBC_right_s0',num2str(speciesnumber),'.mat'];
      outfile = 'fBC_right.mat';
      save(outfile,'-v7.3','fBC')

      if ~absolutelynomessages
%         disp(['ionospheric BC for species ',num2str(speciesnumber),' done!'])
        disp(['Ionospheric BC done!'])
      end
      clear fBC
      delete('./lock.fBC')
    catch
      felfelfel=lasterror;
      disp(felfelfel.message)
      delete('./lock.fBC')
    end % end try
  end
end

cd(ccc)