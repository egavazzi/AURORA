dir0 = '/data/to_DATA/BGU001-test/SteadyStateRuns201904';
cd(dir0)
dDir = dir('*Mono*')

for i1 = 1:numel(dDir),
  cd(dir0)
  cd(dDir(i1).name)
  load IeSteadyState-01.mat Ie_zE curr_BW E0dE f_f
  Q_N2p = exc_tz_of_Ie_zE(h_atm,E(1:size(Ie_zE,2)),Ie_zE,nN2,xsN2ion(1:size(Ie_zE,2)));
  Q_O2p = exc_tz_of_Ie_zE(h_atm,E(1:size(Ie_zE,2)),Ie_zE,nO2,xsO2ion(1:size(Ie_zE,2)));
  Q_Op = exc_tz_of_Ie_zE(h_atm,E(1:size(Ie_zE,2)),Ie_zE,nO,xsOion(1:size(Ie_zE,2)));
  Q_all(i1).ionN2 = Q_N2p;
  Q_all(i1).ionO2 = Q_O2p;
  Q_all(i1).ionO = Q_Op;
  Q_all(i1).dOmega = curr_BW;
  Q_all(i1).dirname = dDir(i1).name;
end
idx4_3keV = [19 25 31 37 43 46 49];

char(Q_all(idx4_3keV).dirname)
figure
idx4_3keVFA = [19 25 31 37 43 49 52 55];
for i1 = 1:numel(idx4_3keVFA),
  ph3keV(i1) = semilogx([Q_all(idx4_3keVFA(i1)).ionN2]./max(Q_all(idx4_3keVFA(i1)).ionN2),h_atm/1e3);
  hold on
end
set(ph3keV(end-2:end-1),'linewidth',2)
set(ph3keV(4),'linewidth',2)          

figure
for i1 = 1:numel(idx4_3keVFA),
  ph3keV(i1) = semilogx([Q_all(idx4_3keVFA(i1)).ionN2],h_atm/1e3);
  hold on
end
set(ph3keV(end-2:end-1),'linewidth',2)
set(ph3keV(4),'linewidth',2)          

legend(ph3keV,char(Q_all(idx4_3keVFA).dirname))
figure
idx4_3keVISO = [idx4_3keVFA(1:5)+1,58 61 64];
for i1 = 1:numel(idx4_3keVISO),
  ph3keV(i1) = semilogx([Q_all(idx4_3keVISO(i1)).ionN2]./max(Q_all(idx4_3keVISO(i1)).ionN2),h_atm/1e3);
  hold on
end

axis([1e-6 2 80 300])
legend(ph3keV,char(Q_all(idx4_3keVISO).dirname))

idx4_5keVFA = [21    27    33    39    45    50    53    56];
idx4_5keVISO = [idx4_5keVFA(1:5)+1,59 62 65];
idx4_7keVFA = [idx4_5keVISO(1:5)+1,51 54 57];
idx4_7keVISO = [idx4_7keVFA(1:5)+1,60 63 66];

idx4_3keVFA = [19 25 31 37 43 49 55 58 61];
idx4_3keVISO = [20 26 32 42 44 50 64 67 70];
idx4_5keVFA = cell2mat(femfa(:,2)');
idx4_5keVISO = cell2mat(femiso(:,2)');
idx4_7keVFA = cell2mat(sjufa(:,2)');
idx4_7keVISO = cell2mat(sjuiso(:,2)');

idxis4plot = {idx4_3keVFA,idx4_3keVISO,idx4_5keVFA,idx4_5keVISO,idx4_7keVFA,idx4_7keVISO};

for i0 = 1:numel(idxis4plot),
  figure
  i2p = idxis4plot{i0};
  for i1 = 1:numel(i2p),
    ph3keV(i1) = semilogx([Q_all(i2p(i1)).ionN2]./max(Q_all(i2p(i1)).ionN2),h_atm/1e3);
    hold on
  end
  axis([1e-6 2 80 300])
  legend(ph3keV,char(Q_all(i2p).dirname))
  set(ph3keV(end-2:end-1),'linewidth',2)
  set(ph3keV(4),'linewidth',2)          
end
DN = {'MonoSSFA20_1-3000-3keV-0-3',1;
      'MonoSSFA20_1-3000-5keV-0-3',2;
      'MonoSSFA20_1-3000-7keV-0-3',3;
      'MonoSSFA20_2-3000-3keV-0-3',4; 
      'MonoSSFA20_2-3000-5keV-0-3',5;
      'MonoSSFA20_2-3000-7keV-0-3',6; 
      'MonoSSFA9-3000-3keV-0-3',7; 
      'MonoSSFA9-3000-5keV-0-3',8; 
      'MonoSSFA9-3000-7keV-0-3',9; 
      'MonoSSISO20_1-3000-3keV-0-3',10; 
      'MonoSSISO20_1-3000-5keV-0-3',11; 
      'MonoSSISO20_1-3000-7keV-0-3',12; 
      'MonoSSISO20_2-3000-3keV-0-3',13; 
      'MonoSSISO20_2-3000-5keV-0-3',14; 
      'MonoSSISO20_2-3000-7keV-0-3',15; 
      'MonoSSISO9-3000-3keV-0-3',16; 
      'MonoSSISO9-3000-5keV-0-3',17; 
      'MonoSSISO9-3000-7keV-0-3',18; 
      'nMonoSS10B1-3000-3keV-FA',19; 
      'nMonoSS10B1-3000-3keV-ISO',20; 
      'nMonoSS10B1-3000-5keV-FA',21; 
      'nMonoSS10B1-3000-5keV-ISO',22; 
      'nMonoSS10B1-3000-7keV-FA',23; 
      'nMonoSS10B1-3000-7keV-ISO',24; 
      'nMonoSS10B2-3000-3keV-FA',25; 
      'nMonoSS10B2-3000-3keV-ISO',26; 
      'nMonoSS10B2-3000-5keV-FA',27; 
      'nMonoSS10B2-3000-5keV-ISO',28; 
      'nMonoSS10B2-3000-7keV-FA',29; 
      'nMonoSS10B2-3000-7keV-ISO',30; 
      'nMonoSS9B1530-3000-3keV-FA',31; 
      'nMonoSS9B1530-3000-3keV-ISO',32; 
      'nMonoSS9B1530-3000-5keV-FA',33; 
      'nMonoSS9B1530-3000-5keV-ISO',34; 
      'nMonoSS9B1530-3000-7keV-FA',35; 
      'nMonoSS9B1530-3000-7keV-ISO',36;
      'nMonoSS9B1530B-3000-3keV-FA',37;
      'nMonoSS9B1530B-3000-5keV-FA',38;
      'nMonoSS9B1530B-3000-5keV-ISO',39;
      'nMonoSS9B1530B-3000-7keV-FA',40;
      'nMonoSS9B1530B-3000-7keV-ISO',41;
      'nMonoSS9B51530B-3000-3keV-ISO',42;
      'nMonoSS9B54-3000-3keV-FA',43;
      'nMonoSS9B54-3000-3keV-ISO',44; 
      'nMonoSS9B54-3000-5keV-FA',45; 
      'nMonoSS9B54-3000-5keV-ISO',46; 
      'nMonoSS9B54-3000-7keV-FA',47; 
      'nMonoSS9B54-3000-7keV-ISO',48; 
      'nMonoSS9B7p5-3000-3keV-FA',49; 
      'nMonoSS9B7p5-3000-3keV-ISO',50; 
      'nMonoSS9B7p5-3000-5keV-FA',51; 
      'nMonoSS9B7p5-3000-5keV-ISO',52; 
      'nMonoSS9B7p5-3000-7keV-FA',53; 
      'nMonoSS9B7p5-3000-7keV-ISO',54; 
      'nMonoSSFA20_1-3000-3keV-0-3',55; 
      'nMonoSSFA20_1-3000-5keV-0-3',56; 
      'nMonoSSFA20_1-3000-7keV-0-3',57; 
      'nMonoSSFA20_2-3000-3keV-0-3',58; 
      'nMonoSSFA20_2-3000-5keV-0-3',59; 
      'nMonoSSFA20_2-3000-7keV-0-3',60; 
      'nMonoSSFA9-3000-3keV-0-3',61; 
      'nMonoSSFA9-3000-5keV-0-3',62; 
      'nMonoSSFA9-3000-7keV-0-3',63; 
      'nMonoSSISO20_1-3000-3keV-0-3',64; 
      'nMonoSSISO20_1-3000-5keV-0-3',65; 
      'nMonoSSISO20_1-3000-7keV-0-3',66; 
      'nMonoSSISO20_2-3000-3keV-0-3',67; 
      'nMonoSSISO20_2-3000-5keV-0-3',68; 
      'nMonoSSISO20_2-3000-7keV-0-3',69; 
      'nMonoSSISO9-3000-3keV-0-3',70; 
      'nMonoSSISO9-3000-5keV-0-3',71; 
      'nMonoSSISO9-3000-7keV-0-3',72};



trefa = {'nMonoSS10B1-3000-3keV-FA',19; 
         'nMonoSS10B2-3000-3keV-FA',25; 
         'nMonoSS9B1530-3000-3keV-FA',31; 
         'nMonoSS9B1530B-3000-3keV-FA',37;
         'nMonoSS9B54-3000-3keV-FA',43;
         'nMonoSS9B7p5-3000-3keV-FA',49; 
         'nMonoSSFA20_1-3000-3keV-0-3',55; 
         'nMonoSSFA20_2-3000-3keV-0-3',58; 
         'nMonoSSFA9-3000-3keV-0-3',61}; 


treiso = {'nMonoSS10B1-3000-3keV-ISO',20; 
          'nMonoSS10B2-3000-3keV-ISO',26; 
          'nMonoSS9B1530-3000-3keV-ISO',32; 
          'nMonoSS9B51530B-3000-3keV-ISO',42;
          'nMonoSS9B54-3000-3keV-ISO',44; 
          'nMonoSS9B7p5-3000-3keV-ISO',50; 
          'nMonoSSISO20_1-3000-3keV-0-3',64; 
          'nMonoSSISO20_2-3000-3keV-0-3',67; 
          'nMonoSSISO9-3000-3keV-0-3',70}; 

femfa = {'nMonoSS10B1-3000-5keV-FA',21; 
         'nMonoSS10B2-3000-5keV-FA',27; 
         'nMonoSS9B1530-3000-5keV-FA',33; 
         'nMonoSS9B1530B-3000-5keV-FA',38;
         'nMonoSS9B54-3000-5keV-FA',45; 
         'nMonoSS9B7p5-3000-5keV-FA',51; 
         'nMonoSSFA20_1-3000-5keV-0-3',56; 
         'nMonoSSFA20_2-3000-5keV-0-3',59; 
         'nMonoSSFA9-3000-5keV-0-3',62};

femiso = {'nMonoSS10B1-3000-5keV-ISO',22; 
          'nMonoSS10B2-3000-5keV-ISO',28; 
          'nMonoSS9B1530-3000-5keV-ISO',34; 
          'nMonoSS9B1530B-3000-5keV-ISO',39;
          'nMonoSS9B54-3000-5keV-ISO',46; 
          'nMonoSS9B7p5-3000-5keV-ISO',52; 
          'nMonoSSISO20_1-3000-5keV-0-3',65; 
          'nMonoSSISO20_2-3000-5keV-0-3',68; 
          'nMonoSSISO9-3000-5keV-0-3',71}; 

sjufa = {'nMonoSS10B1-3000-7keV-FA',23; 
         'nMonoSS10B2-3000-7keV-FA',29; 
         'nMonoSS9B1530-3000-7keV-FA',35; 
         'nMonoSS9B1530B-3000-7keV-FA',40;
         'nMonoSS9B54-3000-7keV-FA',47; 
         'nMonoSS9B7p5-3000-7keV-FA',53; 
         'nMonoSSFA20_1-3000-7keV-0-3',57; 
         'nMonoSSFA20_2-3000-7keV-0-3',60; 
         'nMonoSSFA9-3000-7keV-0-3',63};

sjuiso = {'nMonoSS10B1-3000-7keV-ISO',24; 
          'nMonoSS10B2-3000-7keV-ISO',30; 
          'nMonoSS9B1530-3000-7keV-ISO',36;
          'nMonoSS9B1530B-3000-7keV-ISO',41;
          'nMonoSS9B54-3000-7keV-ISO',48; 
          'nMonoSS9B7p5-3000-7keV-ISO',54; 
          'nMonoSSISO20_1-3000-7keV-0-3',66; 
          'nMonoSSISO20_2-3000-7keV-0-3',69; 
          'nMonoSSISO9-3000-7keV-0-3',72};



