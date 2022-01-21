% Entferne Daten zu PKM, die bereits aus der Datenbank gelöscht wurden

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

parroblib_gen_bitarrays();
repopath=fileparts(which('parroblib_path_init.m'));
% Gehe Datenbank durch und finde zunächst alle Daten-Ordner
EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
 for j = 1:size(EEFG_Ges,1)
   PNames_Kin = parroblib_filter_robots(EEFG_Ges(j,1:6));
   EEstr = sprintf('%dT%dR', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
   robdirs = dir(fullfile(repopath, ['sym_', EEstr], 'P*'));
   fprintf('Untersuche Roboter mit EE-FG %s. %d Ordner.\n', EEstr, length(robdirs));
   for i = 1:length(robdirs)
     if ~robdirs(i).isdir, continue; end
     RobName = robdirs(i).name;
     if ~any(contains(PNames_Kin, [RobName, 'G']))
       fprintf(['Ordner %s gehört zu keinem Roboter in der csv/mat-Daten', ...
         'bank. Lösche den Ordner.\n'], RobName);
       rmdir(fullfile(robdirs(i).folder, RobName), 's');
     end
   end
 end