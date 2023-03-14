% Entferne verwaiste Einträge aus der Datenbank
% Betrifft PKM aus der Kinematik-Liste, ohne Aktuierung und Aktuierungen
% ohne Kinematik
% Prüfe die CSV-Tabellen, falls parroblib_gen_bitarrays nicht mehr
% funktioniert aufgrund inkonsistenter Daten
% 
% Dieses Skript muss zweimal ausgeführt werden zur Aktualisierung der
% Mat-Dateien

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2023-03
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
dryrun = false;
repopath=fileparts(which('parroblib_path_init.m'));

EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);

for j = 1:size(EEFG_Ges,1)
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  fprintf('Prüfe Konsistenz der Einträge mit %s FG\n', EEstr);
  kintabfile_csv = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list.csv']);
  KinTab = readtable(kintabfile_csv, 'Delimiter', ';');
  PNames_Kin = KinTab.Name;
  fprintf('Aktualisiere Datenbank\n');
  parroblib_gen_bitarrays(EEFG_Ges(j,1:6));
  [~, PNames_Act] = parroblib_filter_robots(EEFG_Ges(j,1:6), 6);
  %% Prüfe verwaiste Kinematiken
  fprintf('Prüfe %d Kinematiken\n', length(PNames_Kin));
  for i = 1:length(PNames_Kin)
    % Lade die Aktuierungen direkt aus der CSV-Tabelle
    [~, ~, ~, ~, ~, ~, ~, ~, PName_Legs_i] = parroblib_load_robot(PNames_Kin{i}, 0);
    acttabfile_csv = fullfile(repopath, ['sym_', EEstr], PName_Legs_i, 'actuation.csv');
    remove = false;
    if ~exist(acttabfile_csv, 'file')
      remove = true;
    else
      ActTab_i = readtable(acttabfile_csv, 'Delimiter', ';', 'PreserveVariableNames', ...
      true, 'ReadVariableNames', true);
      if ~any(contains(ActTab_i.Name, [PNames_Kin{i}, 'A']))
        remove = true;
      end
    end
    if remove
      fprintf('Zu Kinematik %d (%s) gibt es keine Aktuierung. Lösche.\n', ...
        i, PNames_Kin{i});
      if dryrun, continue; end
      success = parroblib_remove_robot(PNames_Kin{i});
      if ~success
        error('Fehler beim Löschen von %s', PNames_Kin{i});
      end
    end
  end
  %% Prüfe verwaiste Aktuierungen
  fprintf('Prüfe %d Aktuierungen\n', length(PNames_Act));
  for i = 1:length(PNames_Act)
    [~, ~, ~, ~, ~, ~, ~, PName_Kin_i] = parroblib_load_robot(PNames_Act{i}, 0);
    if ~any(contains(PNames_Kin, PName_Kin_i))
      fprintf('Zu Aktuierung %d (%s) gibt es keine Kinematik. Lösche.\n', ...
        i, PNames_Act{i});
      if dryrun, continue; end
      success = parroblib_remove_robot(PNames_Act{i});
      if ~success
        error('Fehler beim Löschen von %s', PNames_Act{i});
      end
    end
  end
  fprintf('Aktualisiere Datenbank\n');
  parroblib_gen_bitarrays(EEFG_Ges(j,1:6));
end
