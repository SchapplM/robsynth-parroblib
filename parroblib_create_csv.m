% csv-Datei für alle vorhanden PKM aus Datenbank erzeugen
% 
% Eingabe:
% NLEG: Anzahl der Beinketten
%   3: 3T0R
%   4: 3T1R

% Junnan Li, ljn931022@gmail.com, 2020-05
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_create_csv(NLEG)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
EE_FG_Mask = logical([1 1 1 1 1 1]); % hat keine Wirkung
if NLEG == 3
  csvtable = fullfile(repopath, sprintf('synthesis_result_lists'), sprintf('3T0R.csv'));
  EE_FG0 = logical([1 1 1 0 0 0]);
elseif NLEG == 4
  csvtable = fullfile(repopath, sprintf('synthesis_result_lists'), sprintf('3T1R.csv'));
  EE_FG0 = logical([1 1 1 0 0 1]);
end
if ~exist(csvtable, 'file')
  % Tabelle existiert nicht. Ordner erstellen.
  mkdirs(fileparts(csvtable));
end

%% Kinematik-Tabelle durchsuchen
[PNames_Kin, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots( ...
  NLEG, EE_FG0, EE_FG_Mask, 6);

%% PKM von Datenbank zu csv eingeben
tableheading = {'Beinkette','Anzahl_Beinketten','G_Methode','P_Methode','Status', 'Anzahl_Akt'};
All_PKM = {};
for i = 1:length(PNames_Kin)
  % Finde alle Aktuierungen zu dieser Kinematik
  AnzahlErfolgAkt_i = 0;
  for j = 1:length(PNames_Akt)
    if strcmp( PNames_Kin{i}, PNames_Akt{j}(1:length(PNames_Kin{i})) )
      if AdditionalInfo_Akt(j) == 0 % kein Rangverlust
        AnzahlErfolgAkt_i = AnzahlErfolgAkt_i + 1;
      end
    end
  end
  Status = 0; % Da die PKM in der Datenbank ist, konnte sie offensichtlich vorher geprüft werden.
  [~,LEG_Names,~, Coupling] = parroblib_load_robot(sprintf('%s%s', PNames_Kin{i}, 'A1'));
  Robcell = {LEG_Names{1},NLEG,Coupling(1),Coupling(2),Status,AnzahlErfolgAkt_i}; % PKM Zeile 
  All_PKM = [All_PKM;Robcell]; %#ok<AGROW>
end
All_PKM_sort = sortrows(All_PKM,[1 3 4]);
tabledata = cell2table(All_PKM_sort);
tabledata.Properties.VariableNames = tableheading;
writetable(tabledata,csvtable, 'Delimiter', ';');
fprintf('Datei %s angelegt.\n', csvtable);
