% CSV-Datei für alle vorhanden PKM aus Datenbank erneuen
% Erzeugt/Aktualisiert Datei synthesis_result_lists/aTbR.csv ({a,b}={0,1,2,3})
% 
% Eingabe:
% SName:
%   Name der seriellen Führungskette (symmetrische PKM)
% Coupling [1x2]:
%   Nummern der Gestell- und Plattform-Koppelpunkt-Varianten
% EE_FG0:
%   Endeffektor-FG (Geschwindigkeitskomponenten)
% Status: Status der PKM in csv-datei zu speichern
%   0=in Ordnung, Rang konnte geprüft werden
%   1=Ausschluss wegen Anzahl der Schubgelenke
%   2=Ausschluss wegen parrob_structsynth_check_leg_dof
%   3=Keine Lösbarkeit der Einzelpunkt-IK in Maßsynthese
%   4=keine Lösbarkeit der Traj.-IK in Maßsynthese
%   5=Parasitäre Bewegung
%   6=noch nicht geprüft
%   7=Unbehandelter Fall (Fehler)
% Rangerfolg: 
%   0: Rangverlust in Jacobi-Matrix. PKM verliert FG und funktioniert nicht.
%   1: Kein Rangverlust in Jacobi-Matrix. Gültige PKM.
% 
% Für erstmalige Generierung der Datei aus den Dateien actuation.csv:
% parroblib_update_csv('', [], logical([1 1 1 1 1 1]))

% Junnan Li, ljn931022@gmail.com, 2020-05
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_update_csv(SName, Coupling, EE_FG0, Status, Rangerfolg)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
if nargin <= 4
  Rangerfolg = 0;
end
assert(isa(EE_FG0, 'logical'), 'Eingabe EE_FG0 muss 1x6 logical sein')

if all(EE_FG0 == logical([1 1 0 0 0 1]))
  csvtable = fullfile(repopath, 'synthesis_result_lists', sprintf('2T1R.csv'));
elseif all(EE_FG0 == logical([1 1 1 0 0 0]))
  csvtable = fullfile(repopath, 'synthesis_result_lists', sprintf('3T0R.csv'));
elseif all(EE_FG0 == logical([1 1 1 0 0 1]))
  csvtable = fullfile(repopath, 'synthesis_result_lists', sprintf('3T1R.csv'));
elseif all(EE_FG0 == logical([1 1 1 1 1 1]))
  csvtable = fullfile(repopath, 'synthesis_result_lists', sprintf('3T3R.csv'));
else
  error('Fall noch nicht implementiert');
end
if ~exist(csvtable, 'file')
  warning('csv-Datei existiert nicht. Erstelle: %s', csvtable);
  % Tabelle existiert nicht. Erstellen
  parroblib_create_csv(EE_FG0, csvtable);
end
if isempty(SName)
  fprintf('Keine Beinkette gegeben. Aufruf nur zum Erstellen der CSV. Abbruch.\n');
  return
end
NLEG = sum(EE_FG0);
%% PKM in csv-Datei eingeben
T = readtable(csvtable, 'ReadVariableNames', true);
% Passenden Index der Eingabedaten in der Tabelle herausfinden
I_name = strcmp(table2cell(T(:,1)), SName);
I_coupl = table2array(T(:,3))==Coupling(1) & ...
          table2array(T(:,4))==Coupling(2);
i = find(I_name&I_coupl);
if length(i) > 1
  error('Inkonsistenz in Tabelle %s. Eintrag doppelt: %s, %d, %d', csvtable, ...
    SName, Coupling(1), Coupling(2));
elseif ~isempty(i)
  LegName = table2cell(T(i,1));
  Coupling_i = table2array(T(i,3:4));
  Status_i = table2array(T(i,5));
  if ~strcmp(LegName,SName{:}) || ~all(Coupling_i == Coupling)
    error('Inkonsistenz beim Laden der Datei %s', csvtable);
  end
  if Status ~= 0  % 0=keine Aktuierung erfolgreich getestet
    if Status_i == Status
      % Keine Änderung an bestehendem Inhalt der Tabelle.
      return
    elseif Status == 6
      % Wurde als "nicht geprüft" aufgerufen. Da schon ein Eintrag
      % vorliegt, wurde die PKM aber bereits mit einem anderem Status
      % geprüft, der mehr Information beinhält.
      % Möglich, wenn die Struktursynthese erneut durchgeführt wird (mit
      % Option "dryrun"). Daher hier Abbruch.
      return
    end
    % Status aktualisieren.
    T(i,5) = {Status}; % Status eintragen
    T(i,6) = {0}; % Null erfolgreiche Aktuierungen per Definition
    writetable(T,csvtable,'Delimiter', ';');
  else % PKM Erfolgreich (schon in Datenbank) 
    SName_str = cell2str(SName); % Name der Beinkette (aus Gesamttabelle)
    % Name der PKM-Kinematik zusammenstellen
    PName = sprintf('P%d%sG%dP%d', NLEG, SName_str(3:end),Coupling(1),Coupling(2));
    % Finde alle Roboter in der Datenbank mit dieser Anzahl Beinen
    [~, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots( ...
      NLEG, EE_FG0, [1 1 1 1 1 1], 6);
    % Finde alle Aktuierungen zu dieser Kinematik
    AnzahlErfolgAkt_i = 0;
    for j = 1:length(PNames_Akt)
      if strcmp( PName, PNames_Akt{j}(1:length(PName)) )
        if AdditionalInfo_Akt(j) == 0 % kein Rangverlust
          AnzahlErfolgAkt_i = AnzahlErfolgAkt_i + 1;
        end
      end
    end
    T(i,5) = {num2str(0)}; 
    T(i,6) = {num2str(AnzahlErfolgAkt_i)};
    writetable(T,csvtable,'Delimiter', ';');
  end
  return % hinzuzufügende PKM ist schon in Datenbank. Abbruch.
end

% Bis hierhin gekommen. Also ist die PKM noch nicht drin gewesen.
%% PKM in Tabelle einfügen
% Neue Tabellenzeile als Tabelle definieren
NewRow = cell2table({SName, NLEG, Coupling(1),Coupling(2),Status,Rangerfolg});
% Überschriften für Zeilen-Tabelle setzen, damit Tabellen kombinierbar werden
NewRow.Properties.VariableNames = T.Properties.VariableNames;
T_new = [T;NewRow];
% Tabelle sortieren und speichern
T_sort = sortrows(T_new,[1 3 4]);
writetable(T_sort,csvtable, 'Delimiter', ';');

end

% csv-Datei für alle vorhanden PKM aus Datenbank erzeugen
% 
% Eingabe:
% EE_FG0: Endeffektor-FG (Geschwindigkeitskomponenten)

% Junnan Li, ljn931022@gmail.com, 2020-05
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_create_csv(EE_FG0, csvtable)

%% Kinematik-Tabelle durchsuchen
EE_FG_Mask = logical([1 1 1 1 1 1]); % hat keine Wirkung
NLEG = sum(EE_FG0);
[PNames_Kin, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots( ...
  NLEG, EE_FG0, EE_FG_Mask, 6);

%% PKM von Datenbank zu csv eingeben
tableheading = {'Beinkette','Anzahl_Beinketten','G_Methode','P_Methode','Status', 'Anzahl_Akt'};
All_PKM = {};
for i = 1:length(PNames_Kin)
  % Finde alle Aktuierungen zu dieser Kinematik
  AnzahlErfolgAkt_i = 0;
  for j = 1:length(PNames_Akt)
    if length(PNames_Kin{i}) > length(PNames_Akt{j})
      % Namen passen nicht zueinander. j-Aktuierung kann nicht zu
      % i-Kinematik passen.
      continue
    end
    if strcmp( PNames_Kin{i}, PNames_Akt{j}(1:length(PNames_Kin{i})) )
      % j-Aktuierung passt zu i-Kinematik
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
if ~exist(csvtable, 'file')
  % Tabelle existiert nicht. Ordner erstellen.
  mkdirs(fileparts(csvtable));
end
writetable(tabledata, csvtable, 'Delimiter', ';');
fprintf('Datei %s angelegt.\n', csvtable);
end
