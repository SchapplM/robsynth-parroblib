% csv-Datei für alle vorhanden PKM aus Datenbank erneuen
% 
% Eingabe:
% NLEG: Anzahl der Beinketten
%   3: 3T0R
%   4: 3T1R
% Status: Status der PKM in csv-datei zu speichern
%   0=in Ordnung, Rang konnte geprüft werden
%   1=Ausschluss wegen Anzahl der Schubgelenke
%   2=Ausschluss wegen parrob_structsynth_check_leg_dof
%   3=Keine Lösbarkeit der Einzelpunkt-IK in Maßsynthese
%   4=keine Lösbarkeit der Traj.-IK in Maßsynthese
%   5=noch nicht geprueft
% Rangerfolg: 
%   0: Rangverlust
%   1: Rangerfolgreich

% Junnan Li, ljn931022@gmail.com, 2020-05
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-07
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_update_csv(SName, Coupling, NLEG, Status, Rangerfolg)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
if nargin == 4
  Rangerfolg = 0;
end

if NLEG == 3
  csvtable = fullfile(repopath, sprintf('synthesis_result_lists'), sprintf('3T0R.csv'));
  EE_FG0 = logical([1 1 1 0 0 0]);
elseif NLEG == 4
  csvtable = fullfile(repopath, sprintf('synthesis_result_lists'), sprintf('3T1R.csv'));
  EE_FG0 = logical([1 1 1 0 0 1]);
end
if ~exist(csvtable, 'file')
  warning('csv-Datei existiert nicht. Erstelle.');
  % Tabelle existiert nicht. Erstellen
  parroblib_create_csv(NLEG);
end

%% PKM in csv-Datei eingeben
T = readtable(csvtable, 'ReadVariableNames', true);
for i = 1:size(T,1)
  LegName = table2cell(T(i,1));
  Coupling_i = table2array(T(i,3:4));
  Status_i = table2array(T(i,5));
  if strcmp(LegName,SName{:}) && all(Coupling_i == Coupling)
    if Status ~= 0  % 0=keine Aktuierung erfolgreich getestet
      if Status_i == Status
        % Keine Änderung an bestehendem Inhalt der Tabelle.
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

