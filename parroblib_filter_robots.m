% Finde PKM-Robotermodelle in der Datenbank mit gegebenen Filtern
% 
% Eingabe:
% NLEG
%   Anzahl der Beine der symmetrischen PKM
% EE_FG [1x6]
%   Vektor mit 1/0 für Belegung ob EE-FG aktiv ist
% EE_FG_Mask [1x6]
%   Maske, die festlegt ob die FG exakt wie in `EE_FG` sind, oder ob auch
%   gesperrte FG wirklich nicht beweglich sind
% 
% Rückgabe:
% PNames_Kin
%   Cell-Array mit Namen aller Roboterkinematiken ohne Hinweis auf
%   Aktuierung
% PNames_Akt
%   Cell-Array aller Roboter mit Aktuierung als Namensbestandteil.
%   Diese Namen können durch die Funktion parroblib_load_robot abgerufen
%   werden
% 
% TODO: Aktuell sind nur symmetrische PKM berücksichtigt.
% 
% Siehe auch: serroblib_filter_robots

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-01
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [PNames_Kin, PNames_Akt] = parroblib_filter_robots(NLEG, EE_FG0, EE_FG_Mask)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
%% Kinematik-Tabelle durchsuchen
PNames_Kin = {};
PNames_Akt = {};
kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
fid = fopen(kintabfile);
if fid == -1
  warning('Datei existiert nicht. Sollte aber im Repo enthalten sein.');
  return
end
% Tabelle zeilenweise durchgehen
tline = fgetl(fid);
while ischar(tline)
  % Zeile spaltenweise als Cell-Array
  csvline = strsplit(tline, ';');
  tline = fgetl(fid); % nächste Zeile
  if isempty(csvline) || strcmp(csvline{1}, '')
    continue
  end
  if length(csvline) ~= 8
    % warning('Zeile %s sieht ungültig aus', tline);
    continue % nicht genug Spalten: Ungültiger Datensatz oder Überschrift
  end
  
  % Gültige Tabellenzeile. Daten auslesen
  PName_Kin = csvline{1};
  LegName = csvline{2}; %#ok<NASGU>
  EEFG0_cell = csvline(3:8);
  EE_FG_csv = zeros(1,6);
  for i = 1:6
    if strcmp(EEFG0_cell{i}, '1')
      EE_FG_csv(i) = 1;
    end
  end
  
  % Daten filtern
  if ~all( (EE_FG_csv & EE_FG_Mask) == (EE_FG0 & EE_FG_Mask))
    continue
  end
  
  % Roboter passt zu Filter. Auswahl in Liste
  PNames_Kin = {PNames_Kin{:}, PName_Kin}; %#ok<CCAT>
end
fclose(fid);

%% Aktuierung heraussuchen
PNames_Akt = {};
for i = 1:length(PNames_Kin)
  PName_Kin = PNames_Kin{i};
  acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, 'actuation.csv');
  fid = fopen(acttabfile);
  if fid == -1
    % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
    return
  end
  tline = fgetl(fid);
  firstline = true;
  while ischar(tline)
    % Zeile spaltenweise als Cell-Array
    csvline = strsplit(tline, ';');
    tline = fgetl(fid); % nächste Zeile
    if firstline % erste Zeile ist Überschrift. Ignorieren
      firstline = false;
      continue % TODO: Robustere Methode zur Erkennung richtiger Zeilen
    end
    if isempty(csvline) || strcmp(csvline{1}, '')
      continue
    end
    % Daten auslesen
    Name_Akt = csvline{1};
    % Roboter in Liste aufnehmen (es gibt keinen Filter für die Aktuierung)
    PNames_Akt = {PNames_Akt{:}, Name_Akt}; %#ok<CCAT>
  end
  fclose(fid);
end