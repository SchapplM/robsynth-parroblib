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
% max_rankdeficit [1x1]
%   Grad des Rangverlustes der PKM-Jacobi-Matrix (Standard: 0)
%   Wird ein Wert von 1 bis 6 (max.) eingesetzt, werden auch PKM
%   ausgegeben, deren gewählte Aktuierung keine Mobilität ergibt.
% 
% Rückgabe:
% PNames_Kin
%   Cell-Array mit Namen aller Roboterkinematiken ohne Hinweis auf
%   Aktuierung
% PNames_Akt
%   Cell-Array aller Roboter mit Aktuierung als Namensbestandteil.
%   Diese Namen können durch die Funktion parroblib_load_robot abgerufen
%   werden
% AdditionalInfo_Akt
%   Array mit zusätzlichen Infos für alle Strukturen aus PNames_Akt (in den Zeilen).
%   Spalten:
%   1: Rangverlust der Jacobi-Matrix (in den vorgesehenen FG der PKM)
% 
% TODO: Aktuell sind nur symmetrische PKM berücksichtigt.
% 
% Siehe auch: serroblib_filter_robots

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-01
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [PNames_Kin, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots(NLEG, EE_FG0, EE_FG_Mask, max_rankdeficit)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
if nargin < 4
  max_rankdeficit = 0;
end
%% Kinematik-Tabelle durchsuchen
PNames_Kin = {};
PNames_Akt = {};
AdditionalInfo_Akt = [];
kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
fid = fopen(kintabfile);
if fid == -1
  warning('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile);
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
  % Extrahiere Beinketten-Name der PKM (ohne Ausrichtungs-Nummern G/P)
  % Unter diesem Namen ist die Aktuierungstabelle gespeichert
  expression = 'P(\d)([RP]+)(\d+)[V]?(\d*)[G]?(\d*)[P]?(\d*)'; % Format "P3RRRG1P11A1" oder "P3RRR1V1G1P1A1"
  [tokens, ~] = regexp(PName_Kin,expression,'tokens','match');
  if isempty(tokens)
    error('Eintrag %s aus der Tabelle stimmt nicht mit Namensschema überein.', PName_Kin);
  end
  res = tokens{1};
  if isempty(res{4}) % serielle Kette ist keine abgeleitete Variante
    PName_Legs = ['P', res{1}, res{2}, res{3}];
  else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
    PName_Legs = ['P', res{1}, res{2}, res{3}, 'V', res{4}];
  end
  acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Legs, 'actuation.csv');
  fid = fopen(acttabfile);
  if fid == -1
    % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
    warning('Zu Kinematik %s gibt es keine Aktuierungstabelle %s', PName_Kin, acttabfile);
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
    
    % Prüfe, ob Kinematik die gleiche ist (es stehen mehrere
    % Gestell-Varianten in der gleichen Aktuierungstabelle)
    if ~contains(Name_Akt, PName_Kin)
      continue
    end
    % Wert für den Rangverlust (die Datenbank enthält auch PKM, deren
    % Aktuierung oder Gestell-Anordnung nicht sinnvoll ist)
    rankdef = str2double(csvline{end});
    if rankdef > max_rankdeficit
      % Die Struktur soll nicht ausgegeben werden: Rangverlust ist zu groß
      % Wenn in Tabelle ein "?" steht, dann NaN. Soll aber ausgegeben werden
      % (Bedeutet, dass noch nicht geprüft ist)
      continue
    end
    AdditionalInfo_Akt = [AdditionalInfo_Akt; rankdef]; %#ok<AGROW>
    
    % Roboter in Liste aufnehmen (es gibt keinen Filter für die Aktuierung)
    PNames_Akt = {PNames_Akt{:}, Name_Akt}; %#ok<CCAT>
  end
  fclose(fid);
end