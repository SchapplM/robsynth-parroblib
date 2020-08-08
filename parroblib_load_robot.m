% Daten für mit Namen gegebenes Robotermodell laden.
% Abgrenzung zu parroblib_find_robot.m: Dort werden die Eigenschaften des
% Roboters gegeben und der Name zurückgegeben. Hier umgekehrt
% 
% Eingabe:
% Name
%   Name des Roboters in der Datenbank (Kinematik und Aktuierung)
%   Format: P6RRPRRR14V3G1P1A1
% 
% Ausgabe:
% NLEG [1x1]
%   Anzahl der Beine
% Leg_Names [1xNLEG cell-Array]
%   (falls nur ein Wert: Symmetrisch)
% Actuation [1xNLEG cell-Array]
%   Nummern der aktuierten Gelenke jeder Beinkette.
% Coupling [1x2]
%   Nummern des Gestell- und Plattform-Koppel-Typs (entsprechend ParRob)
% ActNr
%   Laufende Nummer dieser Aktuierungsmöglichkeit der Roboterkinematik
% symrob [1x1 logical]
%   true, wenn es eine symmetrische PKM ist
% EEdof0 [1x6]
%   Vektor mit beweglichen EE-FG des Roboters (Geschw. und Winkelgeschw. im
%   Basis-KS. Entspricht Vorgabe in der Struktursynthese von Ramirez)
%   1="Komponente durch Roboter beeinflussbar"; 0="nicht beeinflussbar"
% PName_Kin
%   Eindeutiger Name der PKM-Kinematik (ohne Aktuierung)
% PName_Legs
%   Eindeutiger Name der Führungsketten-Zusammensetzung (ohne Aktuierung
%   und Gestell- oder Plattformausrichtung)
% AdditionalInfo_Akt
%   Zusätzliche Infos. Spalten:
%   1: Rangverlust der Jacobi-Matrix (in den vorgesehenen FG der PKM)
%   2: Möglicher Wertebereich für den ersten frei wählbaren Strukturwinkel
%      (z.B. für theta1)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [NLEG, LEG_Names, Actuation, Coupling, ActNr, symrob, EE_dof0, ...
  PName_Kin, PName_Legs, AdditionalInfo_Akt] = ...
  parroblib_load_robot(Name)
%% Initialisierung
NLEG = 0;
LEG_Names = {};
Actuation = {};
ActNr = 0;
symrob = true;
EE_dof0 = NaN(1,6);
AdditionalInfo_Akt = NaN(1,2);
repopath=fileparts(which('parroblib_path_init.m'));

% Name der Kinematischen Struktur von Aktuierung trennen
% Namensschema für symmetrische, serielle PKM als regulären Ausdruck
% TODO: Anpassung an nicht-serielle, nicht-symmetrische PKM
expression = 'P(\d)([RP]+)(\d+)[V]?(\d*)[G]?(\d*)[P]?(\d*)A(\d+)'; % Format "P3RRR1G1P1A1" oder "P3RRR1V1G1P1A1"
[tokens, ~] = regexp(Name,expression,'tokens','match');
if isempty(tokens)
  error('Eingegebener Name %s entspricht nicht dem Namensschema', Name);
end
res = tokens{1};
NLEG = str2double(res{1});
if isempty(res{4}) % serielle Kette ist keine abgeleitete Variante
  PName_Legs = ['P', res{1}, res{2}, res{3}];
else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
  PName_Legs = ['P', res{1}, res{2}, res{3}, 'V', res{4}];
end
if isempty(res{5})
  % Für Kompatibilität zu alten Aufrufen der Funktion. Akzeptiere auch
  % Eingaben der Form
  warning('Eingegebener Robotername %s entspricht altem Format ohne G-/P-Nummer', Name);
  Coupling = [1 1];
else
  Coupling = [str2double(res{5}), str2double(res{6})];
end
PName_Kin = [PName_Legs, sprintf('G%dP%d', Coupling(1), Coupling(2))];

ActNr = str2double(res{7});
PName_Akt = [PName_Kin, sprintf('A%d', ActNr)];
%% csv-Tabelle öffnen: KinematikEE_dof0 = NaN(6,1);
% Ergebnis: Tabellenzeile csvline_kin für den gesuchten Roboter
kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
fid = fopen(kintabfile);
if fid == -1
  error('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile);
end
% Tabelle zeilenweise durchgehen
tline = fgetl(fid);
found = false;
while ischar(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';');
  tline = fgetl(fid); % nächste Zeile
  if isempty(csvline) || strcmp(csvline{1}, '')
    continue
  end
  if strcmp(csvline{1}, PName_Kin)
    % gefunden
    csvline_kin = csvline;
    found = true;
    break;
  end
end
fclose(fid);

if ~found
  error('Roboter %s wurde nicht in der Tabelle %s gefunden.', PName_Kin, kintabfile);
end

%% Ausgabe für Kinematik speichern
LEG_Names = csvline_kin(2);
symrob = true; % TOOD: Berücksichtigung nicht-kinematisch-symmetrischer PKM
% identische Beinketten nochmal in die Variable schreiben
for i = 1:NLEG
  LEG_Names{i} = LEG_Names{1};
end
%% csv-Tabelle öffnen: Aktuierung
% Ergebnis: Tabellenzeile csvline_act für den gesuchten Roboter
acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Legs, 'actuation.csv');
csvline_act = [];
if ActNr ~= 0
  fid = fopen(acttabfile);
  if fid == -1
    warning('Aktuierungstabelle %s existiert nicht', acttabfile);
    % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
    return
  end
  % Tabelle zeilenweise durchgehen
  tline = fgetl(fid);
  while ischar(tline)
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
    tline = fgetl(fid); % nächste Zeile
    if isempty(csvline) || strcmp(csvline{1}, '')
      continue
    end
    if strcmp(csvline{1}, PName_Akt)
      % gefunden
      csvline_act = csvline;
      break;
    end
  end
  fclose(fid);
end
%% Aktuierung abspeichern
Actuation = cell(NLEG,1);
if ~isempty(csvline_act)
  % Aktuierung in Zahl-Index-Format umwandeln. Siehe parroblib_find_robot.m
  LegJointDOF = str2double(LEG_Names{1}(2)); % Format SxRRPR... (nehme zweites Zeichen)
  k=2;
  for iL = 1:NLEG
    ActSel = NaN(LegJointDOF,1);
    for j = 1:LegJointDOF
      ActSel(j) = str2double(csvline_act{k});
      k=k+1;
    end
    Actuation{iL} = find(ActSel);
  end
  % Zusätzliche Informationen kommen am Ende der Zeile
  % Rangverlust der Jacobi-Matrix
  AdditionalInfo_Akt(1) = str2double(csvline_act{k});
  % Mögliche Zahlenwerte für frei einstellbare Winkel wie theta1.
  AdditionalInfo_Akt(2) = 0; % nicht definiert
  if length(csvline_act) == k+1
    % Die Werte sind gesetzt (ist vorerst optional).
    if strcmp(csvline_act{k+1}, '0')
      AdditionalInfo_Akt(2) = 1;
    elseif strcmp(csvline_act{k+1}, '90')
      AdditionalInfo_Akt(2) = 2;
    elseif strcmp(csvline_act{k+1}, '090')
      AdditionalInfo_Akt(2) = 3;
    elseif strcmp(csvline_act{k+1}, '*') || strcmp(csvline_act{k+1}, '')
      AdditionalInfo_Akt(2) = 4;
    else
      error('Wert "%s" ist nicht definiert für Wertebereich eines freien Winkels', csvline_act{k+1});
    end
  else % TODO: Ermögliche Wertebereiche für mehrere Winkel vorzugeben
    error('Länge %d der Zeile ist nicht definiert', length(csvline_act));
  end
elseif ActNr ~= 0
  warning('Aktuierung %d (%s) nicht in Aktuierungstabelle %s gefunden', ActNr, PName_Kin, acttabfile);
end

%% EE-FG abspeichern
% EE-FG sind in Kinematik-Tabelle enthalten
for i = 1:6
  EE_dof0(i) = str2double(csvline_kin{2+i});
end
