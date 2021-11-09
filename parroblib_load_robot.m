% Daten für mit Namen gegebenes Robotermodell laden.
% Abgrenzung zu parroblib_find_robot.m: Dort werden die Eigenschaften des
% Roboters gegeben und der Name zurückgegeben. Hier umgekehrt
% 
% Eingabe:
% Name
%   Name des Roboters in der Datenbank (Kinematik und Aktuierung)
%   Format: P6RRPRRR14V3G1P1A1
% Modus
%   Schalter für Umfang des Zugriffs auf die Datenbank
%   0: Sehr schnell. Keine Dateien öffnen (ermöglicht nur wenige
%      Informationen)
%   1: Alle Informationen. Längerer Lesezugriff
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
% StructuralDHParam
%   Information, welcher der freien strukturbestimmenden DH-Parameter
%   (alpha und theta) auf welchen Wert gelegt werden darf. Variable Anzahl
%   von Einträgen, je nachdem, welche Kombinationen möglich sind.
%   Jeder Parameter kann 4 Zustände haben (1-4).
%   Die Stelle der Zahl gibt die Reihenfolge des Parameters an. Kodierung:
%      0: "u" - Nicht definiert
%      1: "p" - Nur 0°
%      2: "o" - nur 90°
%      3: "b" - 0° oder 90°
%      4: "a" - Alle Werte möglich.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [NLEG, LEG_Names, Actuation, Coupling, ActNr, symrob, EE_dof0, ...
  PName_Kin, PName_Legs, AdditionalInfo_Akt, StructuralDHParam] = ...
  parroblib_load_robot(Name, Modus)
%% Initialisierung
if nargin < 2
  Modus = 1;
end
assert(isa(Name, 'char'), 'Name muss char sein');
NLEG = 0;
LEG_Names = {};
Actuation = {};
Coupling = [0 0];
ActNr = 0;
symrob = true;
EE_dof0 = NaN(1,6);
PName_Kin = '';
AdditionalInfo_Akt = NaN(1,1);
StructuralDHParam = {};
repopath=fileparts(which('parroblib_path_init.m'));

% Name der Kinematischen Struktur von Aktuierung trennen
% Namensschema für symmetrische, serielle PKM als regulären Ausdruck
% TODO: Anpassung an nicht-serielle, nicht-symmetrische PKM
expression = 'P(\d)([RP]+)(\d+)[V]?(\d*)[G]?(\d*)[P]?(\d*)A(\d+)'; % Format "P3RRR1G1P1A1" oder "P3RRR1V1G1P1A1"
[tokens, ~] = regexp(Name,expression,'tokens','match');
expression_kin = 'P(\d)([RP]+)(\d+)[V]?(\d*)';
[tokens_kin, ~] = regexp(Name,expression_kin,'tokens','match');
if isempty(tokens_kin)
  error('Eingegebener Name %s entspricht nicht dem Namensschema', Name);
elseif isempty(tokens)
  % Eingabe entspricht dem Beinketten-Namen ohne G-/P-Nummer und Aktuierung
  Name = [Name, 'G1P1A0']; % Fülle mit Platzhaltern auf
  [tokens, ~] = regexp(Name,expression,'tokens','match');
end
res_kin = tokens_kin{1};
NLEG = str2double(res_kin{1});
if isempty(res_kin{4}) % serielle Kette ist keine abgeleitete Variante
  PName_Legs = ['P', res_kin{1}, res_kin{2}, res_kin{3}];
else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
  PName_Legs = ['P', res_kin{1}, res_kin{2}, res_kin{3}, 'V', res_kin{4}];
end
LEG_Names = repmat({['S', sprintf('%d',length(res_kin{2})), PName_Legs(3:end)]}, 1, NLEG);
if isempty(tokens)
  error('Eingabe ist ungültig (Fall darf nicht auftreten)');
end
res = tokens{1};
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
if Modus == 0
  % Keine Tabellen öffnen. Nur Plausibilitätsprüfung und Extraktion von
  % Informationen aus dem Roboternamen.
  return
end
%% csv-Tabelle öffnen: Kinematik
% Da Anzahl Beinketten gegeben ist, aber die PKM nach FG gespeichert sind,
% müssen mehrere Tabellen durchsucht werden.
EEFG_Ges = logical(...
  [1 1 0 0 0 0; 1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
EEstr = ''; % Platzhalter, wird im folgenden belegt.
for jj = 1:size(EEFG_Ges,1)
  if sum(EEFG_Ges(jj,:)) ~= NLEG, continue; end % PKM-FG passen nicht zu Beinketten
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(jj,1:3)), sum(EEFG_Ges(jj,4:6)));
  % Ergebnis: Tabellenzeile csvline_kin für den gesuchten Roboter
  kintabfile = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list.csv']);
  fid = fopen(kintabfile);
  if fid == -1
    error('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile);
  end
  % Tabelle zeilenweise durchgehen
  tline = fgetl(fid);
  found = false;
  while ischar(tline) && ~isempty(tline)
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
  if found
    break; % In aktueller Tabelle gefunden. Abbruch.
  end
end
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
acttabfile = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
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
  while ischar(tline) && ~isempty(tline)
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
  % Die einzelnen Möglichkeiten sind Komma-getrennte Strings, bei der
  % jede Stelle einen Parameter kodiert. In den meisten Fällen dürften
  % nicht viele Kombinationen möglich sein.
  StructuralDHParam = strsplit(csvline_act{k+1}, ',', 'CollapseDelimiters', false);
  if length(csvline_act) ~= k+1
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
