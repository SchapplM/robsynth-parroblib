% Finde das gesuchte Robotermodell in der Datenbank (csv-Tabelle)
% 
% Eingabe:
% NLEG [1x1]
%   Anzahl der Beine
% Leg_Names [1xNLEG cell-Array]
%   (falls nur ein Wert: Symmetrisch)
% Actuation [1xNLEG cell-Array]
%   Nummern der aktuierten Gelenke jeder Beinkette.
% Coupling [1x2]
%   Nummern des Gestell- und Plattform-Koppel-Typs (entsprechend ParRob)
% symrob [1x1 logical]
%   true, wenn es eine symmetrische PKM ist
% 
% Ausgabe:
% found [1x2 double]
%  Marker für Kinematik und Aktuierung: Falls true liegen die Daten für den
%  Roboter mit den gegebenen Eigenschaften in den csv-Tabellen
% Name
%   Name der PKM (Kinematikstuktur und Aktuierung) in der Datenbank
% PName_Kin
%   Name der Kinematikstruktur der PKM (ohne Betrachtung der Aktuierung)
% csvline_kin 
%   Zeile des Roboters in der csv-Tabelle für die Kinematik (...list.csv)
% csvline_act
%   Zeile in der csv-Tabelle für Aktuierung (actuation.csv)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function [found, Name, PName_Kin, csvline_kin, csvline_act] = parroblib_find_robot(NLEG, LEG_Names, Actuation, Coupling, symrob)

%% Initialisierung
found = false(1,2);
Name = '';
csvline_kin = {};
csvline_act = {};
repopath=fileparts(which('parroblib_path_init.m'));

if ~symrob
  error('Nicht kinematisch symmetrische Roboter noch nicht implementiert');
end
%% Roboterkinematik in Tabelle suchen
% Da Anzahl Beinketten gegeben ist, aber die PKM nach FG gespeichert sind,
% müssen mehrere Tabellen durchsucht werden.
EEFG_Ges = logical(...
  [1 1 0 0 0 0; 1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
EEstr = ''; % Platzhalter, wird im folgenden belegt.
PName_Kin = '';
for jj = 1:size(EEFG_Ges,1)
  if sum(EEFG_Ges(jj,:)) ~= NLEG, continue; end % PKM-FG passen nicht zu Beinketten
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(jj,1:3)), sum(EEFG_Ges(jj,4:6)));
  kintabfile = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list.csv']);
  fid = fopen(kintabfile);
  if fid == -1
    warning('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile);
    return
  end
  tline = fgetl(fid);
  while ischar(tline) && ~isempty(tline)
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';');
    tline = fgetl(fid); % nächste Zeile
    if isempty(csvline) || strcmp(csvline{1}, '')
      continue
    end
    if length(csvline) ~= 3
      % warning('Zeile %s sieht ungültig aus', tline);
      continue % nicht genug Spalten: Ungültiger Datensatz oder Überschrift
    end
    % Vergleich mit Eingabedaten
    if ~strcmp(csvline{2}, LEG_Names{1})
      continue
    end
    if ~strcmp(csvline{1}, sprintf('P%d%sG%dP%d', NLEG, LEG_Names{1}(3:end), Coupling(1), Coupling(2)))
      continue
    end
    % Bis hier hin gekommen: Roboter-Kinematik wurde gefunden
    PName_Kin = csvline{1};
    csvline_kin = csvline;
    found(1) = true;
    break;
  end
  fclose(fid);
  if found(1), break; end % In aktueller Tabelle gefunden. Abbruch.
end

if isempty(PName_Kin)
  % Der Kinematikname fehlt bislang in der Datenbank
  return
end

% Anzahl der Beingelenke aus Namen der seriellen Beinkette herausfinden
% (funktioniert nur bei seriellen Ketten. Andernfalls muss man die Abfrage
% über die Datenbank machen)
LegJointDOF = str2double(LEG_Names{1}(2)); % Format SxRRPR...
%% Roboteraktuierung in Tabelle suchen
% Bestimme Beinketten-Name durch Entfernen der G- und P-Angabe
expression = '(P[\d][RP]+[\d]+[V]?[\d]*)G[\d]+P[\d]+'; % Zielformat "P3RRR1A1" oder "P3RRR1V1A1" (ohne G1P1)
[tokens, ~] = regexp(PName_Kin,expression,'tokens','match');
if isempty(tokens)
  error('Name %s aus Tabelle passt nicht zum Namensschema', PName_Kin);
end
PName_Legs = tokens{1}{1};
acttabfile = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
fid = fopen(acttabfile);
if fid == -1
  % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
  return
end
tline = fgetl(fid);
while ischar(tline) && ~isempty(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
  tline = fgetl(fid); % nächste Zeile
  if strcmp(csvline{1}, '')
    continue
  end
  if length(csvline) ~= (1+NLEG*LegJointDOF+2)
    warning('Zeile %s... sieht ungültig aus (Länge: %d)', csvline{1}, length(csvline));
    continue % nicht genug Spalten: Ungültiger Datensatz oder Überschrift
  end
  if ~contains(csvline(1), [PName_Kin, 'A'])
    continue % Kinematik in Aktuierungs-Tabelle ist andere G-/P-Nummer
  end
  % Vergleich der Tabellenzeile mit den Eingabedaten
  k=2;
  abort=false;
  for iL = 1:NLEG
    csv_act_vect_leg = NaN( LegJointDOF, 1 );
    for j = 1:LegJointDOF
      csv_act_vect_leg(j) = str2double(csvline{k});
      k=k+1;
    end
    ActSel = false(LegJointDOF,1);
    ActSel(Actuation{iL}) = true;
    if any(csv_act_vect_leg(ActSel)~=1) || any(csv_act_vect_leg(~ActSel)~=0)
      % Die Aktuierung stimmt nicht
      abort=true;
      break;
    end
  end
  if abort
    continue; % die Zeile enthält nicht den gesuchten Roboter. Probiere den nächsten.
  end
  % Bis hierhin gekommen: Die gewünschte Aktuierung ist schon in der
  % Tabelle
  Name = csvline{1};
  found(2) = true;
  csvline_act = csvline;
  break;
end
fclose(fid);
