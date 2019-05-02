% Finde das gesuchte Robotermodell in der Datenbank (csv-Tabelle)
% 
% Eingabe:
% NLEG [1x1]
%   Anzahl der Beine
% Leg_Names [1xNLEG cell-Array]
%   (falls nur ein Wert: Symmetrisch)
% Actuation [1xNLEG cell-Array]
%   Nummern der aktuierten Gelenke jeder Beinkette.
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
% (C) Institut für mechatronische Systeme, Universität Hannover

function [found, Name, PName_Kin, csvline_kin, csvline_act] = parroblib_find_robot(NLEG, LEG_Names, Actuation, symrob)

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
PName_Kin = '';
kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
fid = fopen(kintabfile);
if fid == -1
  warning('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile);
  return
end
tline = fgetl(fid);
while ischar(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';');
  tline = fgetl(fid); % nächste Zeile
  if isempty(csvline) || strcmp(csvline{1}, '')
    continue
  end
  if length(csvline) ~= 8
    % warning('Zeile %s sieht ungültig aus', tline);
    continue % nicht genug Spalten: Ungültiger Datensatz oder Überschrift
  end
  % Vergleich mit Eingabedaten
  if ~strcmp(csvline{2}, LEG_Names{1})
    continue
  end
  % Bis hier hin gekommen: Roboter-Kinematik wurde gefunden
  PName_Kin = csvline{1};
  csvline_kin = csvline;
  found(1) = true;
  break;
end
fclose(fid);

% Anzahl der Beingelenke aus Namen der seriellen Beinkette herausfinden
% (funktioniert nur bei seriellen Ketten. Andernfalls muss man die Abfrage
% über die Datenbank machen)
LegJointDOF = str2double(LEG_Names{1}(2)); % Format SxRRPR...
%% Roboteraktuierung in Tabelle suchen
acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, 'actuation.csv');
fid = fopen(acttabfile);
if fid == -1
  % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
  return
end
tline = fgetl(fid);
while ischar(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';');
  tline = fgetl(fid); % nächste Zeile
  if isempty(csvline) || strcmp(csvline{1}, '')
    continue
  end
  if length(csvline) ~= (NLEG*LegJointDOF+1)
    % warning('Zeile %s sieht ungültig aus', tline);
    continue % nicht genug Spalten: Ungültiger Datensatz oder Überschrift
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
