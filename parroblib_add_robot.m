% Füge ein Robotermodell zur Datenbank (csv-Tabelle) hinzu
% 
% Eingabe:
% NLEG [1x1]
%   Anzahl der Beine
% Leg_Names [1xNLEG cell-Array]
%   (falls nur ein Wert: Symmetrisch)
% Actuation [1xNLEG cell-Array]
%   Nummern der aktuierten Gelenke jeder Beinkette.
% EEdof0
%   Vektor mit beweglichen EE-FG des Roboters (Geschw. und Winkelgeschw. im
%   Basis-KS. Entspricht Vorgabe in der Struktursynthese von Ramirez)
%   1="Komponente durch Roboter beeinflussbar"; 0="nicht beeinflussbar"
% 
% Ausgabe:
% Name
%   Name der PKM (Kinematikstuktur und Aktuierung) in der Datenbank
% new
%   true, wenn Roboter neu ist und der Datenbank hinzugefügt wurde.
% 
% Schreibt Dateien:
% sym3leg/sym3leg_list.csv (Liste aller Roboterstrukturen)
% actuation.csv (Tabelle der Aktuierungen)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function [Name, new] = parroblib_add_robot(NLEG, LEG_Names, Actuation, EEdof0)
%% Init
repopath=fileparts(which('parroblib_path_init.m'));
% Prüfe, ob kinematisch symmetrisch (durch Vergleich der Namen der
% Beinketten)
symrob = true;
if length(LEG_Names) == 1
  symrob = true;
  name = LEG_Names{1};
  LEG_Names = cell(NLEG,1);
  for i = 1:NLEG, LEG_Names{i} = name; end
else
  Leg1Name = LEG_Names{1};
  for i = 2:NLEG
    if ~strcmp(Leg1Name, LEG_Names{i})
      symrob = false;
      break;
    end
  end
end

% Suche in Datenbank nach dem Roboter
[found, Name, PName_Kin] = parroblib_find_robot(NLEG, LEG_Names, Actuation, symrob);
if found(2)
  new = false;
  return
end

% Es wird auf jeden Fall ein neuer Eintrag erstellt.
new = true;
%% Erstelle neuen Eintrag: Kinematik-Tabelle
% TODO: Aktuell wird die Nummerierung nicht angepasst. Für jede serielle
% Beinkette gibt es nur eine zugehörige PKM, egal wie die Ausrichtung der
% Gestell-Koppelgelenke ist.
if ~found(1)
  % Für diese Kinematik gibt es noch keinen Eintrag
  kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
  if ~exist(kintabfile, 'file')
    % Tabelle existiert nicht. Erstellen
    mkdirs(fileparts(kintabfile));
    csvline_head1 = {'Name','Beinkette','EE-FG (Basis-KS)','','','','',''};
    csvline_head2 = {'','','vx0','vy0','vz0','wx0','wy0','wz0'};
    % String aus Cell-Array erzeugen
    line_head1 = csvline_head1{1};
    line_head2 = csvline_head2{1};
    for i = 2:length(csvline_head1)
      line_head1 = sprintf('%s;%s', line_head1, csvline_head1{i});
      line_head2 = sprintf('%s;%s', line_head2, csvline_head2{i});
    end
    % Kopfzeile in csv-Tabelle schreiben
    fid = fopen(kintabfile, 'a');
    fwrite(fid, [line_head1, newline]);
    fwrite(fid, [line_head2, newline]);
    fclose(fid);
  end
  % Namen der Kinematik-Struktur generieren
  PName_Kin = ['P', num2str(NLEG), LEG_Names{1}(3:end)]; % Namensschema PxRRPRyy
  
  % Roboterstruktur an Tabelle anhängen
  csvline_robkin = {PName_Kin, LEG_Names{1}};
  % Spalten für EE-Freiheitsgrade
  c=2;
  for i = 1:6
    c = c+1; 
    if ~isempty(EEdof0)
      csvline_robkin{c} = sprintf('%d',EEdof0(i));
    else
      csvline_robkin{c} = '';
    end
  end
  % Cell-Array in csv-Zeile umwandeln
  line_robot = csvline_robkin{1};
  for i = 2:length(csvline_robkin)
    line_robot = sprintf('%s;%s', line_robot, csvline_robkin{i});
  end
  % Zeile in Datei anhängen
  fid = fopen(kintabfile, 'a');
  fwrite(fid, [line_robot, newline]);
  fclose(fid);
else
  % Roboter ist in Kinematik-Tabelle, aber noch nicht in
  % Aktuierungs-Tabelle
end

%% Erstelle neuen Eintrag: Aktuierungs-Tabelle
% Anzahl der Bein-FG aus Namen der Beinkette feststellen. TODO: Das ist auf
% symmetrische serielle PKM begrenzt
LegJointDOF = str2double(LEG_Names{1}(2)); % Format SxRRPR...

acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, 'actuation.csv');
if ~exist(acttabfile, 'file')
  % Tabelle existiert nicht. Erstellen
  % Kopfzeile als Cell-Array
  csvline_head1 = {'Name'};
  % Für jedes Gelenk der Beinkette angeben, ob dieses Aktuiert ist oder
  % nicht
  % Anzahl der Beingelenke aus Namen der seriellen Beinkette herausfinden
  c=1;
  for iL = 1:NLEG
    for j = 1:LegJointDOF
      c=c+1;
      if j == 1
        csvline_head1{c} = sprintf('Aktuierung Bein %d', iL);
      else
        % Überschrift geht über mehrere Spalten
        csvline_head1{c} = '';
      end
    end
  end
  % String aus Cell-Array erzeugen
  line_head1 = csvline_head1{1};
  for i = 2:length(csvline_head1)
    line_head1 = sprintf('%s;%s', line_head1, csvline_head1{i});
  end
  % Kopfzeile in csv-Tabelle schreiben
  mkdirs(fileparts(acttabfile));
  fid = fopen(acttabfile, 'w');
  fwrite(fid, [line_head1, newline]);
  fclose(fid);
end

% Suche letzten Eintrag in der Tabelle
A_lfdNr = 0;
fid = fopen(acttabfile);
tline = fgetl(fid);
while ischar(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';');
  tline = fgetl(fid); % nächste Zeile
  if isempty(csvline) || strcmp(csvline{1}, '')
    continue
  end
  expression = [PName_Kin, 'A(\d+]']; % Format "P3RRR1A1"
  [tokens, ~] = regexp(csvline{1},expression,'tokens','match');
  if isempty(tokens)
    continue
  end
  % Aktuelle laufende Nummer aus der Tabelle holen
  A_lfdNr = str2double(tokens{1}{3});
end
fclose(fid);
% neuer Eintrag hat eine höhere Nummer als der vorherige
A_lfdNr_neu = A_lfdNr + 1;

% Neuen Namen erstellen
Name = sprintf('%sA%d',PName_Kin,A_lfdNr_neu);

% Zeile zusammenstellen
csv_act_rob = {Name};
c=1;
for iL = 1:NLEG
  % Binär-Indizes für Aktuierung, wie in die CSV-Tabelle geschrieben werden
  % soll
  ActSel = false(LegJointDOF,1);
  ActSel(Actuation{iL}) = true;
  for j = 1:LegJointDOF
    c=c+1;
    csv_act_rob{c} = sprintf('%d', ActSel(j));
  end
end
% Cell-Array in csv-Zeile umwandeln
line_act = csv_act_rob{1};
for i = 2:length(csv_act_rob)
  line_act = sprintf('%s;%s', line_act, csv_act_rob{i});
end
% Füge neuen Eintrag hinzu
fid = fopen(acttabfile, 'a');
fwrite(fid, [line_act, newline]);
fclose(fid);

