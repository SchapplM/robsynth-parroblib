% Füge ein Robotermodell zur Datenbank (csv-Tabelle) hinzu
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
% sym_xTyR/sym_xTyR_list.csv (Liste aller Roboterstrukturen)
% actuation.csv (Tabelle der Aktuierungen)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function [Name, new] = parroblib_add_robot(NLEG, LEG_Names, Actuation, Coupling, EEdof0)
%% Init
repopath=fileparts(which('parroblib_path_init.m'));
if length(Actuation) ~= NLEG
  error('Zu jeder Beinkette muss die Aktuierung gegeben sein');
end
assert(isa(Coupling, 'double') && all(size(Coupling) == [1 2]), ...
  'Koppelpunkt-Kennung muss 1x2 double sein');
EEstr = sprintf('%dT%dR', sum(EEdof0(1:3)), sum(EEdof0(4:6)));
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
[found, Name, PName_Kin] = parroblib_find_robot(NLEG, LEG_Names, Actuation, Coupling, symrob);
if found(2)
  new = false;
  return
end

% Es wird auf jeden Fall ein neuer Eintrag erstellt.
new = true;
%% Erstelle neuen Eintrag: Kinematik-Tabelle
if ~found(1)
  % Für diese Kinematik gibt es noch keinen Eintrag
  kintabfile = fullfile(repopath, ['sym_%s', EEstr], ['sym_',EEstr,'_list.csv']);
  kintabtmp1file = fullfile(repopath, ['sym_%s', EEstr], ['sym_',EEstr,'_list.tmp1.csv']);
  kintabtmp2file = fullfile(repopath, ['sym_%s', EEstr], ['sym_',EEstr,'_list.tmp2.csv']);
  % Kopfzeilendatei erstellen (Temp-Datei zum Kopieren)
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
  fid = fopen(kintabtmp1file, 'w');
  fwrite(fid, [line_head1, newline]);
  fwrite(fid, [line_head2, newline]);
  fclose(fid);
  % Namen der Kinematik-Struktur generieren
  % Namensschema PxRRPRyyGuuPvv
  PName_Kin = sprintf('P%d%sG%dP%d',NLEG,LEG_Names{1}(3:end),Coupling(1),Coupling(2));
  
  % Roboterstruktur an Tabelle anhängen
  csvline_robkin = {PName_Kin, LEG_Names{1}};
  % Spalten für EE-Freiheitsgrade
  c=2;
  for i = 1:6
    c = c+1; 
    if ~isempty(EEdof0)
      csvline_robkin{c} = EEdof0(i);
    else
      csvline_robkin{c} = 0; % Setze FG auf Null. Wenn alle 0 sind, ist der unbekannte Status gekennzeichnet.
    end
  end
  % Tabelle in Matlab öffnen
  T = readtable(kintabfile, 'NumHeaderLines', 0, 'ReadVariableNames', 0, 'Delimiter', ';');
  % Zeile ans Ende der temporären Datentabelle einfügen
  newrow = cell2table(csvline_robkin);
  newrow.Properties.VariableNames = T.Properties.VariableNames;
  T = [T;newrow];
  % Daten alphabetisch sortieren
  T_sort = sortrows(T,1);
  % Tabelle temporär schreiben
  writetable(T_sort, kintabtmp2file, 'WriteVariableNames', 0, 'Delimiter', ';');
  % Beide Dateien (Kopfzeilen und Daten) kombinieren
  fid = fopen(kintabtmp1file, 'a') ;
  fwrite(fid, fileread(kintabtmp2file)) ;
  fclose(fid);
  % Kopieren der temporären Dateien und löschen
  copyfile(kintabtmp1file, kintabfile);
  delete(kintabtmp1file);
  delete(kintabtmp2file);
else
  % Roboter ist in Kinematik-Tabelle, aber noch nicht in
  % Aktuierungs-Tabelle
end

%% Erstelle neuen Eintrag: Aktuierungs-Tabelle
% Anzahl der Bein-FG aus Namen der Beinkette feststellen. TODO: Das ist auf
% symmetrische serielle PKM begrenzt
LegJointDOF = str2double(LEG_Names{1}(2)); % Format SxRRPR...
% Extrahiere Beinketten-Name der PKM (ohne Ausrichtungs-Nummern G/P)
expression = '(P[\d][RP]+[\d]+[V]?[\d]*)G[\d]+P[\d]+'; % Format "P3RRR1G1P1A1" oder "P3RRR1V1G1P1A1"
[tokens, ~] = regexp(PName_Kin,expression,'tokens','match');
PName_Legs = tokens{1}{1};

acttabfile = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
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
  c=c+1; csvline_head1{c} = 'Rangverlust Plattform-FG';
  % Überschrift für Wertebereich der Winkel. Kann mit Skript
  % "modify_csv_column_theta1.m" noch aktualisiert werden
  % (Enthält dann die Variablennamen auf die sich die Spalte bezieht).
  c=c+1; csvline_head1{c} = 'Wertebereich freie Winkel';
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
while ischar(tline) && ~isempty(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';');
  tline = fgetl(fid); % nächste Zeile
  if isempty(csvline) || strcmp(csvline{1}, '')
    continue
  end
  expression = [PName_Kin, 'A(\d+)']; % Format "P3RRR1G1P1A1"
  [tokens, ~] = regexp(csvline{1},expression,'tokens','match');
  if isempty(tokens)
    continue
  end
  % Aktuelle laufende Nummer aus der Tabelle holen
  A_lfdNr = str2double(tokens{1}{1});
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
c=c+1; csv_act_rob{c} = '?'; % Keine Information über Rangverlust vorliegend
c=c+1; csv_act_rob{c} = ''; % Keine Information über beschränkten Wertebereich vorliegend
% Cell-Array in csv-Zeile umwandeln
line_act = csv_act_rob{1};
for i = 2:length(csv_act_rob)
  line_act = sprintf('%s;%s', line_act, csv_act_rob{i});
end
% Füge neuen Eintrag hinzu (in alphabetischer Reihenfolge)
acttabfile_copy = [acttabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung
fid = fopen(acttabfile, 'r');
fidc = fopen(acttabfile_copy, 'w');
tline = fgetl(fid);
Name_old = '';
written = false;
i=0;
while ischar(tline) && ~isempty(tline)
  i=i+1;
  csvline = strsplit(tline, ';', 'CollapseDelimiters', false); % Spaltenweise als Cell-Array
  if strcmp(csvline{1}, ''), continue; end
  Name_now = csvline{1};
  if i > 1 && string(Name) > string(Name_old) && string(Name) < string(Name_now) && ~written
    fwrite(fidc, [line_act, newline]);
    written = true;
  end
  fwrite(fidc, [tline, newline]);
  Name_old = Name_now;
  tline = fgetl(fid); % nächste Zeile der ursprünglichen Tabelle
  if ~ischar(tline) && ~written
    % Erster Roboter: Nur die Überschriften wurden geschrieben. Schreibe
    fwrite(fidc, [line_act, newline]);
  end
end
fclose(fid);
fclose(fidc);
% Modifizierte Tabelle zurückkopieren
copyfile(acttabfile_copy, acttabfile);
% Kopie-Tabelle löschen
delete(acttabfile_copy);

