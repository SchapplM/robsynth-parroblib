% Entferne ein Robotermodell aus der Datenbank (csv-Tabelle)
% Die Zeile des Robotermodells (z.B. P3PRP2) wird aus der csv-Tabelle der
% Gelenkstruktur (z.B. sym3leg.csv) entfernt. Zusätzlich werden die
% Matlab-Funktionen aus dem gleichnamigen Unterordner ("P3PRP2") entfernt.
% 
% Eingabe:
% PName
%   Name des Roboters in der Datenbank. Optionen:
%   * PName_Kin (nur Kinematik, ohne Endung für Aktuierung; z.B. "P3RPR1G1P1")
%   * PName_Act (Kinematik mit Aktuierung; z.B. "P3RPR1G1P1A2")
%   Je nachdem welche Art von Name übergeben wird, wird der Eintrag gelöscht
% nofiledelete (optional)
%   true: Keine Dateien (mit generiertem Code löschen). Nur umbenennen.
%   false: Ordner werden als Sicherung nur verschoben anstatt gelöscht
%   (Standard)
% 
% Ausgabe:
% success
%   true, falls Roboter erfolgreich entfernt wurde
% 
% Siehe auch: serroblib_remove_robot

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function success = parroblib_remove_robot(PName_Input, nofiledelete)

if nargin < 2
  nofiledelete = false;
end
if ~isa(PName_Input, 'char')
  error('Der Datentyp von PName_Input muss char sein!');
end
%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));

Name_Typ = 0; % 0=Fehler, 1=Kinematik, 2=Aktuierung
ActNr = -1; % Initialisierung als "nicht belegt"
Coupling = [-1, -1]; % Initialisierung
% Namensschema für symmetrische, serielle PKM als regulären Ausdruck
% Nehme nicht parroblib_load_robot um beide Arten von Namen zuzulassen
% TODO: Anpassung an nicht-serielle, nicht-symmetrische PKM
expression_kin = 'P(\d)([RP]+)(\d+)[V]?(\d*)G(\d+)P(\d+)'; % Format "P3RRR1G1P1" oder "P3RRR1V1G1P1"
expression_act = [expression_kin,'A(\d+)']; % Format "P3RRR1A1"
[tokens_kin, ~] = regexp(PName_Input,[expression_kin,'$'],'tokens','match');
[tokens_act, ~] = regexp(PName_Input,[expression_act,'$'],'tokens','match');
if ~isempty(tokens_kin)
  res = tokens_kin{1};
  Name_Typ = 1;
  Coupling = [str2double(res{5}), str2double(res{6})];
elseif ~isempty(tokens_act)
  res = tokens_act{1};
  Name_Typ = 2;
  Coupling = [str2double(res{5}), str2double(res{6})];
  ActNr = str2double(res{7});
else
  error('Eingegebener Name %s entspricht nicht dem Namensschema', PName_Input);
end
if isempty(res{4}) % serielle Kette ist keine abgeleitete Variante
  PName_Legs = ['P', res{1}, res{2}, res{3}];
else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
  PName_Legs = ['P', res{1}, res{2}, res{3}, 'V', res{4}];
end
PName_Kin = [PName_Legs, 'G', res{5}, 'P', res{6}];
Coupling = [str2double(res{5}), str2double(res{6})];
NLEG = str2double(res{1});

%% Tabelle für Kinematik öffnen und Zeile entfernen
% Prüfen, in welcher Tabelle der Roboter ist
EEFG_Ges = logical(...
  [1 1 0 0 0 0; 1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
EEstr = ''; % Platzhalter, wird im folgenden belegt.
for jj = 1:size(EEFG_Ges,1)
  if sum(EEFG_Ges(jj,:)) ~= NLEG, continue; end % PKM-FG passen nicht zu Beinketten
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(jj,1:3)), sum(EEFG_Ges(jj,4:6)));
  kintabfile=fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_kin.mat']);
  if ~exist(kintabfile, 'file')
    error('Datei %s existiert nicht. parroblib_gen_bitarrays ausführen!', kintabfile);
  end
  tmp = load(kintabfile);
  KinTab = tmp.KinTab;
  if any(strcmp(KinTab.Name, PName_Kin))
    break; % Der Roboter ist in der aktuellen Tabelle. Variable EEstr wird übernommen
  end
end
  

if Name_Typ == 1
  kintabfile = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list.csv']);
  kintabfile_copy = [kintabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung

  fid = fopen(kintabfile, 'r');
  if fid == -1
    warning('Tabelle %s konnte nicht geöffnet werden', kintabfile);
    success = false;
    return
  end
  fidc = fopen(kintabfile_copy, 'w');
  tline = fgetl(fid);
  found = false;
  while ischar(tline) && ~isempty(tline)
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
    if strcmp(csvline{1}, PName_Input)
      % Zu löschenden Roboter gefunden. Diese Zeile nicht in Dateikopie schreiben
      found = true;
    else
      % Zeile in die Dateikopie schreiben
      fwrite(fidc, [tline, newline()]);
    end
    tline = fgetl(fid); % nächste Zeile
  end
  fclose(fid);
  fclose(fidc);
  if ~found
    success = false;
    warning('Zu löschendes Modell %s nicht in %s gefunden', PName_Input, kintabfile);
    return
  end
  % Modifizierte Tabelle zurückkopieren
  copyfile(kintabfile_copy, kintabfile);
  % Kopie-Tabelle löschen
  delete(kintabfile_copy);
end
%% Tabelle für Aktuierung öffnen und Zeile entfernen
if Name_Typ == 2
  acttabfile = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
  acttabfile_copy = [acttabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung

  fid = fopen(acttabfile, 'r');
  if fid == -1
    warning('Tabelle %s konnte nicht geöffnet werden', acttabfile);
    success = false;
    return
  end
  fidc = fopen(acttabfile_copy, 'w');
  tline = fgetl(fid);
  found = false;
  lines_written = 0; % Anzahl der geschriebenen Zeilen für diese Führungs-Beinkette
  lines_written_kin = 0; % Anzahl der geschriebenen Zeilen mit dieser Kinematik (inkl G/P-Nummer)
  while ischar(tline) && ~isempty(tline)
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
    if strcmp(csvline{1}, PName_Input)
      % Zu löschenden Roboter gefunden. Diese Zeile nicht in Dateikopie schreiben
      found = true;
    else
      % Zeile in die Dateikopie schreiben
      fwrite(fidc, [tline, newline()]);
      lines_written = lines_written + 1;
      if contains(csvline{1}, PName_Kin)
        lines_written_kin = lines_written_kin + 1;
      end
    end
    tline = fgetl(fid); % nächste Zeile
  end
  fclose(fid);
  fclose(fidc);
  if ~found
    success = false;
    warning('Zu löschendes Modell %s nicht in %s gefunden', PName_Input, acttabfile);
    return
  end
  % Modifizierte Tabelle zurückkopieren
  copyfile(acttabfile_copy, acttabfile);
  % Kopie-Tabelle löschen
  delete(acttabfile_copy);
  
  % Falls jetzt keine Aktuierungen mehr für das Kinematik-Modell da sind,
  % kann der entsprechende Eintrag auch gelöscht werden
  removed_Kin = false;
  if lines_written == 1
    % Es wurde nur die Überschrift in die Datei geschrieben. Keine
    % Aktuierung mehr übrig.
    delete(acttabfile);
  end
  if lines_written_kin == 0
    removed_Kin = parroblib_remove_robot(PName_Kin);
  end
end

%% Restliche Dateien entfernen
if Name_Typ == 1
  % Ordner mit Code und Parameter-Modellen löschen (auch alle Aktuierungen)
  % Aktuell muss noch manuell geprüft werden, ob dabei wichtige Daten verloren gehen
  robdirs = dir(fullfile( repopath, ['sym_', EEstr], PName_Legs, ...
    sprintf('hd_G%dP%d*',Coupling(1),Coupling(2)) ));
  for i = 1:length(robdirs)
    robdir = fullfile(robdirs(i).folder, robdirs(i).name);
    if nofiledelete
      movefile(robdir, [robdir, '_delete']);
    else
      rmdir(robdir, 's');
    end
  end
  % Prüfe, ob die Tabelle "actuation.csv" gelöscht werden kann
  acttabfile = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
  fid = fopen(acttabfile, 'r');
  if fid ~= -1
    tmp = textscan(fid, '%s','delimiter','\n');
    numlines = size(tmp{1},1);
    fclose(fid);
    if numlines < 2 && ~nofiledelete
      % Falls nur die Überschrift steht ist die Länge 0. Bei zweiter Zeile ist sie 2.
      delete(acttabfile);
    end
  end
  % Falls der Ordner jetzt leer ist, wird er auch gelöscht
  tpldir = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'tpl');
  if ~exist(acttabfile, 'file') && ... % gelöschte PKM war die letzte. Aktuierungs-Datei ist gelöscht
      exist(tpldir, 'file') % es existieren Vorlagen-Funktionen, die gelöscht werden können.
    rmdir(tpldir, 's');
  end
  PName_Dir = fullfile(repopath, ['sym_', EEstr], PName_Legs);
  PName_DirContent = dir(fullfile(PName_Dir, '*'));
  PName_DirContent=PName_DirContent(~ismember({PName_DirContent.name},{'.','..'}));
  if isempty(PName_DirContent)
    rmpath_genpath(fullfile(PName_Dir));
    rmdir(fullfile(PName_Dir));
  end
elseif Name_Typ == 2 && ~removed_Kin
  % Aktuierung: Lösche nur den Unterordner, der zu dieser Aktuierung gehört
  % Nur löschen, wenn der Hauptordner des Kinematikmodells noch da ist.
  codedir = fullfile(repopath, ['sym_', EEstr], PName_Legs, ...
    sprintf('hd_G%dP%dA%s', Coupling(1), Coupling(2), ActNr));
  if exist(codedir, 'file') % Das Verzeichnis wird erst durch Code-Generierung angelegt. Ist eventuell noch nicht erfolgt.
    if nofiledelete
      movefile(codedir, [codedir, '_delete']);
    else
      rmdir(codedir, 's');
    end
  end
end
success = true;
return
