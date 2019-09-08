% Entferne ein Robotermodell aus der Datenbank (csv-Tabelle)
% Die Zeile des Robotermodells (z.B. P3PRP2) wird aus der csv-Tabelle der
% Gelenkstruktur (z.B. sym3leg.csv) entfernt. Zusätzlich werden die
% Matlab-Funktionen aus dem gleichnamigen Unterordner ("P3PRP2") entfernt.
% 
% Eingabe:
% PName
%   Name des Roboters in der Datenbank. Optionen:
%   * PName_Kin (nur Kinematik, ohne Endung für Aktuierung; z.B. "P3RPR1")
%   * PName_Act (Kinematik mit Aktuierung; z.B. "P3RPR1A2")
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

% Namensschema für symmetrische, serielle PKM als regulären Ausdruck
% TODO: Anpassung an nicht-serielle, nicht-symmetrische PKM
expression_kin = 'P(\d)([RP]+)(\d+)[V]?(\d*)'; % Format "P3RRR1" oder "P3RRR1V1"
expression_act = [expression_kin,'A(\d+)']; % Format "P3RRR1A1"
[tokens_kin, ~] = regexp(PName_Input,[expression_kin,'$'],'tokens','match');
[tokens_act, ~] = regexp(PName_Input,[expression_act,'$'],'tokens','match');
if ~isempty(tokens_kin)
  PName_Kin = PName_Input;
  Name_Typ = 1;
  res = tokens_kin{1};
elseif ~isempty(tokens_act)
  res = tokens_act{1};
  if isempty(res{4}) % serielle Kette ist keine abgeleitete Variante
    PName_Kin = ['P', res{1}, res{2}, res{3}];
  else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
    PName_Kin = ['P', res{1}, res{2}, res{3}, 'V', res{4}];
  end
  Name_Typ = 2;
else
  error('Eingegebener Name %s entspricht nicht dem Namensschema', PName_Input);
end
NLEG = str2double(res{1});

%% Tabelle für Kinematik öffnen und Zeile entfernen
if Name_Typ == 1
  kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
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
  while ischar(tline)
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
  acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, 'actuation.csv');
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
  lines_written = 0;
  while ischar(tline)
    % Spaltenweise als Cell-Array
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
    if strcmp(csvline{1}, PName_Input)
      % Zu löschenden Roboter gefunden. Diese Zeile nicht in Dateikopie schreiben
      found = true;
    else
      % Zeile in die Dateikopie schreiben
      fwrite(fidc, [tline, newline()]);
      lines_written = lines_written + 1;
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
    removed_Kin = parroblib_remove_robot(PName_Kin);
  end
end

%% Restliche Dateien entfernen
if Name_Typ == 1
  % Ordner mit Code und Parameter-Modellen löschen (auch alle Aktuierungen)
  % Aktuell muss noch manuell geprüft werden, ob dabei wichtige Daten verloren gehen
  robdir = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Input);
  if nofiledelete
    move(robdir, [robdir, '_delete']);
  else
    rmdir(robdir, 's');
  end
elseif Name_Typ == 2 && ~removed_Kin
  % Aktuierung: Lösche nur den Unterordner, der zu dieser Aktuierung gehört
  % Nur löschen, wenn der Hauptordner des Kinematikmodells noch da ist.
  actnr = PName_Input(length(PName_Kin)+2:end);
  codedir = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, sprintf('hd_A%s', actnr));
  if exist(codedir, 'file') % Das Verzeichnis wird erst durch Code-Generierung angelegt. Ist eventuell noch nicht erfolgt.
    if nofiledelete
      move(codedir, [codedir, '_delete']);
    else
      rmdir(codedir, 's');
    end
  end
end
success = true;
return
