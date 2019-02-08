% Entferne ein Robotermodell aus der Datenbank (csv-Tabelle)
% Die Zeile des Robotermodells (z.B. P3PRP2) wird aus der csv-Tabelle der
% Gelenkstruktur (z.B. sym3leg.csv) entfernt. Zusätzlich werden die
% Matlab-Funktionen aus dem gleichnamigen Unterordner ("P3PRP2") entfernt.
% 
% Eingabe:
% PName_Kin
%   Name des Roboters in der Datenbank (nur Kinematik, ohne Endung für
%   Aktuierung; z.B. "P3PRP1")
% 
% Ausgabe:
% success
%   true, falls Roboter erfolgreich entfernt wurde
% 
% Siehe auch: serroblib_remove_robot

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

function success = parroblib_remove_robot(PName_Kin)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));

% Namensschema für symmetrische, serielle PKM als regulären Ausdruck
% TODO: Anpassung an nicht-serielle, nicht-symmetrische PKM
expression = 'P(\d)([RP]+)(\d+)'; % Format "P3RRR1"
[tokens, ~] = regexp(PName_Kin,expression,'tokens','match');
if isempty(tokens)
  error('Eingegebener Name %s entspricht nicht dem Namensschema', PName_Kin);
end
res = tokens{1};
NLEG = str2double(res{1});

%% Tabelle öffnen und Zeile entfernen
kintabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
kintabfile_copy = [kintabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung

fid = fopen(kintabfile);
fidc = fopen(kintabfile_copy, 'w');
tline = fgetl(fid);
found = false;
while ischar(tline)
  % Spaltenweise als Cell-Array
  csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
  if strcmp(csvline{1}, PName_Kin)
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
  warning('Zu löschendes Modell %s nicht in %s gefunden', ...
    Name, kintabfile);
  return
end
% Modifizierte Tabelle zurückkopieren
copyfile(kintabfile_copy, kintabfile);
% Kopie-Tabelle löschen
delete(kintabfile_copy);

%% Restliche Dateien entfernen
% Ordner mit Code und Parameter-Modellen löschen
% Aktuell muss noch manuell geprüft werden, ob dabei wichtige Daten verloren gehen
robdir = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin);
rmdir(robdir, 's');

success = true;
return
