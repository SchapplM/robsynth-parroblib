% Ändere die zusätzlich gespeicherten Eigenschaften von PKM in der Datenbank
% Es wird die Datei actuation.csv des jeweiligen Roboters bearbeitet.
% 
% Eingabe:
% PName_Akt
%   Name des mit Aktuierung angegebenen Robotermodells
% varargin
%   Name und Wert für zusätzliche PKM-Eigenschaften:
%   * rankloss: Anzahl der Dimensionen des Rangverlustes

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-09
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_change_properties(PName_Akt, varargin)
p = inputParser;
addRequired(p,'PName_Akt');
addParameter(p,'rankloss',[])
parse(p,PName_Akt,varargin{:});

repopath=fileparts(which('parroblib_path_init.m'));
[NLEG, ~, ~, ~, ~, ~, ~, ~, PName_Legs] = parroblib_load_robot(PName_Akt);

%% Durchsuche die Aktuierungstabelle und ändere die entsprechende Zeile
acttabfile = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Legs, 'actuation.csv');
acttabfile_copy = [acttabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung
fid = fopen(acttabfile, 'r');
fidc = fopen(acttabfile_copy, 'w');
tline = fgetl(fid);
while ischar(tline)
  csvline = strsplit(tline, ';', 'CollapseDelimiters', false); % Spaltenweise als Cell-Array
  if isempty(csvline) || strcmp(csvline{1}, '')
    % ursprüngliche Zeile direkt schreiben
  elseif strcmp(csvline{1}, PName_Akt)
    % Zeile gehört zum gesuchten Roboter. Verändere die Eigenschaft und
    % speichere die Zeile
    if ~isempty(p.Results.rankloss)
      csvline{end} = p.Results.rankloss;
    end
    % Zeile neu generieren
    tline = csvline{1};
    for i = 2:length(csvline)
      tline = sprintf('%s;%s', tline, csvline{i});
    end
  else
    % Anderer Roboter, ursprüngliche Zeile direkt schreiben
  end
  fwrite(fidc, [tline, newline]);
  tline = fgetl(fid); % nächste Zeile
end
fclose(fid);
fclose(fidc);
copyfile(acttabfile_copy, acttabfile); % Modifizierte Tabelle zurückkopieren
delete(acttabfile_copy); % Kopie-Tabelle löschen

