% Ändere die zusätzlich gespeicherten Eigenschaften von PKM in der Datenbank
% Es wird die Datei actuation.csv des jeweiligen Roboters bearbeitet.
% 
% Eingabe:
% PName_Akt
%   Name des mit Aktuierung angegebenen Robotermodells
% varargin
%   Name und Wert für zusätzliche PKM-Eigenschaften:
%   * rankloss: Anzahl der Dimensionen des Rangverlustes
%   * values_angles: Mögliche Werte für die freien Winkel (alpha und theta).
%     Kodierung durch Komma-getrennte Zeichen für jeden möglichen Zustand.
%     Reihenfolge: Aufsteigende Gelenknummer, erst alpha, dann theta
%     Beispiel: z.B. "op,po". Siehe parroblib_load_robot.
%   * joint_parallelity: Parallelität der Gelenke innerhalb der Beinketten

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-09
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_change_properties(PName_Akt, varargin)
p = inputParser;
addRequired(p,'PName_Akt');
addParameter(p,'rankloss',[])
addParameter(p,'values_angles',[])
addParameter(p,'joint_parallelity',[])
parse(p,PName_Akt,varargin{:});
if any(strcmp(varargin, 'values_angles'))
  values_angles_set = true; % Sonst nicht erkennbar, falls leerer String übergeben wird
else
  values_angles_set = false;
end
repopath=fileparts(which('parroblib_path_init.m'));
[~, ~, ~, ~, ~, ~, EEdof0, PName_Kin, PName_Legs] = parroblib_load_robot(PName_Akt);
EEstr = sprintf('%dT%dR', sum(EEdof0(1:3)), sum(EEdof0(4:6)));
%% Durchsuche die Kinematik-Tabelle und ändere die entsprechende Zeile
if ~isempty(p.Results.joint_parallelity)
  kintabfile = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list.csv']);
  T = readtable(kintabfile, 'Delimiter', ';');
  assert(all(size(p.Results.joint_parallelity) <= [1 6]), 'Dimension der Eingabe joint_parallelity stimmt nicht');
  T.Gelenkgruppen{strcmp(T.Name, PName_Kin)} = strrep(disp_array( ...
    p.Results.joint_parallelity, '%d'), ', ', '-'); % Format: 1-1-2-2-3
  writetable(T, kintabfile, 'Delimiter', ';');
end
%% Durchsuche die Aktuierungstabelle und ändere die entsprechende Zeile
if ~isempty(p.Results.rankloss) || values_angles_set
  acttabfile = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
  acttabfile_copy = [acttabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung
  fid = fopen(acttabfile, 'r');
  fidc = fopen(acttabfile_copy, 'w');
  tline = fgetl(fid);
  while ischar(tline) && ~isempty(tline)
    csvline = strsplit(tline, ';', 'CollapseDelimiters', false); % Spaltenweise als Cell-Array
    if strcmp(csvline{1}, '')
      % ursprüngliche Zeile direkt schreiben
    elseif strcmp(csvline{1}, PName_Akt)
      % Zeile gehört zum gesuchten Roboter. Verändere die Eigenschaft und
      % speichere die Zeile
      if ~isempty(p.Results.rankloss)
        csvline{end-1} = p.Results.rankloss;
      end
      if values_angles_set
        csvline{end} = p.Results.values_angles;
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
end
