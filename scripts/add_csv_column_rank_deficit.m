% Füge Spalten für Mobilitätsverlust in Aktuierungs-Tabellen hinzu
% (actuation.csv)
% 
% Setze voraus, dass die neue Spalte noch in keiner Tabelle drin ist.
% 
% Siehe auch: parroblib_add_robot.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-09
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

repopath=fileparts(which('parroblib_path_init.m'));
actfilelist = dir(fullfile(repopath, '**', 'actuation.csv'));

for i = 1:length(actfilelist)
  acttabfile = fullfile(actfilelist(i).folder, actfilelist(i).name);
  fprintf('Bearbeite %s\n', acttabfile);
  acttabfile_copy = [acttabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung
  fid = fopen(acttabfile, 'r');
  fidc = fopen(acttabfile_copy, 'w');
  tline = fgetl(fid);
  iline=0;
  while ischar(tline)
    iline=iline+1;
    if iline == 1
      tline=[tline, ';Rangverlust Plattform-FG']; %#ok<AGROW>
    else
      % Gehe davon aus, dass alle bestehenden Systeme keinen Rangverlust
      % haben
      tline=[tline, ';0']; %#ok<AGROW>
    end
    fwrite(fidc, [tline, newline]);
    tline = fgetl(fid); % nächste Zeile
  end
  fclose(fid);
  fclose(fidc);
  copyfile(acttabfile_copy, acttabfile);
  delete(acttabfile_copy);
end