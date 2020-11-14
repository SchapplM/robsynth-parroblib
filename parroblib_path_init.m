% Diese Datei dient nur zum Auffinden dieses Ordners mit dem Befehl
% `which('parroblib_path_init.m')` und zum Hinzufügen des 
% Hauptverzeichnisses des Repos

this_tb_path = fileparts( mfilename('fullpath') );
addpath(this_tb_path);

% Alle csv-Tabellen (versionsverwaltet) nach .mat konvertieren (nicht
% versionsverwaltet, da binär-Format). Nur machen, wenn Dateien fehlen.
% Falls Dateien veraltet sind, muss das manuell durchgeführt werden.
% Dauert relativ lange.
for NLEG = 3:6
  kintabfile_mat = fullfile(this_tb_path, sprintf('sym%dleg', NLEG), ...
    sprintf('sym%dleg_list_kin.mat', NLEG));
  acttabfile_mat = fullfile(this_tb_path, sprintf('sym%dleg', NLEG), ...
    sprintf('sym%dleg_list_act.mat', NLEG));
  if ~exist(kintabfile_mat, 'file') || ~exist(acttabfile_mat, 'file')
    fprintf('mat-Datei für PKM mit %d Beinketten existiert nicht. Erstelle.\n', NLEG);
    parroblib_gen_bitarrays(NLEG);
  end
end
