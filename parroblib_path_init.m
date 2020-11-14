% Diese Datei dient nur zum Auffinden dieses Ordners mit dem Befehl
% `which('parroblib_path_init.m')` und zum Hinzufügen des 
% Hauptverzeichnisses des Repos

this_tb_path = fileparts( mfilename('fullpath') );
addpath(this_tb_path);

% Alle csv-Tabellen (versionsverwaltet) nach .mat konvertieren (nicht
% versionsverwaltet, da binär-Format). Nur machen, wenn Dateien fehlen.
% Falls Dateien veraltet sind, muss das manuell durchgeführt werden.
% Dauert relativ lange.
EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
for jj = 1:size(EEFG_Ges,1)
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(jj,1:3)), sum(EEFG_Ges(jj,4:6)));
  kintabfile_mat = fullfile(this_tb_path, ['sym_', EEstr], ...
    ['sym_',EEstr,'_list_kin.mat']);
  acttabfile_mat = fullfile(this_tb_path, ['sym_', EEstr], ...
    ['sym_',EEstr,'_list_act.mat']);
  if ~exist(kintabfile_mat, 'file') || ~exist(acttabfile_mat, 'file')
    fprintf('mat-Datei für PKM mit %s FG existiert nicht. Erstelle:\n\t%s\n\t%s\n', ...
      EEstr, kintabfile_mat, acttabfile_mat);
    parroblib_gen_bitarrays(EEFG_Ges(jj,:));
  end
end
