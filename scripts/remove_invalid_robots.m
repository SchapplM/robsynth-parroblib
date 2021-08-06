% Entferne PKM, die ungültig sind, z.B. weil die Angabe "Rangverlust" fehlt

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
repopath=fileparts(which('parroblib_path_init.m'));
PKM_List_invalid = {};
EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
% Aktualisiere .mat-Dateien. Vor Aufruf dieses Skripts wird üblicherweise
% die Struktursynthese durchgeführt und die Datenbank hat sich geändert.
parroblib_gen_bitarrays();
for j = 1:size(EEFG_Ges,1)
  fprintf('Prüfe PKM mit FG %dT%dR Beinketten\n', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  [~, PNames_Akt] = parroblib_filter_robots(EEFG_Ges(j,1:6));
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  acttabfile=fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_act.mat']);
  tmp = load(acttabfile); % siehe parroblib_gen_bitarrays
  ActTab = tmp.ActTab;
  I_invalid = isnan(ActTab.Rankloss_Platform);
  PKM_List_invalid = [PKM_List_invalid(:)', ActTab.Name(I_invalid)'];
end
fprintf('Lösche %d PKM:\n', length(PKM_List_invalid));
disp(PKM_List_invalid);
for i = 1:length(PKM_List_invalid)
  fprintf('Lösche %d/%d: %s\n', i, length(PKM_List_invalid), PKM_List_invalid{i});
  success = parroblib_remove_robot(PKM_List_invalid{i});
  if ~success
    error('Fehler beim Löschen von %s', PKM_List_invalid{i});
  end
end
parroblib_gen_bitarrays();