% Entferne Roboter aus der Datenbank, die auf Beinketten aufbauen, die
% aktualisiert wurden (z.B. bezüglich Basis-Orientierung)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

% Modus für Skript: Entweder löschen ('delete'), oder Dateien
% aktualisieren ('update')
mode = 'delete';

% Vorher PKM-Datenbank (mat-Datei) aktualisieren, damit Vorgang auch bei
% mehrfacher und unterbrochenr Durchführung konsistent ist.
parroblib_gen_bitarrays();
% Lade Liste mit Beinketten, die aktualisiert wurden. Wird in SerRobLib
% Mit Skript remove_base_alignment_from_MDH.m erzeugt.
papath = fileparts(which('robsynth_projektablage_path.m'));
listfile = fullfile(papath, '03_Entwicklung', 'Struktursynthese', ...
  '20220121_Aktualisierung_Beinketten_Basisausrichtung', ...
  'base_alignment_changed_list20220120_091120_update.txt');
LegNames_updated = readlines(listfile);

% Gehe alle PKM durch und gleiche die Beinketten gegen die Liste ab
EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
for j = 1:size(EEFG_Ges,1)
  PNames_Kin = parroblib_filter_robots(EEFG_Ges(j,1:6));
  PKM_deleted = {};
  for i = 1:length(PNames_Kin)
    fprintf('%d/%d: Prüfe PKM %s\n', i, length(PNames_Kin), PNames_Kin{i});
    [~, LEG_Names, ~, Coupling, ~, ~, ~, ~, PName_Legs] = parroblib_load_robot(PNames_Kin{i}, 0);
    % Prüfe, ob Name der Beinkette in Liste mit aktualisierten Ketten ist
    I_find = strcmp(LegNames_updated, LEG_Names{1});
    if all(~I_find)
      % Alles in Ordnung: Die Beinkette ist nicht in der Liste
      continue
    end
    if strcmp(mode, 'update'), OutStr = 'Aktualisiere';
    else, OutStr = 'Lösche'; end
    fprintf('Beinkette %s wurde in SerRob-DB geändert. %s PKM %s.\n', ...
      LEG_Names{1}, OutStr, PNames_Kin{i});
    if strcmp(mode, 'update')
      % Lösche aus Vorlagen-Generierte Funktionen
      EE_FG_str = sprintf('%dT%dR', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
      parroblib_path = fileparts(which('parroblib_path_init.m'));
      tpl_dir = fullfile(parroblib_path, sprintf('sym_%s', EE_FG_str), ...
        PName_Legs, 'tpl');
      if exist(tpl_dir, 'file')
        rmdir(tpl_dir, 's');
        fprintf('Verzeichnis %s gelöscht\n', tpl_dir);
      end
      parroblib_update_template_functions(PNames_Kin(i));
    else
      success = parroblib_remove_robot(PNames_Kin{i});
      if ~success
        error('Fehler beim Löschen von %s', PNames_Kin{i});
      else
        % Aktualisiere Struktursynthese-Datenbank (setze auf gelöscht)
        parroblib_update_csv(LEG_Names{1}, Coupling, EEFG_Ges(j,1:6), 9);
      end
    end
  end
end
% Nachher PKM-Datenbank aktualisieren, damit gelöschte in mat-Dateien auch
% entfernt werden.
parroblib_gen_bitarrays();