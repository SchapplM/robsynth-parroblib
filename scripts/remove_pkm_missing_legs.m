% Entferne Roboter aus der Datenbank, die auf Beinketten aufbauen, die
% mittlerweile gelöscht wurden (z.B. da sie Isomorphismen sind)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-02
% (C) Institut für Mechatronische Systeme, Universität Hannover

% Falls direkt vorher die Beinketten entfernt wurden:
serroblib_gen_bitarrays();

clc
clear

serroblibpath=fileparts(which('serroblib_path_init.m'));
parroblibpath=fileparts(which('parroblib_path_init.m'));
EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
for j = 1:size(EEFG_Ges,1)
  PNames_Kin = parroblib_filter_robots(EEFG_Ges(j,1:6));

  for i = 1:length(PNames_Kin)
    fprintf('%d/%d: Prüfe PKM %s\n', i, length(PNames_Kin), PNames_Kin{i});
    [~, LEG_Names] = parroblib_load_robot([PNames_Kin{i},'A0'], 0);
    % Prüfe, ob Name der Beinkette in Modellbibliothek serieller Roboter ist
    NFG_Leg = str2double(LEG_Names{1}(2));
    mdllistfile_Ndof = fullfile(serroblibpath, sprintf('mdl_%ddof', NFG_Leg), sprintf('S%d_list.mat',NFG_Leg));
    l = load(mdllistfile_Ndof, 'Names_Ndof');
    if any(strcmp(l.Names_Ndof,LEG_Names{1}))
      % Alles in Ordnung: Die Beinkette gibt es
      continue
    end
    fprintf('Beinkette %s nicht in Datenbank serieller Roboter. Lösche PKM %s.\n', LEG_Names{1}, PNames_Kin{i});
    success = parroblib_remove_robot(PNames_Kin{i});
    if ~success
      error('Fehler beim Löschen von %s', PNames_Kin{i});
    end
  end
  % Gehe durch Liste der Synthese-Ergebnisse
  for k = 1:2
    if k == 1, varstr = '';
    else,      varstr = '_var';
    end
    synthresfile = fullfile(parroblibpath, 'synthesis_result_lists', ...
      sprintf('%dT%dR%s.csv', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)), varstr));
    synthrestab = readtable(synthresfile);
    I_del = false(size(synthrestab,1), 1);
    for i = 1:size(synthrestab,1)
      if ~any(strcmp(l.Names_Ndof, synthrestab.Beinkette{i}))
        I_del(i) = true;
      end
    end
    fprintf('Lösche %d Einträge aus Tabelle %s\n', sum(I_del), synthresfile);
    if sum(I_del) > 0
      synthrestab = synthrestab(~I_del,:);
      writetable(synthrestab, synthresfile, 'Delimiter', ';');
    end
  end
end
