% Entferne PKM aus der Kinematik-Liste, für die keine Aktuierung
% gespeichert ist

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2023-03
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear


repopath=fileparts(which('parroblib_path_init.m'));

EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
for j = 1:size(EEFG_Ges,1)
  parroblib_gen_bitarrays(EEFG_Ges(j,1:6));
  [PNames_Kin, PNames_Act] = parroblib_filter_robots(EEFG_Ges(j,1:6));
  for i = 1:length(PNames_Kin)
    if ~any(contains(PNames_Act, PNames_Kin{i}))
      fprintf('Zu Kinematik %d (%s) gibt es keine Aktuierung. Lösche.\n', ...
        i, PNames_Kin{i});
      success = parroblib_remove_robot(PNames_Kin{i});
      if ~success
        error('Fehler beim Löschen von %s', PNames_Kin{j});
      end
    end
  end
end
