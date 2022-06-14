% Suche Isomorphismen in der PKM-Datenbank und entferne sie. Meistens fällt
% erst später nach der Struktursynthese auf, dass manche Strukturen gleich
% aussehen. Fälle von Isomorphismen:
% 1: Gestell-Modus 1 und 10 bei Beinketten, die mit PR beginnen

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

iso_case = 1;

serroblibpath=fileparts(which('serroblib_path_init.m'));
serroblist = load(fullfile(serroblibpath, 'serrob_list.mat'));

EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
for j = 1:size(EEFG_Ges,1)
  PNames_Kin = parroblib_filter_robots(EEFG_Ges(j,1:6));
  for i = 1:length(PNames_Kin)
    match = false;
    [~, LEG_Names, ~, Coupling] = parroblib_load_robot([PNames_Kin{i},'A0'], 0);
    Chain = LEG_Names{1};
    I_Legchain = strcmp(serroblist.Names, Chain);
    % Prüfe Fall 1: 
    if iso_case == 1
      if Coupling(1)~=10, continue; end % Suche nur G10-PKM
      if ~strcmp(Chain(3:4), 'PR'), continue; end % Suche nur mit PR-Kette
      % Prüfe die Parallelität der Gelenke
    end
    [~,PS] = serroblib_create_robot_class(Chain, '', true);
    if iso_case == 1
      if PS.alpha(2) ~= 0
        % a-Parameter wird nicht betrachtet, da er in der Maßsynthese auf
        % Null gesetzt wird (entspricht Plattform-Größe)
        continue; % Keine Parallelität von P und R. Kein Isomorphismus
      end
    end
    % Bis hier gekommen. PKM löschen.
    success = parroblib_remove_robot(PNames_Kin{i});
    assert(success, sprintf('Fehler beim Löschen von %s', PNames_Kin{i}));
    fprintf('PKM %s gelöscht (Isomorphismus)\n', PNames_Kin{i});
  end
end