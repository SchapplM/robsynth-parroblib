% Suche parallele Roboter mit Beinketten, die es nicht mehr gibt. 
% Dieser Fall kann bei Löschungen in der SerRobLib auftreten.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-01
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc
serroblibpath=fileparts(which('serroblib_path_init.m'));
%% Durchsuche ParRobLib
Delete_List = {};
for N = 3:6
  fprintf('Prüfe PKM mit %d FG\n', N);
  % Lade Namen aller möglicher Beinketten (mit verschiedener FG-Anzahl
  l = struct('Names_Ndof', []);
  for i = N:6 
    mdllistfile_Ndof = fullfile(serroblibpath, sprintf('mdl_%ddof', i), sprintf('S%d_list.mat',i));
    l_tmp = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo', 'BitArrays_EEdof0');
    l.Names_Ndof = [l.Names_Ndof,l_tmp.Names_Ndof];
  end
  % Lade alle Roboter die es gibt (keine Filter wirksam)
  EEFG_Ges = logical(...
    [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
     1 1 1 1 1 0; 1 1 1 1 1 1]);
  PNames_Kin = {};
  for j = 1:size(EEFG_Ges,1)
    PNames_Kin = [PNames_Kin, parroblib_filter_robots(EEFG_Ges(j,:))]; %#ok<AGROW>
  end
  for i = 1:length(PNames_Kin)
    PName_kin_i = PNames_Kin{i};
    % Bestimme Namen der Beinkette
    [~, LEG_Names] = parroblib_load_robot([PName_kin_i,'A0']);
    % Prüfe, ob die Beinketten in der SerRobLib existieren
    I_srl = strcmp(LEG_Names{1}, l.Names_Ndof);
    if ~any(I_srl)
      % Beinkette existiert nicht
      fprintf('%s: Serielle Beinkette %s existiert nicht\n', PName_kin_i, LEG_Names{1});
      Delete_List = [Delete_List(:); PName_kin_i];
    end
  end
end

%% PKM löschen
fprintf('%d PKM haben keine Beinkette mehr in der Datenbank. Entfernung manuell!\n', length(Delete_List));
disp(Delete_List);
return

for i = 1:length(Delete_List)
  PName = Delete_List{i};
  success = parroblib_remove_robot(PName);
  if ~success
    warning('Löschen fehlgeschlagen');
  end
end