% Entferne PKM, die ungültig sind, z.B. weil die Angabe "Rangverlust" fehlt

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
PKM_List_invalid = {};
EEFG_Ges = logical(...
  [1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
   1 1 1 1 1 0; 1 1 1 1 1 1]);
for j = 1:size(EEFG_Ges,1)
  fprintf('Prüfe PKM mit FG %dT%dR Beinketten\n', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  [~, PNames_Akt] = parroblib_filter_robots(EEFG_Ges(j,1:6));

  for i = 1:length(PNames_Akt)
    fprintf('%d/%d: Prüfe PKM %s\n', i, length(PNames_Akt), PNames_Akt{i});
    [~, LEG_Names, Actuation, Coupling, ActNr, symrob, EE_dof0, ...
      PName_Kin, PName_Legs, AdditionalInfo_Akt] = parroblib_load_robot(PNames_Akt{i});
    if isnan(AdditionalInfo_Akt(1))
      warning('Rangverlust hat keinen abgespeicherten Zahlenwert für %s', PNames_Akt{i});
      PKM_List_invalid = [PKM_List_invalid(:)', PNames_Akt{i}];
    end
  end
end
fprintf('Lösche %d PKM:\n', length(PKM_List_invalid));
for i = 1:length(PKM_List_invalid)
  fprintf('Lösche %d/%d: %s\n', i, length(PKM_List_invalid), PKM_List_invalid{i});
  success = parroblib_remove_robot(PKM_List_invalid{i});
  if ~success
    error('Fehler beim Löschen von %s', PKM_List_invalid{i});
  end
end
