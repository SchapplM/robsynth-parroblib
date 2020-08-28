% Entferne PKM, die ungültig sind, z.B. weil die Angabe "Rangverlust" fehlt

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
PKM_List_invalid = {};
for NLEG = [3 4 6]
  fprintf('Prüfe PKM mit %d Beinketten\n', NLEG);
  EE_FG0 = [0 0 0 0 0 0];
  EE_FG_Mask = [0 0 0 0 0 0]; % Maske 0, EE-FG sind egal.
  [PNames_Kin, PNames_Akt] = parroblib_filter_robots(NLEG, EE_FG0, EE_FG_Mask);

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
