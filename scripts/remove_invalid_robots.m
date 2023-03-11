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

%% Suche nach ungültigen Aktuierungen
for j = 1:size(EEFG_Ges,1)
  fprintf('Prüfe PKM mit FG %dT%dR\n', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  [~, PNames_Akt] = parroblib_filter_robots(EEFG_Ges(j,1:6));
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  acttabfile=fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_act.mat']);
  tmp = load(acttabfile); % siehe parroblib_gen_bitarrays
  ActTab = tmp.ActTab;
  I_invalid = isnan(ActTab.Rankloss_Platform);
  PKM_List_invalid = [PKM_List_invalid(:)', ActTab.Name(I_invalid)'];
end

%% Suche auch doppelte Aktuierungs-Nummern kann fälschlicherweise auftreten
PKM_List_invalid = {};
for j = 1:size(EEFG_Ges,1)
  fprintf('Prüfe PKM mit FG %dT%dR\n', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  EEstr = sprintf('%dT%dR', sum(EEFG_Ges(j,1:3)), sum(EEFG_Ges(j,4:6)));
  acttabfile=fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_act.mat']);
  tmp = load(acttabfile); % siehe parroblib_gen_bitarrays
  ActTab = tmp.ActTab;
  PNames_Kin = parroblib_filter_robots(EEFG_Ges(j,1:6));
  for ii = 1:length(PNames_Kin)
    NLEG = parroblib_load_robot(PNames_Kin{ii}, 0);
    ii_act = contains(ActTab.Name, PNames_Kin{ii});
    % Prüfe, ob die Einträge Duplikate enthalten
    if sum(ii_act) == 1, continue; end
    for k = fliplr(find(ii_act)') % bei letzten Einträgen anfangen
%       fprintf('Prüfe Aktuierung %s\n', ActTab.Name{k});
      for l = find(ii_act)'
        if l == k, continue; end % nicht mit sich selbst prüfen
        % nicht gegen bereits entfernte PKM prüfen. Sonst werden beide
        % doppelte entfernt statt nur eine. Falls Einträge dreifach
        % vorkommen, muss die Prüfung mehrfach laufen
        if any(strcmp(PKM_List_invalid, ActTab.Name{k}))
          continue;
        end
        if any(strcmp(PKM_List_invalid, ActTab.Name{l}))
          continue; % ist schon auf Lösch-Liste
        end
        duplicate = true; % Hypothese: l ist Duplikat von k
        for m = 1:NLEG
          Col_m = ActTab.(sprintf('Act_Leg%d', m));
          if any(Col_m(l,:) ~= Col_m(k,:)) % Aktuierung ist unterschiedlich
            duplicate = false; % Hypothese widerlegt
            break;
          end
        end
        if duplicate
          fprintf('Duplikat: %s vs %s\n', ActTab.Name{k}, ActTab.Name{l});
          ActTab_dupl = ActTab(ii_act,:);
          PKM_List_invalid = [PKM_List_invalid(:)', ActTab.Name{k}];
        end
      end % for l
    end % for k
  end % for ii
end

%% Änderungen übernehmen
return
fprintf('Lösche %d PKM:\n', length(PKM_List_invalid));
disp(PKM_List_invalid);
for NLEG = 1:6 % gehe die Einträge gemäß der Anzahl Beinketten durch
  % Finde alle Einträge zu dieser Anzahl Beinketten
  II = find(~cellfun(@isempty, regexp(PKM_List_invalid, sprintf('(^P%d)', NLEG), 'match')));
  for i = 1:length(II)
    fprintf('Lösche %d/%d: %s\n', i, length(II), PKM_List_invalid{II(i)});
    success = parroblib_remove_robot(PKM_List_invalid{II(i)}, false, ...
      i==length(II)); % nur beim letzten CSV aktualisieren
    if ~success
      continue
      error('Fehler beim Löschen von %s', PKM_List_invalid{i});
    end
  end
end
parroblib_gen_bitarrays();
