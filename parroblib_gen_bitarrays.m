% Generiere .mat-Dateien zum schnelleren Durchsuchen der Datenbank
% 
% Eingabe:
% EEFG_update (optional)
%   Array mit EE-FG, für die neu generiert werden soll. Zeilenweise FG im
%   Format [1 1 0 0 0 1]
% 
% Ausgabe:
% KinTab_All (optional)
%   Tabelle mit Kinematik-Eigenschaften aller PKM in der Datenbank
% 
% Schreibt Dateien:
% sym_xTyR_list_kin.mat (x,y aus EEFG aus Eingabe). Enthält Variable
% KinTab. Entspricht symxleg_list.csv. Matlab-Tabelle mit Feldern:
%   Name
%     Name der PKM-Kinematik (mit Koppelgelenk): "P6PRRRRR1G1P5"
%   Beinkette
%     Name der seriellen Beinkette (siehe serroblib): "S6PRRRRR1"
%   Beinkette_Tech
%     Name der technischen Gelenke der seriellen Beinkette: "UPS"
%   EEFG [1x6]
%     Einträge mit 1/0, je nachdem welcher EE-FG bewegt werden kann. FG
%     werden im Basis-KS gezählt (Plattform-Geschw. und -winkelgeschw.)
% parrob_list_kin.mat
%   Enthält alle Dateiinhalte von sym_xTyR_list_kin.mat für alle x/y.
%   
% symxleg_list_act.mat (x=Anzahl Beinketten aus Eingabe). Enthält Variable
% ActTab. Entspricht actuation.csv in den Unterordnern. Matlab-Tabelle:
%   Name
%     Name der aktuierten PKM: "P6PRRRRR1G1P5A1"
%   Act_Leg1 [1xNJ1]
%     Aktuierung der `NJ1` Gelenke von Beinkette 1 (1/0 für aktiv/passiv)
%   Act_Leg2 ...
%     Das gleiche für alle weiteren Beinketten der PKM
%   Rankloss_Platform
%     Rangverlust der Jacobi-Matrix (in den vorgesehenen FG der PKM)
%   Values_Angle_Parameters
%     Mögliche Werte für die freien Winkelparameter der Beinketten. Siehe
%     parroblib_load_robot; Werte 'o', 'a', 'o', 'b' für jeden Parameter.
%   MaxIdxActJoint
%     Nummer des aktuierten Gelenks in allen Gelenk-FG. Eintrag bspw.
%     "3" bei RRPRRR-Beinkette (die auch UPS-Kette ist)
%   MaxIdxActTechJoint
%     Nummer des aktuierten technischen Gelenks. Eintrag bspw. "2" bei UPS.
%     Kardan- und Kugelgelenke zählen als ein technisches Gelenk.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-11
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function KinTab_All = parroblib_gen_bitarrays(EEFG_update)
complete_update = false; % Merker, ob Gesamt-Datei erstellt werden soll
if nargin == 0
  EEFG_update = logical(... % Aktualisiere alle Roboter
    [1 1 0 0 0 0; 1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
     1 1 1 1 1 0; 1 1 1 1 1 1]);
end
if size(EEFG_update,1) == 6
  complete_update = true; % Erstelle auch eine Datei für eine Gesamt-Datenbank
end
if size(EEFG_update,2)~=6
  error('Eingabe EEFG_update muss 0/1-Matrix mit 6 Spalten sein');
end

repopath=fileparts(which('parroblib_path_init.m'));
serroblibpath=fileparts(which('serroblib_path_init.m'));
SerRob_DB_all = load(fullfile(serroblibpath, 'serrob_list.mat'));
KinTab_All = {};
for iFG = 1:size(EEFG_update,1)
  EEstr = sprintf('%dT%dR', sum(EEFG_update(iFG,1:3)), sum(EEFG_update(iFG,4:6)));
  %% Kinematik-Tabelle durchsuchen
  kintabfile_csv = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list.csv']);
  if ~exist(kintabfile_csv, 'file')
    warning('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile_csv);
    continue % trotzdem mit anderen FG weitermachen. Sowieso getrennte Tabellen.
  end
  % Tabelle lesen (und dabei Überschriften nachbessern)
  KinTab = readtable(kintabfile_csv, 'Delimiter', ';');
  % Entferne Zeilen, die ungültige Daten enthalten. Darf eigentlich gar
  % nicht vorkommen. Mögliche Ursache: Fehler bei Tabellenüberschriften
  I_valid = true(size(KinTab,1),1);
  for i = 1:size(KinTab,1)
    expression_kin = 'P(\d)([RP]+)(\d+)[V]?(\d*)';
    [tokens_kin, ~] = regexp(KinTab.Name{i},expression_kin,'tokens','match');
    if isempty(tokens_kin)
      warning('Zeile "%s" in Tabelle %s fehlerhaft. Überspringe.', ...
        KinTab.Name{i}, kintabfile_csv)
      I_valid(i) = false;
    end
  end
  KinTab = KinTab(I_valid,:);
  % Formatiere die Spalte zu den Gelenk-Parallelitäten wieder als Zahlen
  KinTab.Gelenkgruppen(strcmp(KinTab.Gelenkgruppen,'?')) = {NaN};
  for i = 1:size(KinTab,1)
    if isnan(KinTab.Gelenkgruppen{i}), continue; end
    [~, match] = regexp(KinTab.Gelenkgruppen{i}, '(\d)', 'tokens', 'match');
    if isempty(match)
      warning('unbekannter Eintrag für Nr. %d (%s)', i, KinTab.Name{i});
    end
    KinTab.Gelenkgruppen{i} = str2double(match);
  end
  
  % Füge zusätzliche Spalten hinzu, die nicht in der CSV-Datei der PKM-
  % Datenbank stehen (da die Informationen indirekt in der SerRobLib sind).
  % Zeichenkette für technische Gelenke der Beinkette eintragen
  Beinkette_Tech = cell(length(KinTab.Name), 1);
  for k = 1:length(KinTab.Name)
    I_SR_k = strcmp(SerRob_DB_all.Names, KinTab.Beinkette{k});
    SName_TechJoint = fliplr(regexprep(num2str(SerRob_DB_all.AdditionalInfo(I_SR_k,7)), ...
      {'1','2','3','4','5'}, {'R','P','C','U','S'}));
    Beinkette_Tech{k} = SName_TechJoint;
  end
  KinTab = addvars(KinTab, Beinkette_Tech);
  KinTab.Properties.VariableNames(end) = {'Beinkette_Tech'};
  %% Sortiere die Tabelle mit aufsteigenden Nummern
  % Sonst kommt nach 1 die 10 bei den Beinketten-Modellnummern.
  expression_kin = 'P(\d)([RP]+)(\d+)[V]?(\d*)G(\d+)P(\d+)';
  [tokens_kin, ~] = regexp(KinTab.Name,expression_kin,'tokens','match');
  KinTabNameTab = cell2table(cell(size(KinTab,1),6), 'VariableNames', ...
      {'NLEG', 'Chain', 'ChainNum', 'Var', 'G', 'P'});
  for i = 1:size(KinTab,1)
    row_i_str = tokens_kin{i}{1};
    if isempty(row_i_str{4}), row_i_str{4} = '0'; end % sonst NaN und Sortierung falsch.
    KinTabNameTab{i,:} = {str2double(row_i_str{1}), row_i_str{2}, ...
      str2double(row_i_str{3}), str2double(row_i_str{4}), str2double(row_i_str{5}), ...
      str2double(row_i_str{6})};
  end
  [~,I] = sortrows(KinTabNameTab);
  KinTab = KinTab(I,:);
  %% Kinematik-Tabelle speichern
  % Prüfe, ob die aktuell ausgelesene CSV-Datei der gleiche Stand ist wie
  % die binär einzulesende Datei
  kintabfile_mat = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_kin.mat']);
  write_new_kin = false;
  if exist(kintabfile_mat, 'file')
    tmp = load(kintabfile_mat);
    KinTab_old = tmp.KinTab;
    if ~isequaln(KinTab, KinTab_old)
      write_new_kin = true; % Dateiinhalt wird sich ändern. Schreibe neu.
    end
  else
    mkdirs(fileparts(kintabfile_mat)); % Ordner könnte nicht existieren
    write_new_kin = true; % Schreibe mat-Datei neu, sie existiert noch nicht
  end
  if write_new_kin
    save(kintabfile_mat, 'KinTab');
  end
  % Eintragen in Gesamt-Tabelle
  if complete_update
    KinTab_EEcol = KinTab;
    KinTab_EEcol = addvars(KinTab, repmat(EEFG_update(iFG,:), size(KinTab,1), 1));
    KinTab_EEcol.Properties.VariableNames(end) = {'EEFG'};
    % Füge Spalte mit EE-FG hinzu
    if isempty(KinTab_All)
      KinTab_All = KinTab_EEcol;
    else
      KinTab_All = [KinTab_All; KinTab_EEcol]; %#ok<AGROW> 
    end
  end
  %% Menge der PKM reduzieren (G-/P-Nummer wieder entfernen)
  PName_Legs_all = {};
  for i = 1:size(KinTab,1)
    PName_Kin = KinTab.Name{i};
    % Extrahiere Beinketten-Name der PKM (ohne Ausrichtungs-Nummern G/P)
    % Unter diesem Namen ist die Aktuierungstabelle gespeichert
    expression = 'P(\d)([RP]+)(\d+)[V]?(\d*)[G]?(\d*)[P]?(\d*)'; % Format "P3RRRG1P11A1" oder "P3RRR1V1G1P1A1"
    [tokens, ~] = regexp(PName_Kin,expression,'tokens','match');
    if isempty(tokens)
      warning('Eintrag "%s" aus der Tabelle stimmt nicht mit Namensschema überein.', PName_Kin);
      continue
    end
    res = tokens{1};
    if isempty(res{4}) % serielle Kette ist keine abgeleitete Variante
      PName_Legs = ['P', res{1}, res{2}, res{3}];
    else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
      PName_Legs = ['P', res{1}, res{2}, res{3}, 'V', res{4}];
    end
    PName_Legs_all = [PName_Legs_all, PName_Legs]; %#ok<AGROW>
  end
  PName_Legs_unique = unique(PName_Legs_all, 'stable');
  %% Aktuierungs-Tabelle auslesen
  ActTab = [];
  for i = 1:length(PName_Legs_unique)
    PName_Legs = PName_Legs_unique{i};
    acttabfile_csv = fullfile(repopath, ['sym_', EEstr], PName_Legs, 'actuation.csv');
    if ~exist(acttabfile_csv, 'file')
      % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
      warning('Zu Kinematik "%s" gibt es keine Aktuierungstabelle %s', PName_Kin, acttabfile_csv);
      continue
    end
    ActTab_i = readtable(acttabfile_csv, 'NumHeaderLines', 1);
    ActTab_i_headers = readtable(acttabfile_csv, 'PreserveVariableNames', ...
      true, 'ReadVariableNames', true);
    if size(ActTab_i,2) == size(ActTab_i_headers,2) - 1
      % Im Datenbereich der Tabelle ist die letzte Spalte leer. Füge ein
      % leeres Cell-Array ein. Notwendig, da hier sonst im Cell-Array die
      % Kürzel der Werte für freie Parameter stehen.
      ActTab_i = addvars(ActTab_i, cell(size(ActTab_i,1),1), ...
        'NewVariableNames', 'Values_Angle_Parameters');  
    end
    varnames_act = ActTab_i_headers.Properties.VariableNames;
    % Überschriften aktualisieren
    varnames_act{end-1} = 'Rankloss_Platform'; % Matlab-Konform ohne Leerzeichen
    ActTab_i.Properties.VariableNames = varnames_act;

    % Einzelne Spalten mit 0/1 Inhalt für die Aktuierung zusammenfassen.
    % Dazu Spaltenüberschriften interpretieren. Siehe actuation.csv.
    ileg = 0;
    ilegj = 0;
    NLEG = sum(EEFG_update(iFG,:)); % Annahme: keine Antriebsredundanz, etc.
    vargroups_leg_acts = cell(NLEG,1);
    for jj = 1:length(varnames_act)
      [tokens] = regexp(varnames_act{jj},'Aktuierung Bein (\d)','tokens');
      isactcol = false;
      if strcmp(varnames_act{jj}(1:3), 'Var')
        isactcol = true;
        ilegj = ilegj + 1;
      end
      if ~isempty(tokens)
        ileg = str2double(tokens{1}{1});
        ilegj = 1;
      end
      if isactcol || ~isempty(tokens)
        vargroups_leg_acts{ileg} = [vargroups_leg_acts{ileg}, varnames_act(jj)];
      end
    end
    % Füge neue Spalten am Ende hinzu, die die Aktuierung direkt beinhalten
    for ileg = NLEG:-1:1
      legactinfo = [];
      for jj = 1:length(vargroups_leg_acts{ileg})
        % Lese Spalte aus und hänge an
        legactinfo = [legactinfo, ActTab_i.(vargroups_leg_acts{ileg}{jj})]; %#ok<AGROW>
      end
      % Fülle die Matrix mit Nullen auf (damit Tabellen mit
      % unterschiedlichen Beinketten stapelbar sind)
      legactinfo = [legactinfo, zeros(size(legactinfo,1),6-size(legactinfo,2))]; %#ok<AGROW>
      % Entferne Spalten aus ursprünglicher Datenbank
      ActTab_i = removevars(ActTab_i, vargroups_leg_acts{ileg});
      % Füge neue Spalte am Ende ein
      ActTab_i = addvars(ActTab_i, legactinfo, 'After', 1);
      ActTab_i.Properties.VariableNames(2) = {sprintf('Act_Leg%d', ileg)};
    end
    % Ändere den Namen der Spalte für die frei wählbaren Winkel (sonst
    % nicht stapelbar, weil unterschiedliche Spaltennamen
    ActTab_i.Properties.VariableNames(end) = {'Values_Angle_Parameters'};
    % Prüfe, ob die Datenzeile einen passenden Namen hat. Benutze Ausdruck
    % von oben (ohne A-Nummer)
    tokens = regexp(ActTab_i.Name,expression,'tokens');
    I_invalid = cellfun(@isempty,tokens); % Wenn der Name bspw. leer ist, wird der Ausdruck nicht gefunden
    if any(I_invalid)
      warning('%d ungültige Zeile(n) in %s.', sum(I_invalid), acttabfile_csv);
      ActTab_i = ActTab_i(~I_invalid,:);
    end
    % Nachverarbeitung der Spalte für Rangverlust. Wenn dort ein "?" in der
    % csv-Tabelle steht, muss das in NaN übersetzt werden.
    if isa(ActTab_i.Rankloss_Platform, 'cell')
      % Wenn in der csv ein "?" steht, wird die Spalte als cell erzeugt.
      % Wandle die Daten wieder in double um und setze NaN
      ranklosscoll = NaN(size(ActTab_i,1),1);
      I_NaN = strcmp(ActTab_i.Rankloss_Platform,'?');
      ranklosscoll(~I_NaN)=str2double(ActTab_i.Rankloss_Platform(~I_NaN));
      ActTab_i.Rankloss_Platform = ranklosscoll;
    end
    % Stapele die einzelnen Tabellen
    if isempty(ActTab) % Anfang muss nicht bei i=1 sein, wenn Daten fehlen
      ActTab = ActTab_i;
    else
      ActTab = [ActTab;ActTab_i]; %#ok<AGROW>
    end
  end
  if isempty(ActTab)
    warning('Keine einzige Aktuierungstabelle für FG %s', EEstr)
    continue; % Nächster FG
  end
  %% Tabelle nachbearbeiten: Zusätzliche Informationen eintragen
  MaxIdxActJoint = NaN(length(ActTab.Name),1);
  MaxIdxActTechJoint = NaN(length(ActTab.Name),1);
  for k = 1:length(ActTab.Name)
    PName_k = ActTab.Name{k};
    % Finde die aktuierte PKM in der Kinematik-Tabelle
    I_k = find(strcmp(KinTab.Name, PName_k(1:end-2)));
    if isempty(I_k)
      warning('Kein Eintrag in Kinematik-Tabelle für %s. Überspringe', PName_k(1:end-2));
      continue
    end
    assert(length(I_k)==1, sprintf('Kein eindeutiger Eintrag für "%s" in Kinematik-DB', PName_k(1:end-2)));
    SName_TechJoint = KinTab.Beinkette_Tech{I_k};
    % Einzelne technische Gelenke durchgehen und zählen
    TechJointNumbersInChain = zeros(1,6);
    jj = 0;
    for kk = 1:length(SName_TechJoint)
      switch SName_TechJoint(kk)
        case 'R', jv = 1;
        case 'P', jv = 1;
        case 'U', jv = 2;
        case 'S', jv = 2;
        case 'C', jv = 2;
      end
      jj = jj + jv;
      if kk == 1
        TechJointNumbersInChain(1:jj) = kk;
      else
        TechJointNumbersInChain(jj-jv+1:jj) = kk;%TechJointNumbersInChain(jj-jv)+jv;
      end
    end
    % Bestimme maximalen Index des aktuierten Gelenks (unabhängig von
    % technischer Realisierung, nur Gelenk-FG-Nummer)
    MaxIdxActJoint(k) = max(ActTab.Act_Leg1(k,:).*(1:size(ActTab.Act_Leg1,2)));
    % Bestimme Index des technischen Gelenks für das aktuierte Gelenke.
    % Annahme: Symmetrische Aktuierung
    MaxIdxActTechJoint(k) = TechJointNumbersInChain(MaxIdxActJoint(k));
  end
  ActTab = addvars(ActTab, MaxIdxActJoint);
  ActTab.Properties.VariableNames(end) = {'MaxIdxActJoint'};
  ActTab = addvars(ActTab, MaxIdxActTechJoint);
  ActTab.Properties.VariableNames(end) = {'MaxIdxActTechJoint'};
  %% Aktuierungs-Tabelle speichern
  % Prüfe, ob die aktuell ausgelesene CSV-Datei der gleiche Stand ist wie
  % die binär einzulesende Datei
  acttabfile_mat = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_act.mat']);
  write_new_act = false;
  if exist(acttabfile_mat, 'file')
    tmp = load(acttabfile_mat);
    ActTab_old = tmp.ActTab;
    if ~isequaln(ActTab, ActTab_old)
      write_new_act = true; % Dateiinhalt wird sich ändern. Schreibe neu.
    end
  else
    mkdirs(fileparts(acttabfile_mat)); % Ordner könnte nicht existieren
    write_new_act = true; % Schreibe mat-Datei neu, sie existiert noch nicht
  end
  if write_new_act
    save(acttabfile_mat, 'ActTab');
  end
end

%% Schreibe Gesamt-Tabelle neu
kintaballfile_mat = fullfile(repopath, 'parrob_list_kin.mat');
write_new_kin_all = false;
if exist(kintaballfile_mat, 'file') && ~isempty(KinTab_All)
  tmp = load(kintaballfile_mat);
  KinTab_All_old = tmp.KinTab;
  if ~isequaln(KinTab_All, KinTab_All_old)
    write_new_kin_all = true; % Dateiinhalt wird sich ändern. Schreibe neu.
  end
else
  write_new_kin_all = true; % Schreibe mat-Datei neu, sie existiert noch nicht
end
if write_new_kin_all
  KinTab = KinTab_All;
  save(kintaballfile_mat, 'KinTab');
end
