% Generiere .mat-Dateien zum schnelleren Durchsuchen der Datenbank
% 
% Eingabe:
% N_update (optional)
%   Anzahl der Beinketten der Roboter, für die die .mat-Datei
%   aktualisiert werden soll.
% 
% Schreibt Dateien:
% symxleg_list_kin.mat (x=Anzahl Beinketten aus Eingabe). Enthält Variable
% KinTab. Entspricht symxleg_list.csv. Matlab-Tabelle mit Feldern:
%   Name
%     Name der PKM-Kinematik (mit Koppelgelenk): "P6PRRRRR1G1P5"
%   Beinkette
%     Name der seriellen Beinkette (siehe serroblib): "S6PRRRRR1"
%   EEFG [1x6]
%     Einträge mit 1/0, je nachdem welcher EE-FG bewegt werden kann. FG
%     werden im Basis-KS gezählt (Plattform-Geschw. und -winkelgeschw.)
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

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-11
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_gen_bitarrays(N_update)

if nargin == 0
  N_update = 3:6; % Aktualisiere alle Roboter
end
repopath=fileparts(which('parroblib_path_init.m'));
for NLEG = N_update(:)'
  %% Kinematik-Tabelle durchsuchen
  kintabfile_csv = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list.csv', NLEG));
  if ~exist(kintabfile_csv, 'file')
    warning('Datei %s existiert nicht. Sollte aber im Repo enthalten sein.', kintabfile_csv);
    return
  end
  % Tabelle lesen (und dabei Überschriften nachbessern)
  KinTab = readtable(kintabfile_csv, 'HeaderLines', 2);
  KinTab_headers = readtable(kintabfile_csv, 'PreserveVariableNames', true);
  varnames_kin = KinTab_headers.Properties.VariableNames;
  varnames_kin{1} = 'Name'; % doppelte Überschriftszeile funktioniert ...
  varnames_kin{2} = 'Beinkette'; % ... nicht gut mit readtable
  KinTab.Properties.VariableNames = varnames_kin;
  % Fasse die EE-FG-Spalten zusammen
  EEFG = NaN(size(KinTab,1), 6);
  for i = 1:6
    EEFG(:,i) = KinTab.(varnames_kin{end-6+i});
  end
  % Lösche die Überflüssigen Spalten und erstelle eine neue
  KinTab = removevars(KinTab, varnames_kin(end-5:end));
  KinTab = addvars(KinTab, EEFG);
  KinTab.Properties.VariableNames(end) = {'EEFG'};
  
  %% Kinematik-Tabelle speichern
  % Prüfe, ob die aktuell ausgelesene CSV-Datei der gleiche Stand ist wie
  % die binär einzulesende Datei
  kintabfile_mat = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list_kin.mat', NLEG));
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
  %% Aktuierungs-Tabelle auslesen
  for i = 1:size(KinTab,1)
    PName_Kin = KinTab.Name{i};
    % Extrahiere Beinketten-Name der PKM (ohne Ausrichtungs-Nummern G/P)
    % Unter diesem Namen ist die Aktuierungstabelle gespeichert
    expression = 'P(\d)([RP]+)(\d+)[V]?(\d*)[G]?(\d*)[P]?(\d*)'; % Format "P3RRRG1P11A1" oder "P3RRR1V1G1P1A1"
    [tokens, ~] = regexp(PName_Kin,expression,'tokens','match');
    if isempty(tokens)
      error('Eintrag %s aus der Tabelle stimmt nicht mit Namensschema überein.', PName_Kin);
    end
    res = tokens{1};
    if isempty(res{4}) % serielle Kette ist keine abgeleitete Variante
      PName_Legs = ['P', res{1}, res{2}, res{3}];
    else % serielle Kette ist eine Variante abgeleitet aus Hauptmodell
      PName_Legs = ['P', res{1}, res{2}, res{3}, 'V', res{4}];
    end
    acttabfile_csv = fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Legs, 'actuation.csv');
    if ~exist(acttabfile_csv, 'file')
      % Diese Aktuierung zu der Kinematik ist noch nicht gespeichert / generiert.
      warning('Zu Kinematik %s gibt es keine Aktuierungstabelle %s', PName_Kin, acttabfile_csv);
      return
    end
    ActTab_i = readtable(acttabfile_csv, 'HeaderLines', 1);
    ActTab_i_headers = readtable(acttabfile_csv, 'PreserveVariableNames', ...
      true, 'ReadVariableNames', true);
    varnames_act = ActTab_i_headers.Properties.VariableNames;
    % Überschriften aktualisieren
    varnames_act{end-1} = 'Rankloss_Platform'; % Matlab-Konform ohne Leerzeichen
    ActTab_i.Properties.VariableNames = varnames_act;

    % Einzelne Spalten mit 0/1 Inhalt für die Aktuierung zusammenfassen.
    % Dazu Spaltenüberschriften interpretieren. Siehe actuation.csv.
    ileg = 0;
    ilegj = 0;
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
    % Stapele die einzelnen Tabellen
    if i == 1
      ActTab = ActTab_i;
    else
      ActTab = [ActTab;ActTab_i]; %#ok<AGROW>
    end
  end
  %% Aktuierungs-Tabelle speichern
  % Prüfe, ob die aktuell ausgelesene CSV-Datei der gleiche Stand ist wie
  % die binär einzulesende Datei
  acttabfile_mat = fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list_act.mat', NLEG));
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
