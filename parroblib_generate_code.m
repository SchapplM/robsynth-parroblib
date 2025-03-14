% Generiere Matlab-Code mit Maple-Dynamik-Toolbox für PKM-Roboterstrukturen
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen
% force_par
%   Erzwinge Neu-Kompilierung, egal ob bereits generierter Code (für die
%   PKM) vorliegt 
% force_ser
%   Erzwinge Neu-Kompilierung der seriellen Kette
% 
% Vorher: 
% * Funktion maplerepo_path.m muss vorliegen mit Rückgabe des
%   Repo-Pfades der Maple-Dynamik-Toolbox ("HybrDyn")
% * Maple-Eingabedaten müssen für die Roboterstruktur mit
%   parroblib_generate_mapleinput.m erzeugt werden
% 
% Siehe auch: serroblib_generate_code.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_generate_code(Names, force_par, force_ser)
%% Init
if nargin < 2
  force_par = false;
end
if nargin < 3
  force_ser = false;
end

repopath=fileparts(which('parroblib_path_init.m'));
serrobpath=fileparts(which('serroblib_path_init.m'));
%% Roboterstrukturen durchgehen
for i = 1:length(Names)
  
  %% Vorbereitung und Prüfung
  n = Names{i};

  % Daten über den Roboter zusammenstellen
  [NLEG, LEG_Names, Actuation, Coupling, ActNr, ~, EE_FG0, PName_Kin, PName_Legs, AdditionalInfo_Akt] = parroblib_load_robot(n);
  EEstr = sprintf('%dT%dR', sum(EE_FG0(1:3)), sum(EE_FG0(4:6)));
  % Robotereigenschaften aus dem Namen auslesen.
  % TODO: Einbindung nicht-symmetrischer PKM
  
  % Prüfen, ob der Roboter modelliert werden kann
  for kk = 1:length(Actuation)
    if any(Actuation{kk} == NLEG)
      warning('Der Roboter %s kann nicht modelliert werden. Aktuierte Plattform-Koppelgelenke werden noch nicht unterstützt', n);
      continue
    end
  end
  % Weitere Prüfungen (siehe generate_mapleinput)
  % Robotereigenschaften der Beinketten prüfen. Siehe serroblib_gen_bitarrays
  legdata = load(fullfile(serrobpath, sprintf('mdl_%sdof', LEG_Names{1}(2)), ...
    sprintf('S%s_list', LEG_Names{1}(2))));
  AddInfo_Leg = legdata.AdditionalInfo(strcmp(legdata.Names_Ndof,LEG_Names{1}),:);
  if all(EE_FG0==[1 1 0 0 0 1]) && AddInfo_Leg(1) > 2 || AddInfo_Leg(1) > 3
    warning('Symbolischer Code kann nicht basierend auf %s-Beinkette gebildet werden. Zu viele positionsbeeinflussende Gelenke (%d)', LEG_Names{1}, AddInfo_Leg(1));
    continue
  end
  
  % Pfad zur Maple-Dynamik-Toolbox (muss im Repo abgelegt werden; s.o.)
  mrp = fileparts(which('hybrdyn_path_init.m'));
  if ispc()
    [~, mrp_wsl] = system(sprintf('wsl wslpath -u "%s"', mrp));
      mrp_wsl = strtrim(mrp_wsl);
  else
    mrp_wsl = mrp;
  end
  
  % Maple-Toolbox-Eingabe laden (wurde an anderer Stelle erzeugt)
  % (durch parroblib_generate_mapleinput.m)
  mapleinputfile=fullfile(repopath, ['sym_', EEstr], PName_Legs, ...
    sprintf('hd_G%dP%dA%d',Coupling(1),Coupling(2),ActNr), sprintf('robot_env_par_%s', n));
  if ~exist(mapleinputfile, 'file')
    warning('Datei %s existiert nicht. Wurde `parroblib_generate_mapleinput.m` ausgeführt?', fileparts(mapleinputfile) );
    continue
  end
  % Verzeichnisse für die zu erzeugenden Matlab-Funktionen
  outputdir_tb_par = fullfile(mrp, 'codeexport', n, 'matlabfcn'); % Verzeichnis in der Maple-Toolbox
  outputdir_local = fileparts(mapleinputfile); % Verzeichnis in der Bibliothek
  
  % Prüfe, ob Code schon einmal generiert wurde 
  % (und im Zielverzeichnis vorliegt)
  if ~force_par && length(dir(fullfile(outputdir_local, '*.m'))) > 1
    % Annahme: Wenn bereits Code erstellt wurde, ist dieser vollständig.
    % Für A1... Modelle gibt es nur zwei Matlab-Funktionen
    fprintf('Code existiert bereits in %s\n', outputdir_local);
    continue
  end
  
  %% Code-Generierung serielle Beinkette
  % Code-Erstellung für serielle Beinkette starten, falls diese nicht
  % vorliegt
  LEG_Names_unique = unique(LEG_Names);
  for k = 1:length(LEG_Names_unique)
    % Suche nach temporären Dateien im Arbeitsverzeichnis von HybrDyn
    tmpdir_tb_ser = fullfile(mrp, 'codeexport', LEG_Names_unique{k}, 'tmp'); % tmp-Verzeichnis in der Maple-Toolbox
    if ~force_ser && length(dir(fullfile(tmpdir_tb_ser, '*_maple.m'))) > 30
      % Die serielle Beinkette wurde wahrscheinlich schon generiert (auch
      % mit Dynamik)
      continue
    end
    % Generiere diese Beinkette neu (ohne Rückkopieren der
    % Matlab-Funktionen ins SerRobLib-Repo)
    serroblib_generate_mapleinput(LEG_Names_unique(k));
    serroblib_generate_code(LEG_Names_unique(k), true, true, 3);
  end

  %% Code-Generierung für allgemeine PKM dieses Kinematik-Typs
  % Code-Erstellung für Dynamik nur einmal durchführen, da die Dynamik in
  % Plattformkoordinaten unabhängig von der Aktuierung ist. Gebe
  % gleichwertige Gestell-Methoden an (mit gleichem symbolischen Code)
  switch Coupling(1) % siehe align_base_coupling.m
    case 1
      basecoupling_equiv = [5,4,8];
    case 2
      basecoupling_equiv = [6,4,8];
    case 3
      basecoupling_equiv = [7,4,8];
    case 4
      basecoupling_equiv = 8;
    case 5
      basecoupling_equiv = [1 4 8];
    case 6
      basecoupling_equiv = [2 4 8];
    case 7
      basecoupling_equiv = [3 4 8];
    case 8
      basecoupling_equiv = 4;
  end
  % Der Code für diese Gestell-Kinematik wird zuletzt hinzugefügt (damit
  % wird existierender Code mit einer anderen Nummer bevorzugt)
  basecoupling_equiv = [basecoupling_equiv, Coupling(1)]; %#ok<AGROW>
  G_num = Coupling(1);
  % Suche in den anderen Code-Ordnern nach bestehenden Robotern
  if ~force_par
    for k = basecoupling_equiv
      outputdir_local_A0_k = fullfile(repopath, ['sym_', EEstr], ...
        PName_Legs, sprintf('hd_G%dA0', k));
      if length(dir(fullfile(outputdir_local_A0_k, '*.m'))) > 10
        G_num = k;
        fprintf(['Bereits passendes allgemeines Plattform-Dynamikmodell ', ...
          'für Gestell %d gefunden. Benutze dieses. Keine Neu-Generierung.\n'], k);
        break;
      end
    end
  end
  n_A0 = [PName_Kin(1:end-4),sprintf('G%d',G_num),'A0']; % Entferne P-Nummer wieder. Ist für A0-Modell egal.
  % Roboternamen, Datei- und Ordnernamen für allgemeinen Fall definieren
  outputdir_tb_par_A0 = fullfile(mrp, 'codeexport', n_A0, 'matlabfcn'); % Verzeichnis in der Maple-Toolbox
  mapleinputfile_A0=fullfile(repopath, ['sym_', EEstr], PName_Legs, ...
    sprintf('hd_G%dA0', G_num), sprintf('robot_env_par_%s', n_A0));
  outputdir_local_A0 = fileparts(mapleinputfile_A0);
  % Hier kann noch nicht robust verhindert werden, dass für mehrere
  % Aktuierungen die Dynamik mehrfach für A0 generiert wird (bei force_par)
  if  ... % Nur Generierung, wenn noch kein Code vorhanden oder erzwungen:
      (force_par || ~force_par && length(dir(fullfile(outputdir_local_A0, '*.m'))) < 10) % im A0-Ordner sind mehr Dateien
    % Namen des Roboters in A0-Definition nachbearbeiten
    mkdirs(fileparts(mapleinputfile_A0));
    copyfile(mapleinputfile, mapleinputfile_A0);
    % Allgemeiner Fall heißt "A0", ursprünglich eingegebener Fall heißt "A1"
    if ispc() % Führe den Befehl in WSL aus (mit angepasstem Pfad)
      [~, mapleinputfile_A0_wsl] = system(sprintf('wsl wslpath -u "%s"', mapleinputfile_A0));
      mapleinputfile_A0_wsl = strtrim(mapleinputfile_A0_wsl);
      WSL_prefix = 'wsl ';
    else % Normale Unix-Shell
      WSL_prefix = '';
      mapleinputfile_A0_wsl = mapleinputfile_A0;
    end
    system_wsl(sprintf('sed -i "s/%s/%s/g" "%s"', n, n_A0, mapleinputfile_A0_wsl));

    % aktiviere den Export der Dynamik-Funktionen. Die Definitionsdatei
    % überschreibt die Einstellungen in den Arbeitsblättern, da sie später
    % geladen wird. Daher darf hier das false nur auskommentiert werden.
    % Würde "true" stattdessen gesetzt, würde parallel mehrfach das gleiche
    % gemacht werden.
    system_wsl(sprintf('sed -i "s/codeexport_invdyn := false:/# codeexport_invdyn := false:/g" %s', ...
      mapleinputfile_A0_wsl ));
    system_wsl(sprintf('sed -i "s/codeexport_corvec := false:/# codeexport_corvec := false:/g" %s', ...
      mapleinputfile_A0_wsl ));
    system_wsl(sprintf('sed -i "s/codeexport_grav := false:/# codeexport_grav := false:/g" %s', ...
      mapleinputfile_A0_wsl ));
    system_wsl(sprintf('sed -i "s/codeexport_inertia := false:/# codeexport_inertia := false:/g" %s', ...
      mapleinputfile_A0_wsl ));
    
    % Reduziere den Optimierungsgrad bei der Code-Generierung, da es sonst
    % für 6FG-Systeme zu lange dauert
    if NLEG == 6
      system_wsl(sprintf('sed -i "s/codegen_opt := 2:/codegen_opt := 1:/g" %s', ...
        mapleinputfile_A0_wsl ));
    end

    % Definition kopieren und Code-Erstellung starten. Alle Tests für PKM
    % dort durchführen
    copyfile(mapleinputfile_A0, fullfile(mrp, 'robot_codegen_definitions', 'robot_env_par') );
    fprintf('Starte Dynamik-Code-Generierung %d/%d für %s\n', i, length(Names), n_A0);
    system_wsl( sprintf('cd %s && ./robot_codegen_start.sh --fixb_only --parrob --parallel --not_gen_serial', mrp_wsl) );
    
    % generierten Code zurückkopieren (alle .m-Dateien)
    for f = dir(fullfile(outputdir_tb_par_A0, '*.m'))'
      copyfile(fullfile(outputdir_tb_par_A0, f.name), fullfile(outputdir_local_A0, f.name));
    end
  end

  %% Code-Generierung dieser PKM (Berücksichtigung der Aktuierung)
  % Robotereigenschaft der PKM prüfen.
  if AdditionalInfo_Akt(1) > 0 % Rangverlust
    warning('PKM hat laut Datenbank/Struktursynthese Rangverlust. Symbolischer Code für aktuierte PKM nicht sinnvoll generierbar.');
    % Abbruch wegen der Prüfung erst hier. Dadurch werden PKM ohne
    % Aktuierung noch generiert (zum Testen der Dynamik)
    continue
  end
  
  % Eingabedatei für parallelen Roboter kopieren
  copyfile( mapleinputfile, fullfile(mrp, 'robot_codegen_definitions', 'robot_env_par') );
  
  % Code-Erstellung für parallelen Roboter starten (ohne Generierung der
  % Beinketten; das wurde oben schon gemacht). Daher auch Tests
  % deaktivieren (die Dynamik wird für diesen Roboter nicht generiert)
  fprintf('Starte Kinematik Code-Generierung %d/%d für %s\n', i, length(Names), n);
  system_wsl( sprintf(['cd %s && ./robot_codegen_start.sh --fixb_only ' ...
    '--parrob --not_gen_serial --notest --kinematics_only'], mrp_wsl) );
  
  % generierten Code zurückkopieren (alle .m-Dateien)
  for f = dir(fullfile(outputdir_tb_par, '*.m'))'
    copyfile(fullfile(outputdir_tb_par, f.name), fullfile(outputdir_local, f.name));
  end
  %% Zurückkopierten Code nachverarbeiten
  % Die Dynamik wird symbolisch im "A0"-Modell gespeichert.
  % Für die verschiedenen Aktuierungen wird die Dynamik dann nur noch
  % aufgerufen. Ändere die Funktionsaufrufe in den Funktionsdateien
  for f = dir(fullfile(outputdir_local, '*.m'))'
    if ispc()
      [~, filename_wsl] = system(sprintf('wsl wslpath -u "%s"', ...
        fullfile(outputdir_local, f.name)));
        filename_wsl = strtrim(filename_wsl);
    else
      filename_wsl = fullfile(outputdir_local, f.name);
    end
    % Funktionsaufruf der inversen Dynamik auf A0 beziehen
    system_wsl(sprintf('sed -i "s/%s_invdyn_para_pf/%s_invdyn_para_pf/g" %s', ...
      n, n_A0, filename_wsl ));
  end
  % Dateien löschen, die für alle Aktuierungsvarianten gleich sind, aber
  % trotzdem doppelt erzeugt werden könnten.
  % (wird aktuell nicht doppelt erzeugt wegen "--kinematics_only")
  minparfile = fullfile(outputdir_local, [n, '_minimal_parameter_para.m']);
  if exist(minparfile, 'file')
    delete(minparfile);
  end
end
end

function system_wsl(cmd)
  % Betriebssystem-unabhängiger Aufruf von Befehlen: Windows-Linux-Subsystem
  % Ermöglicht den Systemaufruf aus Windows und Linux
  % Vorher: Wechseln in Verzeichnis aus Matlab heraus
  %
  % Eingabe:
  % cmd:
  %  Befehl, der in der Konsole (WSL oder Unix-Terminal) werden soll.
  %  Im Befehl enthaltene Pfade müssen unter Windows ins WSL-Schema
  %  umgewandelt werden
  if ispc() % Führe den Befehl in WSL aus (mit angepasstem Pfad)
    WSL_prefix = 'wsl ';
    cmd = strrep(cmd, '&&', ';'); % && wird von Windows-Shell gefangen
  else % Normale Unix-Shell
    WSL_prefix = '';
  end
  system([WSL_prefix, cmd]);
end
