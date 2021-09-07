% Rufe alle Matlab-Funktionen einmal auf und prüfe damit die Korrektheit
% der mex-Funktionen. Automatische Neu-Erzeugung, falls Syntax-Fehler.
% Dieses Skript ist dann nützlich, wenn die Schnittstellen im Robotik-Repo
% geändert wurden und noch mex-Dateien vorliegen, die mit der alten Version
% kompiliert wurden. Es werden auch die zugrunde liegenden m-Dateien geprüft.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 0; ...
            1 1 1 1 1 1];
parroblib_path = fileparts(which('parroblib_path_init.m'));
for i_FG = 1:size(EEFG_Ges,1)
  EE_FG = logical(EEFG_Ges(i_FG,:));
  EE_FG_str = sprintf('%dT%dR', sum(EE_FG(1:3)), sum(EE_FG(4:6)));
  [PNames_Kin, PNames_Act] = parroblib_filter_robots(EE_FG, 6);
  % Prüfe nur eine PKM je Gestell- und Plattform-Nummer (die Template-
  % Funktionen sind für alle G-P-Varianten identisch.
  PNames_noGP = cell(size(PNames_Kin));
  for i = 1:length(PNames_Kin)
    PNames_noGP{i} = PNames_Kin{i}(1:end-4); % Annahme: Nur einstellige Nummern
  end
  [~,III] = unique(PNames_noGP);
  fprintf('Prüfe %d PKM mit %s\n', length(III), EE_FG_str);
  for ii = III'
    I_act = find(contains(PNames_Act, PNames_Kin{ii}),1,'first');
    PName = PNames_Act{I_act};
    % Debug: Einzelnen Roboter prüfen.
%     if ~contains(PName, 'P3PRRRR3'), continue; end
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    % Suche nach den m- und mex-Dateien
    tpl_dir = fullfile(parroblib_path, sprintf('sym_%s', EE_FG_str), ...
      PName_Legs, 'tpl');
    mexfilelist = dir(fullfile(tpl_dir, '*_mex.mex*'));
    mfilelist = dir(fullfile(tpl_dir, '*.m'));
    filelist = [mfilelist;mexfilelist];
    if isempty(filelist)
      % Nichts zu prüfen. Es gibt keine mex-Dateien
      continue
    end
    fprintf('Prüfe PKM %d/%d (%s) (%d mex-Dateien und %d m-Dateien liegen in tpl-Ordner %s)\n', ...
      find(III==ii), length(III), PNames_Kin{ii}, length(mexfilelist), length(mfilelist), tpl_dir)
    % Initialisiere Matlab-Klasse und setze auf Nutzung von M-Funktionen
    RP = parroblib_create_robot_class(PName, 1, 0.3);
    % Gehe alle m- und mex-Dateien durch, die da sind. Fange mit m an.
    % Dadurch werden die Vorlagen-Funktionen meistens schon neu generiert.
    RP.fill_fcn_handles(false, false); RP_mex_status = false;
    for kk = 1:length(filelist)
      if contains(filelist(kk).name, '_mex')
        % Prüfe ab jetzt die mex-Dateien. Die m-Dateien sind fertig.
        if RP_mex_status == false
          RP.fill_fcn_handles(true, false);
          RP_mex_status = true;
        end
      end
      for retryiter = 1:3 % mehrfache Neuversuche zur Fehlerkorrektur
        recompile = false;
        % Einzelne Fälle für die mex-Dateien durchgehen und jeweils
        % Dummy-Aufruf der Funktionen, um Syntax-Fehler aufzudecken.
        if contains(filelist(kk).name, 'invkin_mex') || ...
            contains(filelist(kk).name, 'invkin.m')
          try
            % Gebe mehr als einen Startwert vor (neue Schnittstelle seit 2021-06)
            RP.invkin2(zeros(6,1), rand(RP.NJ,3));
          catch err
            if ~strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
              recompile = true;
            end
          end
        end
        if contains(filelist(kk).name, 'invkin3')
          try
            % Prüfe Dateiinhalte auf charakteristische Einträge
            if ~RP_mex_status % Nicht für mex-Dateien
              filetext = fileread(fullfile(tpl_dir, filelist(kk).name));
              if ~contains(filetext, 'installspace_thresh')
                error('Textfragment "installspace_thresh" nicht gefunden. Alte Version.');
              end
            end
            % Prüfe, ob Korrektur von Fehler bei Kollisionsprüfung da ist
            % Behoben ca. 2021-07; max/min mit Eingabe variabler Länge
            s = struct('avoid_collision_finish', true);
            [~,~,~,Stats] = RP.invkin4(zeros(6,1), rand(RP.NJ,3), s);
            % Gebe mehr als einen Startwert vor (neue Schnittstelle seit 2021-06)
            [~,~,~,Stats] = RP.invkin4(zeros(6,1), rand(RP.NJ,3));
            % Prüfe, ob neue Ausgabe (seit 2021-06) da ist.
            tmp = Stats.coll;
          catch err
            if ~strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
              recompile = true;
            end
          end
        end
        if contains(filelist(kk).name, 'invkin_traj')
          try
            % Prüfe Dateiinhalte auf charakteristische Einträge
            if ~RP_mex_status % Nicht für mex-Dateien
              filetext = fileread(fullfile(tpl_dir, filelist(kk).name));
              if ~contains(filetext, 'installspace_thresh')
                error('Textfragment "installspace_thresh" nicht gefunden. Alte Version.');
              end
            end
            % Führe die Funktion aus
            % Prüfe Trajektorie mit nur einem Punkt (Bug, der am 16.08.2021
            % behoben wurde). Prüfung zuerst, damit Syntax-Fehler vor Inf/NaN-Fehler kommt.
            RP.invkin2_traj(zeros(1,6), zeros(1,6), zeros(1,6), 0, zeros(RP.NJ,1));
            % Prüfe mit Dummy-Trajektorie (aus zwei Punkten)
            RP.invkin2_traj(zeros(2,6), zeros(2,6), zeros(2,6), [0;1], zeros(RP.NJ,1));
          catch err
            if ~strcmp(err.identifier, 'MATLAB:svd:matrixWithNaNInf')
              recompile = true;
            end
          end
        end
        % Falls ein Fehler vorliegt, wird neu kompiliert oder generiert.
        if recompile
          warning('Fehler beim Aufruf von Funktion %s. Fehler: %s', ...
            filelist(kk).name, err.message);
          if retryiter == 2
            warning('Zweiter Versuch ohne Erfolg. Voraussichtlich sind die Vorlagen-Funktionen veraltet');
            serroblib_create_template_functions({LEG_Names{1}},  false,false); %#ok<CCAT1>
            parroblib_create_template_functions({PNames_Kin{ii}},false,false); %#ok<CCAT1>
          elseif retryiter == 3
            error('Auch beim dritten Versuch kein Erfolg');
          end
          if contains(filelist(kk).name, '_mex')
            % Mex-Dateien werden neu kompiliert. m-Dateien nicht.
            [~,mexbasename] = fileparts(filelist(kk).name);
            matlabfcn2mex({mexbasename(1:end-4)}); % ohne "_mex"
          end
        else
          % Alles funktioniert. Keine Neuversuche notwendig.
          if retryiter > 1
            fprintf('Datei %s wurde korrigiert\n', filelist(kk).name);
          end
          break
        end % recompile
      end % retryiter
    end % kk (Mex-Dateien)
  end % ii (PKM)
end % i_FG (EE-FG)
    
