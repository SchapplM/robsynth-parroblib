% Rufe alle Matlab-Funktionen einmal auf und prüfe damit die Korrektheit
% der mex-Funktionen. Automatische Neu-Erzeugung, falls Syntax-Fehler.
% Dieses Skript ist dann nützlich, wenn die Schnittstellen im Robotik-Repo
% geändert wurden und noch mex-Dateien vorliegen, die mit der alten Version
% kompiliert wurden.

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
  fprintf('Prüfe %d PKM mit %s\n', length(PNames_Kin), EE_FG_str);
  for ii = 1:length(PNames_Kin)
    I_act = find(contains(PNames_Act, PNames_Kin{ii}),1,'first');
    PName = PNames_Act{I_act};
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    % Suche nach den mex-Dateien
    tpl_dir = fullfile(parroblib_path, sprintf('sym_%s', EE_FG_str), ...
      PName_Legs, 'tpl');
    mexfilelist = dir(fullfile(tpl_dir, '*_mex.mex*'));
    if isempty(mexfilelist)
      % Nichts zu prüfen. Es gibt keine mex-Dateien
      continue
    end
    fprintf('Prüfe PKM %d/%d (%s) (%d mex-Dateien liegen in tpl-Ordner)\n', ...
      ii, length(PNames_Kin), PNames_Kin{ii}, length(mexfilelist))
    % Initialisiere Matlab-Klasse und setze auf Nutzung von Mex-Funktionen
    RP = parroblib_create_robot_class(PName, 1, 0.3);
    list_missing_mex = RP.fill_fcn_handles(true, false);
    % Gehe alle mex-Dateien durch, die da sind.
    for kk = 1:length(mexfilelist)
      for retryiter = 1:3 % mehrfache Neuversuche zur Fehlerkorrektur
        recompile = false;
        % Einzelne Fälle für die mex-Dateien durchgehen und jeweils
        % Dummy-Aufruf der Funktionen, um Syntax-Fehler aufzudecken.
        if contains(mexfilelist(kk).name, 'invkin_mex')
          try
            RP.invkin2(NaN(6,1), zeros(RP.NJ,1));
          catch err
            recompile = true;
          end
        end
        if contains(mexfilelist(kk).name, 'invkin3_mex')
          try
            RP.invkin4(NaN(6,1), zeros(RP.NJ,1));
          catch err
            recompile = true;
          end
        end
        if contains(mexfilelist(kk).name, 'invkin_traj_mex')
          try
            RP.invkin2_traj(NaN(2,6), NaN(2,6), NaN(2,6), [0;1], zeros(RP.NJ,1));
          catch err
            recompile = true;
          end
        end
        % Falls ein Fehler vorliegt, wird neu kompiliert oder generiert.
        if recompile
          warning('Fehler beim Aufruf von Funktion %s. Fehler: %s', ...
            mexfilelist(kk).name, err.message);
          if retryiter == 2
            warning('Zweiter Versuch ohne Erfolg. Voraussichtlich sind die Vorlagen-Funktionen veraltet');
            serroblib_create_template_functions({LEG_Names{1}},  false,false); %#ok<CCAT1>
            parroblib_create_template_functions({PNames_Kin{ii}},false,false); %#ok<CCAT1>
          elseif retryiter == 3
            error('Auch beim dritten Versuch kein Erfolg');
          end
          [~,mexbasename] = fileparts(mexfilelist(kk).name);
          matlabfcn2mex({mexbasename(1:end-4)}); % ohne "_mex"
        else
          % Alles funktioniert. Keine Neuversuche notwendig.
          if retryiter > 1
            fprintf('Mex-Datei %s wurde korrigiert\n', mexfilelist(kk).name);
          end
          break
        end % recompile
      end % retryiter
    end % kk (Mex-Dateien)
  end % ii (PKM)
end % i_FG (EE-FG)
    