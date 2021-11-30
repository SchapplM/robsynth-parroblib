% Prüfe die Template-Funktionen für einen bestimmten Roboter.
% Aktualisiere die Funktionen, falls sie veraltet sind.
% 
% Eingabe:
% Names {1 x n} cell array
%   Namen der Roboter, für die die Funktionen aktualisiert werden.
%   Falls leer gelassen: Alle.
% verbosity
%   Grad der Textausgabe (0=aus, 1=Fortschritt)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_update_template_functions(Names, verbosity)

%% Prüfe Eingabe
if nargin < 1 || isempty(Names) % keine Eingabe einer Liste. Nehme alle.
  % Stelle Liste aller Roboter zusammen
  Names = {};
  EEFG_Ges = logical(...
    [1 1 0 0 0 0; 1 1 0 0 0 1; 1 1 1 0 0 0;  1 1 1 0 0 1; ...
     1 1 1 1 1 0; 1 1 1 1 1 1]);
  for j = 1:size(EEFG_Ges,1)
    [PNames_Kin, ~, ~] = parroblib_filter_robots(EEFG_Ges(j,:), 6); 
    Names = [Names, PNames_Kin]; %#ok<AGROW>
  end
end
if nargin < 2
  verbosity = 0;
end
%% Bestimme aktuelle Version der jeweiligen Vorlagen-Funktion
filelist = {'pkm_invkin.m.template', 'pkm_invkin3.m.template', ...
  'pkm_invkin_traj.m.template'};
fileversions = struct('pkm_invkin', 0, 'pkm_invkin3', 0, 'pkm_invkin_traj', 0);
kinematics_dir = fullfile(fileparts(which('robotics_toolbox_path_init.m')), ...
  'kinematics');
for i = 1:length(filelist)
  [tokens,~] = regexp(fileread(fullfile(kinematics_dir, filelist{i})), ...
    '''version'', (\d+)', 'tokens', 'match');
  key = filelist{i}(1:end-11); % Entferne Endung .m.template
  fileversions.(key) = str2double(tokens{1}{1});
end

%% Gehe alle Dateien durch und prüfe die Version
parroblib_path = fileparts(which('parroblib_path_init.m'));
% Erstelle Namen der Roboter
PNames_noGP = cell(size(Names));
for i = 1:length(Names)
  [~, ~, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(Names{i}, 0);
  PNames_noGP{i} = PName_Legs; % Annahme: Nur einstellige Nummern
end
[~,III] = unique(PNames_noGP);
for ii = III'
  [~, LEG_Names, ~, ~, ~, ~, EE_FG, PName_Kin, PName_Legs, ~] = parroblib_load_robot(Names{ii}, 1);
  PName = [PName_Kin, 'A0'];
  % Suche nach den m- und mex-Dateien
  EE_FG_str = sprintf('%dT%dR', sum(EE_FG(1:3)), sum(EE_FG(4:6)));
  tpl_dir = fullfile(parroblib_path, sprintf('sym_%s', EE_FG_str), ...
    PName_Legs, 'tpl');
  mexfilelist = dir(fullfile(tpl_dir, '*_mex.mex*'));
  mfilelist = dir(fullfile(tpl_dir, '*.m'));
  filelist = [mfilelist;mexfilelist];
  if isempty(filelist)
    % Nichts zu prüfen. Es gibt keine mex-Dateien
    continue
  end
  if verbosity
    fprintf('Prüfe PKM %d/%d (%s) (%d mex-Dateien und %d m-Dateien liegen in tpl-Ordner %s)\n', ...
      find(III==ii), length(III), PName_Kin, length(mexfilelist), length(mfilelist), tpl_dir)
  end
  % Initialisiere Matlab-Klasse und setze auf Nutzung von M-Funktionen
  RP = parroblib_create_robot_class(PName, 1, 0.3);
  % Zufallswerte für Parameter setzen. Sonst kommen NaN-Fehler und
  % überdecken die Syntax-Fehler
  for i = 1:RP.NLEG % Test-Werte sind nicht symmetrisch. Ist aber egal.
    RP.Leg(i).gen_testsettings(true, true);
    RP.Leg(i).qlim = repmat([-1,1],RP.Leg(i).NQJ,1); % Damit die IK-Neuversuche nicht davon abhängen
  end
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
    % Mehrfache Neuversuche zur Fehlerkorrektur. Manchmal wird in v1 noch
    % die mex-Datei aus einer falschen Vorlage erstellt. Daher 3 Versuche
    % und nicht 2.
    for retryiter = 1:3
      recompile = false;
      % Einzelne Fälle für die mex-Dateien durchgehen und jeweils
      % Dummy-Aufruf der Funktionen, um Syntax-Fehler aufzudecken.
      if contains(filelist(kk).name, 'invkin_mex') || ...
          contains(filelist(kk).name, 'invkin.m')
        try
          s = struct('n_max', 1, 'retry_limit', -1); % keine Ausführung
          [~,~,~,Stats]=RP.invkin2(zeros(6,1), rand(RP.NJ,3), s);
          if Stats.version < fileversions.pkm_invkin % hier wird die aktuelle Version eingetragen
            error('Version der Datei ist zu alt (%d). Aktuell %d', ...
              Stats.version, fileversions.pkm_invkin);
          end
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin3')
        try
          s = struct('n_max', 1, 'retry_limit', -1); % Keine Durchführung des Algorithmus
          [~,~,~,Stats] = RP.invkin4(zeros(6,1), rand(RP.NJ,3), s);
          if Stats.version < fileversions.pkm_invkin3 % hier wird die aktuelle Version eingetragen
            error('Version der Datei ist zu alt (%d). Aktuell: %d', ...
              Stats.version, fileversions.pkm_invkin3);
          end
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin_traj')
        try
          [~, ~, ~, ~, ~, ~, ~, Stats] = RP.invkin2_traj(zeros(1,6), zeros(1,6), zeros(1,6), 0, zeros(RP.NJ,1));
          if Stats.version < fileversions.pkm_invkin_traj % hier wird die aktuelle Version eingetragen
            error('Version der Datei ist zu alt (%d). Aktuell: %d', ...
              Stats.version, fileversions.pkm_invkin_traj);
          end
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
        if verbosity
          fprintf('Fehler beim Aufruf von Funktion %s. Fehler: %s\n', ...
            filelist(kk).name, err.message);
        end
        if retryiter == 2
          if verbosity
            fprintf(['Zweiter Versuch ohne Erfolg. Voraussichtlich sind ', ...
              'die Vorlagen-Funktionen veraltet\n']);
          end
          serroblib_create_template_functions({LEG_Names{1}},  false,false); %#ok<CCAT1>
          parroblib_create_template_functions({PName_Kin},false,false);
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
        if retryiter > 1 && verbosity
          fprintf('Datei %s wurde korrigiert\n', filelist(kk).name);
        end
        break
      end % recompile
    end % retryiter
  end % kk (Mex-Dateien)
end % ii (PKM)