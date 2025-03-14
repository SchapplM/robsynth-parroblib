% Prüfe die Template-Funktionen für einen bestimmten Roboter.
% Aktualisiere die Funktionen, falls sie veraltet sind.
% 
% Eingabe:
% Names {1 x n} cell array
%   Namen der Roboter, für die die Funktionen aktualisiert werden.
%   Falls leer gelassen: Alle.
% verbosity
%   Grad der Textausgabe (0=aus, 1=Fortschritt)
% ignore_mex
%   Bei true werden die mex-Dateien ignoriert und nur M-Funktionen
%   aktualisiert. Standard: false

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_update_template_functions(Names, verbosity, ignore_mex)
orig_state = warning('off', 'all'); % Warnungen temporär unterdrücken
assert(isa(Names, 'cell'), 'Variable "Names" muss cell-Array sein');
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
if nargin < 3
  ignore_mex = false;
end
%% Bestimme aktuelle Version der jeweiligen Vorlagen-Funktion
filelist = {'pkm_invkin.m.template', 'pkm_invkin3.m.template', ...
  'pkm_invkin_traj.m.template'};
fileversions = struct('fkine_coll', 0, 'pkm_invkin', 0, 'pkm_invkin3', 0, ...
  'pkm_invkin_traj', 0);
files_found = fileversions; % Merke, welche Dateien da sind
kinematics_dir = fullfile(fileparts(which('robotics_toolbox_path_init.m')), ...
  'kinematics');
for i = 1:length(filelist)
  [tokens,~] = regexp(fileread(fullfile(kinematics_dir, filelist{i})), ...
    '''version'', (\d+)', 'tokens', 'match');
  key = filelist{i}(1:end-11); % Entferne Endung .m.template
  if isempty(tokens) % Tritt evtl. bei Dateisystemfehler auf
    warning('Versionsnummer nicht in Datei %s gefunden. Unerwartet.');
    fileversions.(key) = 0; % Dadurch wird eine Neugenerierung erzwungen.
  else
    fileversions.(key) = str2double(tokens{1}{1});
  end
  
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
  if isempty(mfilelist) && isempty(mexfilelist)
    % Nichts zu prüfen. Es gibt keine mex-Dateien
    continue
  end
  % Füge Dummy-Einträge hinzu für invkin_ser-Funktion
  mexfilelist(length(mexfilelist)+1,1).name = 'invkin_ser_mex';
  mfilelist(length(mfilelist)+1,1).name = 'invkin_ser';
  if ignore_mex
    filelist = mfilelist;
  else
    filelist = [mfilelist;mexfilelist];
  end
  if verbosity
    fprintf('Prüfe PKM %d/%d (%s) (%d mex-Dateien und %d m-Dateien liegen in tpl-Ordner %s)\n', ...
      find(III==ii), length(III), PName_Kin, length(mexfilelist), length(mfilelist), tpl_dir)
  end
  % Initialisiere Matlab-Klasse und setze auf Nutzung von M-Funktionen
  RP = parroblib_create_robot_class(PName, '', 1, 0.3);
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
      if contains(filelist(kk).name, 'fkineEE_traj_mex')
        try
          RP.fkineEE2_traj(zeros(1,RP.NJ), zeros(1,RP.NJ), zeros(1,RP.NJ));
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin_mex') || ...
          contains(filelist(kk).name, 'invkin.m')
        if ~RP_mex_status % markiere die Datei als vorhanden (ohne Mex)
          files_found.pkm_invkin = true;
        end
        try
          s = struct('n_max', 1, 'retry_limit', -1); % keine Ausführung
          q0 = rand(RP.NJ,3);
          [~,~,~,Stats]=RP.invkin2(zeros(6,1), q0, s);
          if Stats.version < fileversions.pkm_invkin % hier wird die aktuelle Version eingetragen
            error('Version der Datei ist zu alt (%d). Aktuell %d', ...
              Stats.version, fileversions.pkm_invkin);
          end
          % Debug: Prüfe ob Kinematik übereinstimmt
%           x0 = RP.fkineEE_traj(q0(:,1)')';
%           [q1,Phi1] = RP.invkin2(x0, q0(:,1), struct('retry_limit', 0, ...
%             'n_max', 20));
%           if any(abs(q0(1:RP.Leg(1).NJ,1)-q1(1:RP.Leg(1).NJ))>1e-6) || ...
%               any(abs(Phi1(1:sum(RP.I_EE_Task)))>1e-6)
%             error('InvKin und DirKin stimmen nicht überein');
%           end
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin_ser_mex') || ...
          contains(filelist(kk).name, 'invkin_ser.m')
        % Dummy-Eintrag in Dateiliste. Diese Datei gibt es nicht. Diese
        % Funktion ruft aber kompilierte SerRob-Funktionen auf.
        try
          s = struct('n_max', 1, 'retry_limit', -1); % keine Ausführung
          RP.invkin_ser(zeros(6,1), rand(RP.NJ,3), s);
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin3')
        if ~RP_mex_status % markiere die Datei als vorhanden (ohne Mex)
          files_found.pkm_invkin3 = true;
        end
        try
          s = struct('n_max', 1, 'retry_limit', -1); % Keine Durchführung des Algorithmus
          [~,~,~,Stats,~] = RP.invkin4(zeros(6,1), rand(RP.NJ,3), s);
          if Stats.version < fileversions.pkm_invkin3 % hier wird die aktuelle Version eingetragen
            error('Version der Datei ist zu alt (%d). Aktuell: %d', ...
              Stats.version, fileversions.pkm_invkin3);
          end
        catch err
          recompile = true;
        end
      end
      if contains(filelist(kk).name, 'invkin_traj')
        if ~RP_mex_status % markiere die Datei als vorhanden (ohne Mex)
          files_found.pkm_invkin_traj = true;
        end
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
      if contains(filelist(kk).name, 'fkine_coll')
        if ~RP_mex_status % markiere die Datei als vorhanden (ohne Mex)
          files_found.fkine_coll = true;
        end
        Tc_stack_cls = RP.fkine_coll(zeros(RP.NJ,1));
        try
          Tc_stack_tpl = RP.fkine_coll2(zeros(RP.NJ,1));
          assert(all(size(Tc_stack_cls)==size(Tc_stack_tpl)), ...
            'Dimension der Ausgabe von fkine_coll stimmt nicht');
        catch err
          recompile = true;
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
          serroblib_create_template_functions({LEG_Names{1}}, false, false); %#ok<CCAT1>
          parroblib_create_template_functions({PName_Kin},false,false);
          % Funktionen neu eintragen, falls vorher welche gefehlt haben
          RP.fill_fcn_handles(RP_mex_status, false);
        elseif retryiter == 3 && ...  % Kein Fehler für Dummy-Eintrag. 
            ~isempty(filelist(kk).bytes) % Dann fehlen einfach die Mex-Dateien der SerRob-Klasse. Ist nicht so schlimm.
          warning(orig_state);
          % Kompiliere nochmal und zeige auch den Bericht dazu zum Debuggen
          if contains(filelist(kk).name, '_mex') && ~isempty(filelist(kk).bytes)
            [~,mexbasename] = fileparts(filelist(kk).name);
            matlabfcn2mex({mexbasename(1:end-4)}, true, false, true);
          end
          error('Auch beim dritten Versuch kein Erfolg');
        end
        if contains(filelist(kk).name, '_mex')
          if ~isempty(filelist(kk).bytes) % nicht für Dummy-Eintrag machen
            % Mex-Dateien werden neu kompiliert. m-Dateien nicht.
            [~,mexbasename] = fileparts(filelist(kk).name);
            matlabfcn2mex({mexbasename(1:end-4)}); % ohne "_mex"
          else
            % Funktion für SerRob-Klasse neu generieren. Fehler kommt bei
            % Dummy-Eintrag (invkin_ser) daher.
            serroblib_update_template_functions(LEG_Names(1), verbosity, ignore_mex);
            % Funktions-Handles neu eintragen. Wenn SerRob-InvKin beim
            % ersten Versuch fehlte, sonst keine Aktualisierung.
            RP.fill_fcn_handles(RP_mex_status, false);
          end
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
  % Prüfe die Funktionen auf Vollständigkeit. Es gibt immer min. eine Datei
  for f = fields(files_found)'
    if ~all(files_found.(f{1}))
      % Es gibt aus Vorlagen generierte Dateien, aber manche Nicht- 
      % Mex-Dateien fehlen. Diese wurden manuell gelöscht. Neu erstellen
      if verbosity
        fprintf('Datei %s fehlt, aber %d Dateien vorhanden.\n', f{1}, length(filelist));
      end
      parroblib_create_template_functions({PName_Kin}, false, false);
      break;
    end
  end
end % ii (PKM)
warning(orig_state);
