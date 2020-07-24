% Aufruf der IK-Funktionen aller existierende PKM in der Bibliothek
% Bestimme Kinematikparameter mit Maßsynthese,
% dann IK testen, Klassenmethode gegen Templatemethode.
% (Jede Methode ist auf zwei Arten implementiert)

% Junnan Li (WiHi bei Moritz Schappler), 2020-06
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-07
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Initialisierung
% Schalter zur Neugenerierung der Template-Funktionen
tpl_fcn_neu = true;
test_mex = true; % Benutze mex-Funktionen beim Testen (schneller, aber kein Debuggen möglich)
recompile_mex = true; % nur notwendig, falls Template-Dateien geändert wurden.
max_num_pkm = 3; % Reduziere die Anzahl der geprüften PKM pro FG
shuffle_pkm_selection = true; % Zufällige Auswahl der PKM
max_single_points = 50; % Anzahl der zu prüfenden Einzelpunkte
% Speicherort der Parameter
rob_path = fileparts(which('robotics_toolbox_path_init.m'));
tmpdir_params = fullfile(rob_path, 'examples_tests', 'tmp_ParRob', 'param_dimsynthres');
mkdirs(tmpdir_params);
% Alle FG-Kombinationen durchgehen
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 0; ...
            1 1 1 1 1 1];
EE_FG_Mask = [1 1 1 1 1 1];

for i_FG = 1:size(EEFG_Ges,1)
  EE_FG = EEFG_Ges(i_FG,:);
  [PNames_Kin, ~] = parroblib_filter_robots(sum(EE_FG), EE_FG, EE_FG_Mask, 6);
  if isempty(PNames_Kin)
    continue % Es gibt keine PKM mit diesen FG.
  end
  III = 1:length(PNames_Kin);
  if shuffle_pkm_selection
    III = III(randperm(length(III)));
  end
%   II = find(strcmp(PNames_Kin, 'P5RPRRR5G1P8'));
  for ii = III(1:min(max_num_pkm, length(III))) % Debug: find(strcmp(PNames_Kin, 'P6RRPRRR14V3G1P4'));
    PName = [PNames_Kin{ii},'A1']; % Nehme nur die erste Aktuierung (ist egal)
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    paramfile_robot = fullfile(tmpdir_params, sprintf('%s_params.mat', PName));
    fprintf('Untersuche PKM %s\n', PName);
    serroblib_create_template_functions({LEG_Names{1}},  ~tpl_fcn_neu,recompile_mex) %#ok<CCAT1>
    parroblib_create_template_functions({PNames_Kin{ii}},~tpl_fcn_neu,recompile_mex); %#ok<CCAT1>
    RP = parroblib_create_robot_class(PName, 1, 0.3);
    % Initialisierung der Funktionen: Kompilierte Funktionen nehmen
    % (notwendig für Struktursynthese)
    RP.fill_fcn_handles(true, true); % zur Prüfung, ob kompilierte Funktionen vorhanden sind
    RP.fill_fcn_handles(test_mex, true); % zur Einstellung mex ja/nein
    
    %% Kinematikparameter durch Optimierung erzeugen (oder gespeichert laden)
    Set = cds_settings_defaults(struct('DoF', EE_FG));
    Set.task.Ts = 1e-2;
    Set.task.Tv = 1e-1;
    Set.task.profile = 1; % Zeitverlauf mit Geschwindigkeit
    Set.task.maxangle = 5*pi/180;
    Traj_W = cds_gen_traj(EE_FG, 1, Set.task);
    % Reduziere Punkte (geht dann schneller, aber auch schlechtere KinPar.)
    % Traj = timestruct_select(Traj, [1, 2]);
    params_success = false;
    % Prüfe, ob die Kinematikparameter schon optimiert wurden.
    if exist(paramfile_robot, 'file')
      params = load(paramfile_robot);
      q0 = params.q0;
      for il = 1:RP.NLEG
        RP.Leg(il).update_mdh(params.pkin); 
        RP.Leg(il).qlim = params.qlim(RP.I1J_LEG(il):RP.I2J_LEG(il),:);
      end
      RP.update_base(params.r_W_0, params.phi_W_0);
      RP.align_base_coupling(params.DesPar_ParRob.base_method, params.DesPar_ParRob.base_par);
      RP.align_platform_coupling(params.DesPar_ParRob.platform_method, params.DesPar_ParRob.platform_par(1:end-1));
      Traj_0 = cds_rotate_traj(Traj_W, RP.T_W_0);
      % Prüfe die Lösbarkeit der IK
      [q_test,Phi]=RP.invkin_ser(Traj_0.X(1,:)', q0);
      if all(abs(Phi)<1e-6) && ~any(isnan(Phi))
        fprintf('IK erfolgreich mit abgespeicherten Parametern gelöst\n');
        params_success = true; % Parameter für erfolgreiche IK geladen.
      else
        warning('IK mit abgespeicherten Parametern nicht lösbar.');
      end
    else
      fprintf('Es gibt keine abgespeicherten Parameter. Führe Maßsynthese durch\n');
    end
    if ~params_success
      % Führe Maßsynthese neu aus. Parameter nicht erfolgreich geladen
      Set.optimization.objective = 'valid_act';
      Set.optimization.ee_rotation = false;
      Set.optimization.ee_translation = false;
      Set.optimization.movebase = false;
      Set.optimization.base_size = false;
      Set.optimization.platform_size = false;
      Set.structures.use_parallel_rankdef = 6;
      Set.structures.whitelist = {PName}; % nur diese PKM untersuchen
      Set.structures.nopassiveprismatic = false; % Für Dynamik-Test egal 
      Set.structures.maxnumprismatic = 6; % Für Dynamik-Test egal wie viele Schubgelenke
      Set.general.noprogressfigure = true;
      Set.general.verbosity = 3;
      Set.general.nosummary = true;
      Traj = Traj_W;
      cds_start
      resmaindir = fullfile(Set.optimization.resdir, Set.optimization.optname);
      resfile = fullfile(resmaindir, sprintf('Rob%d_%s_Endergebnis.mat', 1, PName));
      load(resfile, 'RobotOptRes');
      if isempty(Structures) || RobotOptRes.fval > 1000
        % Die Methode valid_act nimmt die erstbeste bestimmbare Kinematik.
        % Die Wahl der aktuierten Gelenke muss nicht zu vollem Rang führen.
        % Kriterium ist daher nur die Bestimmbarkeit des Rangs (fval <1000)
        warning('Etwas ist bei der Maßsynthese schiefgelaufen');
        load('D:\HiWi_imes\Repos\parrob_mdlbib\examples_tests\failed_PKM.mat')
        failed_pkm = {failed_pkm{:},RobotOptRes.R.mdlname};
        save('failed_PKM.mat','failed_pkm');
        continue
      end
      RP = RobotOptRes.R;
      r_W_0 = RP.r_W_0;
      phi_W_0 = RP.phi_W_0;
      pkin = RP.Leg(1).pkin;
      DesPar_ParRob = RP.DesPar;
      q0 = RobotOptRes.q0;
      qlim = cat(1, RP.Leg.qlim); % Wichtig für Mehrfach-Versuche der IK
      save(paramfile_robot, 'pkin', 'DesPar_ParRob', 'q0', 'r_W_0', 'phi_W_0', 'qlim');
      fprintf('Maßsynthese beendet\n');
      Traj_0 = cds_rotate_traj(Traj_W, RP.T_W_0);
    end
    
    % Klassenmethode gegen Templatemethode
    s = struct('Phit_tol', 1e-9, 'Phir_tol', 1e-9, 'retry_limit', 0, ...
      'normalize', false, 'n_max', 5000);
    %% Teste inverse Kinematik für Einzelpunkte nach erster Methode
    % Wähle nur wenige Punkte der Trajektorie zum Testen der Einzelpunkt-IK
    II_traj = 1:size(Traj_0.X,1);
    II_traj = II_traj(randperm(length(II_traj)));
    max_single_points = min(max_single_points, length(II_traj));
    for jj = 1:3
      s_jj = s;
      q0_jj = q0;
      % Führe die IK mehrfach mit unterschiedlichen Einstelungen durch.
      % Dadurch werden möglichst viele unterschiedliche Testfälle in den
      % Funktionen abgedeckt
      switch jj
        case 1
          % Standard-Einstellungen
        case 2
          % Nehme IK-Ergebnis der ersten Beinkette als Anfangswert für folgende
          q0_jj(RP.I1J_LEG(2):end) = NaN;
        case 3
          if all(EE_FG == [1 1 1 1 1 1])
            s_jj.I_EE = logical([1 1 1 1 1 0]);
          else
            % 3T2R-Aufruf ergibt nur für 3T3R-PKM Sinn.
            continue
          end
      end
      calctimes = NaN(max_single_points,2); % zum Abspeichern der Rechenzeit
      for i = 1:max_single_points
        k = II_traj(i);
        t1=tic();
        [q_kls, Phi_kls, Tc_stack_kls]=RP.invkin_ser(Traj_0.X(k,:)', q0_jj, s); % Klassen
        calctimes(i,1)=toc(t1);
        t1=tic();
        [q_tpl, Phi_tpl, Tc_stack_tpl]=RP.invkin2(Traj_0.X(k,:)', q0_jj, s); % Template
        calctimes(i,2)=toc(t1);
        ik_res_ik2 = (all(abs(Phi_tpl(RP.I_constr_t_red))<s.Phit_tol) && ...
                      all(abs(Phi_tpl(RP.I_constr_r_red))<s.Phir_tol));% IK-Status Funktionsdatei
        ik_res_iks = (all(abs(Phi_kls(RP.I_constr_t_red))<s.Phit_tol) && ... 
                      all(abs(Phi_kls(RP.I_constr_r_red))<s.Phir_tol)); % IK-Status Klassenmethode
        if ik_res_ik2 ~= ik_res_iks % Vergleiche IK-Status (Erfolg / kein Erfolg)
          error('invkin_ser vs invkin2: IK Status nicht gleich');
        elseif ~ik_res_iks  % beide IK-Status false
          error('invkin_ser vs invkin2: IK Status beide nicht erfolgreich');
        end
        ik_test_Tc_stack = Tc_stack_kls - Tc_stack_tpl;
        if max(abs(ik_test_Tc_stack(:))) > 1e-3
          error('invkin_ser vs invkin2: Ausgabe Tc_stack stimmt nicht überein');
        end
        test_q = q_kls-q_tpl;
        if any(abs(test_q)>1e-6)
          warning('invkin_ser vs invkin2: Gelenkwinkel aus kls und tpl stimmen nicht überein');
        end
      end
      fprintf(['Klassen- und Template-Methode für Einzelpunkt-IK stimmen überein (Fall %d, %d Punkte).\n', ...
        'Zeiten: invkin_ser: mean %1.1fms (std %1.1fms), invkin2: mean %1.1fms (std %1.1fms)\n'], ...
        jj, max_single_points, 1e3*mean(calctimes(:,1)), ...
        1e3*std(calctimes(:,1)), 1e3*mean(calctimes(:,2)), 1e3*std(calctimes(:,2)));
    end
    %% Teste inverse Kinematik für Einzelpunkte nach zweiter Methode
    for jj = 1:2
      calctimes = NaN(max_single_points,2);
      % Geht vorerst nur für 3T3R
      if ~all(EE_FG == [1 1 1 1 1 1])
        break % Die folgenden Tests alle nur für 3T3R.
      end
      % Standard-FG-Auswahl
      RP.update_EE_FG(logical([1 1 1 1 1 1]), logical([1 1 1 1 1 1]));
      s_jj = s;
      switch jj
        case 1
          % Teste IK mit vollständigen Freiheitsgraden
          % Standard-Einstellungen: 3T3R-IK
        case 2
          % Teste IK mit Aufgabenredundanz
          % 3T2R-IK
          s_jj.I_EE = logical([1 1 1 1 1 0]);
          RP.update_EE_FG(logical([1 1 1 1 1 1]), logical([1 1 1 1 1 0]));
      end
      for i = 1:max_single_points
        k = II_traj(i);
        t1=tic();
        [q_kls_3, Phi_kls_3, Tc_stack_kls_3]=RP.invkin3(Traj_0.X(k,:)', q0, s_jj);
        calctimes(i,1)=toc(t1);
        t1=tic();
        [q_tpl_3, Phi_tpl_3, Tc_stack_tpl_3]=RP.invkin4(Traj_0.X(k,:)', q0, s_jj);
        calctimes(i,2)=toc(t1);
        ik_res2_ik2 = (all(abs(Phi_tpl_3(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_tpl_3(RP.I_constr_r_red))<s.Phir_tol));% IK-Status Funktionsdatei
        ik_res2_iks = (all(abs(Phi_kls_3(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_kls_3(RP.I_constr_r_red))<s.Phir_tol)); % IK-Status Klassenmethode
        if ik_res2_ik2 ~= ik_res2_iks % Vergleiche IK-Status (Erfolg / kein Erfolg)
          error('invkin3 vs invkin4: IK Status nicht gleich');
        elseif ~ik_res2_iks  % beide IK-Status false
          % Gebe keinen Fehler aus. In der Struktursynthese wird die
          % Einzelpunkt-IK nur für die Eckpunkte benutzt.
          warning(['invkin3 vs invkin4: Pose %d/%d: IK Status beide nicht ', ...
            'erfolgreich (Fall %d)'], k, size(Traj_0.X,1), jj);
        end
        ik_test_Tc_stack_2 = Tc_stack_kls_3 - Tc_stack_tpl_3;
        if max(abs(ik_test_Tc_stack_2(:))) > 5e-2
          warning(['invkin3 vs invkin4: Ausgabe Tc_stack stimmt nicht ', ...
            'überein. Max. Abweichung %1.1e'], max(abs(ik_test_Tc_stack_2(:))));
        end
        test_q_2 = q_kls_3-q_tpl_3;
        if any(abs(test_q_2)>5e-2)
          warning(['invkin3 vs invkin4: Gelenkwinkel aus kls und tpl stimmen ', ...
            'nicht überein. Max. Abweichung %1.1e'], max(abs(test_q_2)));
        end        
      end
      fprintf(['Klassen- und Template-Methode für Einzelpunkt-IK Variante 2 stimmen überein (Fall %d, %d Punkte).\n', ...
        'Zeiten: invkin3: mean %1.1fms (std %1.1fms), invkin4: mean %1.1fms (std %1.1fms)\n'], ...
        jj, max_single_points, 1e3*mean(calctimes(:,1)), 1e3*std(calctimes(:,1)), 1e3*mean(calctimes(:,2)), 1e3*std(calctimes(:,2)));
    end
    %% Teste inverse Kinematik für Trajektorie
    for jj = 1:5
      calctimes = NaN(1,2);
      % Zurücksetzen der Aufgaben-FG auf Standard-Wert
      if all(RP.I_EE_Task == logical([1 1 1 1 1 0])) && RP.NJ == 25
        RP.update_EE_FG(RP.I_EE, RP.I_EE, repmat(RP.I_EE, 5, 1));
      else
        RP.update_EE_FG(RP.I_EE, RP.I_EE);
      end
      s_jj = s;
      switch jj
        case 1
          % Teste IK mit vollständigen Freiheitsgraden
          % Standard-Einstellungen: 3T3R-IK
          s_jj.mode_IK = 1;
        case 2
          % Standard-Einstellungen: 3T3R-IK
          if all(EE_FG == [1 1 1 1 1 1])
            s_jj.mode_IK = 2;
          else
            continue
          end
        case 3
          % Vereinfachte Berechnung der Beschleunigung, vollständige FG
          s_jj.simplify_acc = true;
        case 4
          % Teste IK mit Aufgabenredundanz
          % 3T2R-IK
          if ~all(EE_FG == [1 1 1 1 1 1]), continue; end
          s_jj.mode_IK = 1;
          s_jj.I_EE_Task = logical([1 1 1 1 1 0]);
          RP.update_EE_FG(logical([1 1 1 1 1 1]), logical([1 1 1 1 1 0]));
        case 5
          % Vereinfachte Berechnung der Beschleunigung, Aufgabenredundanz
          if ~all(EE_FG == [1 1 1 1 1 1]), continue; end
          s_jj.simplify_acc = true;
          s_jj.I_EE_Task = logical([1 1 1 1 1 0]);
          RP.update_EE_FG(logical([1 1 1 1 1 1]), logical([1 1 1 1 1 0]));
      end
      t1=tic();
      [Q_t_kls, QD_t_kls, ~, Phi1_t_kls,~,~, JointPos_all_kls] = ...
        RP.invkin_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0, s_jj);% Klassen
      calctimes(1)=toc(t1);
      t1=tic();
      [Q_t_tpl, QD_t_tpl, ~, Phi1_t_tpl,~,~, JointPos_all_tpl] = ...
        RP.invkin2_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0, s_jj);% Template
      calctimes(2)=toc(t1);
      test_JP = JointPos_all_kls - JointPos_all_tpl;
      test_Q = Q_t_kls - Q_t_tpl;
      if any(abs(test_JP(:))>1e-6)
        warning('invkin_traj vs invkin2_traj: JointPos aus kls und tpl stimmen nicht überein. Max. Abweichung %1.1e', ...
          max(abs(test_JP(:))));
      end
      if any(abs(test_Q(:))>1e-6)
        warning('invkin_traj vs invkin2_traj: Q1 aus kls und tpl stimmen nicht überein. Max. Abweichung %1.1e', ...
          max(abs(test_Q(:))));
      end
      fprintf(['Klassen- und Template-Methode für Traj.-IK stimmen überein (Fall %d).\n', ...
        'Zeiten: invkin_traj: %1.2fs, invkin2_traj: %1.2fs\n'], ...
        jj, calctimes(1), calctimes(2));
    end
  end
end
