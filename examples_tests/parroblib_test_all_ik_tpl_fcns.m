% Aufruf der IK-Funktionen aller existierende PKM in der Bibliothek
% Bestimme Kinematikparameter mit Maßsynthese,
% dann IK testen, Klassenmethode gegen Templatemethode.

% Junnan Li, WiHi bei moritz.schappler@imes.uni-hannover.de 2020-06
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Initialisierung
% Schalter zur Neugenerierung der Template-Funktionen
tpl_fcn_neu = true;
max_num_pkm = 20; % Reduziere die Anzahl der geprüften PKM pro FG
% Speicherort der Parameter
rob_path = fileparts(which('robotics_toolbox_path_init.m'));
tmpdir_params = fullfile(rob_path, 'examples_tests', 'tmp_ParRob', 'param_dimsynthres');
% Alle FG-Kombinationen durchgehen
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 1];
EE_FG_Mask = [1 1 1 1 1 1];

for i_FG = 1:size(EEFG_Ges,1)
  EE_FG = EEFG_Ges(i_FG,:);
  [PNames_Kin, ~] = parroblib_filter_robots(sum(EE_FG), EE_FG, EE_FG_Mask, 6);
  if isempty(PNames_Kin)
    continue % Es gibt keine PKM mit diesen FG.
  end
  III = 1:length(PNames_Kin);
  III = III(randperm(length(III))); % shuffle
  for ii = III(1:min(max_num_pkm, length(III)))
    PName = [PNames_Kin{ii},'A1']; % Nehme nur die erste Aktuierung (ist egal)
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    paramfile_robot = fullfile(tmpdir_params, sprintf('%s_params.mat', PName));
    fprintf('Untersuche PKM %s\n', PName);
    if tpl_fcn_neu 
      serroblib_create_template_functions({LEG_Names{1}},false,true)
      parroblib_create_template_functions({PNames_Kin{ii}},false,true);
    end
    RP = parroblib_create_robot_class(PName, 1, 0.3);
    % Initialisierung der Funktionen: Kompilierte Funktionen nehmen
    files_missing = RP.fill_fcn_handles(true, true);
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
    % Teste inverse Kinematik für Einzelpunkte
    for k = 1:size(Traj_0.X,1)
      % Teste inverse Kinematik für Einzelpunkte nach erster Methode
      [q_kls, Phi_kls, Tc_stack_kls]=RP.invkin_ser(Traj_0.X(k,:)', q0, s); % Klassen
      [q_tpl, Phi_tpl, Tc_stack_tpl]=RP.invkin2(Traj_0.X(k,:)', q0, s); % Template
      ik_res_ik2 = (all(abs(Phi_tpl(RP.I_constr_t_red))<s.Phit_tol) && ...
                    all(abs(Phi_tpl(RP.I_constr_r_red))<s.Phir_tol));% IK-Status Funktionsdatei
      ik_res_iks = (all(abs(Phi_kls(RP.I_constr_t_red))<s.Phit_tol) && ... 
                    all(abs(Phi_kls(RP.I_constr_r_red))<s.Phir_tol)); % IK-Status Klassenmethode
      if ik_res_ik2 ~= ik_res_iks % Vergleiche IK-Status (Erfolg / kein Erfolg)
        error('IK Status nicht gleich');
      elseif ~ik_res_iks  % beide IK-Status false
        error('IK Status beide nicht erfolgreich');
      end
      ik_test_Tc_stack = Tc_stack_kls - Tc_stack_tpl;
      if max(abs(ik_test_Tc_stack(:))) > 1e-3
        error('Ausgabe Tc_stack stimmt nicht überein');
      end
      test_q = q_kls-q_tpl;
      if any(abs(test_q)>1e-6)
        warning('Gelenkwinkel aus kls und tpl stimmen nicht überein');
      end
      % Teste inverse Kinematik für Einzelpunkte nach zweiter Methode
      % Geht eventuell nur für 3T3R

      if all(EE_FG == [1 1 1 1 1 1])
        [q_kls_3, Phi_kls_3, Tc_stack_kls_3]=RP.invkin3(Traj_0.X(k,:)', q0, s);
        [q_tpl_3, Phi_tpl_3, Tc_stack_tpl_3]=RP.invkin4(Traj_0.X(k,:)', q0, s);
        ik_res2_ik2 = (all(abs(Phi_tpl_3(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_tpl_3(RP.I_constr_r_red))<s.Phir_tol));% IK-Status Funktionsdatei
        ik_res2_iks = (all(abs(Phi_kls_3(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_kls_3(RP.I_constr_r_red))<s.Phir_tol)); % IK-Status Klassenmethode
        if ik_res2_ik2 ~= ik_res2_iks % Vergleiche IK-Status (Erfolg / kein Erfolg)
          error('IK Status nicht gleich');
        elseif ~ik_res2_iks  % beide IK-Status false
          error('IK Status beide nicht erfolgreich');
        end
        ik_test_Tc_stack_2 = Tc_stack_kls_3 - Tc_stack_tpl_3;
        if max(abs(ik_test_Tc_stack_2(:))) > 1e-3
          error('Ausgabe Tc_stack stimmt nicht überein');
        end
        test_q_2 = q_kls_3-q_tpl_3;
        if any(abs(test_q_2)>1e-6)
          warning('Gelenkwinkel aus kls und tpl stimmen nicht überein');
        end
      end
    end
    s_traj_ik = s;
    % Teste inverse Kinematik für Trajektorie nach erster Methode
    if all(EE_FG == [1 1 1 1 1 1])
      ik_modes = [1 2]; % hier geht die invkin3- und invkin_ser-Methode
    else
      ik_modes = 1; % hier nur invkin_ser-Methode
    end
    for ik = ik_modes
      s_traj_ik.mode_IK = ik;
      [Q1_t_kls, QD1_t_kls, ~, Phi1_t_kls,~,~, JointPos_all_kls] = RP.invkin_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0, s_traj_ik);% Klassen
      [Q1_t_tpl, QD1_t_tpl, ~, Phi1_t_tpl,~,~, JointPos_all_tpl] = RP.invkin2_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0, s_traj_ik);% Template
      test_joint = JointPos_all_kls - JointPos_all_tpl;
      test_Q1 = Q1_t_kls - Q1_t_tpl;
      if any(abs(test_joint)>1e-6)
        warning('JointPos aus kls und tpl stimmen nicht überein');
      end
      if any(abs(test_Q1)>1e-6)
        warning('Q1 aus kls und tpl stimmen nicht überein');
      end
    end
    fprintf('Klassen- und Template-Methode stimmen überein\n');
  end
end
