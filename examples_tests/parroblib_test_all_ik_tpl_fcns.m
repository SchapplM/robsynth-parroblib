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
%  % Debug:
%   if i_FG == 1
%     III = find(strcmp(PNames_Kin, 'P3RRR1G1P1'));
%   elseif i_FG == 2
%     III = find(strcmp(PNames_Kin, 'P3RRRRR10V1G2P2'));
%   elseif i_FG == 4
%     III = find(strcmp(PNames_Kin, 'P5RPRRR8V1G9P8'));
%   end
  for ii = III(1:min(max_num_pkm, length(III))) % Debug: find(strcmp(PNames_Kin, 'P6RRPRRR14V3G1P4'));
    PName = [PNames_Kin{ii},'A1']; % Nehme nur die erste Aktuierung (ist egal)
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    paramfile_robot = fullfile(tmpdir_params, sprintf('%s_params.mat', PName));
    fprintf('Untersuche PKM %s\n', PName);
    serroblib_create_template_functions({LEG_Names{1}},  ~tpl_fcn_neu,recompile_mex) %#ok<CCAT1>
    parroblib_create_template_functions({PNames_Kin{ii}},~tpl_fcn_neu,recompile_mex); %#ok<CCAT1>
%     matlabfcn2mex({[PName_Legs,'_invkin'], [PName_Legs,'_invkin3'], ...
%       [PName_Legs,'_invkin_traj']});
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
      Traj_0 = cds_transform_traj(RP, Traj_W);
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
      Set.optimization.objective = 'condition';
      Set.optimization.ee_rotation = false;
      Set.optimization.ee_translation = false;
      Set.optimization.movebase = false;
      Set.optimization.base_size = false;
      Set.optimization.platform_size = false;
      Set.optimization.obj_limit = 1e3; % Sofort abbrechen, falls Ergebnis gut
      Set.structures.use_parallel_rankdef = 6;
      Set.structures.whitelist = {PName}; % nur diese PKM untersuchen
      Set.structures.nopassiveprismatic = false; % Für Dynamik-Test egal 
      Set.structures.maxnumprismatic = 6; % Für Dynamik-Test egal wie viele Schubgelenke
      Set.general.noprogressfigure = true;
      Set.general.verbosity = 3;
      Set.general.nosummary = true;
      Traj = Traj_W;
      cds_start
      if isempty(Structures)
        error('PKM %s wurde erst aus Datenbank gewählt und danach nicht mehr gefunden. Fehler.', PName);
      end
      resmaindir = fullfile(Set.optimization.resdir, Set.optimization.optname);
      i_select = 0;
      for i = 1:length(Structures) % alle Ergebnisse durchgehen (falls mehrere theta-Varianten)
        resfile = fullfile(resmaindir, sprintf('Rob%d_%s_Endergebnis.mat', Structures{i}.Number, PName));
        tmp = load(resfile, 'RobotOptRes');
        if tmp.RobotOptRes.fval < 1000
          i_select = i;
          RobotOptRes = tmp.RobotOptRes;
          break;
        end
      end
      if isempty(Structures) || i_select == 0
        % Die Methode valid_act nimmt die erstbeste bestimmbare Kinematik.
        % Die Wahl der aktuierten Gelenke muss nicht zu vollem Rang führen.
        % Kriterium ist daher nur die Bestimmbarkeit des Rangs (fval <1000)
        warning('Etwas ist bei der Maßsynthese schiefgelaufen. Keine Lösung.');
        continue
      end
      RP = RobotOptRes.R;
      r_W_0 = RP.r_W_0;
      phi_W_0 = RP.phi_W_0;
      phi_P_E = RP.phi_P_E;
      r_P_E = RP.r_P_E;
      pkin = RP.Leg(1).pkin;
      DesPar_ParRob = RP.DesPar;
      q0 = RobotOptRes.q0;
      qlim = cat(1, RP.Leg.qlim); % Wichtig für Mehrfach-Versuche der IK
      save(paramfile_robot, 'pkin', 'DesPar_ParRob', 'q0', 'r_W_0', 'phi_W_0', 'qlim', 'r_P_E', 'phi_P_E');
      fprintf('Maßsynthese beendet\n');
      Traj_0 = cds_transform_traj(RP, Traj_W);
    end
    
    % Klassenmethode gegen Templatemethode
    s = struct('Phit_tol', 1e-9, 'Phir_tol', 1e-9, 'retry_limit', 0, ...
      'normalize', false, 'n_max', 5000);
    %% Teste inverse Kinematik für Einzelpunkte nach erster Methode
    fprintf('Vergleiche Klassenmethode (invkin_ser) und Template-Methode (invkin2) für Einzelpunkt-IK Variante 1\n');
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
      num_niO = 0;
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
          error('invkin_ser vs invkin2 (Punkt %d/%d): IK Status nicht gleich', i, max_single_points);
        elseif ~ik_res_iks  % beide IK-Status false
          num_niO = num_niO + 1;
          warning('invkin_ser vs invkin2 (Punkt %d/%d): IK Status beide nicht erfolgreich', i, max_single_points);
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
      if num_niO / max_single_points > 0.80 % bei ungünstigen Kinematikparametern evtl schlechte Konvergenz
        error('invkin_ser vs invkin2: Mehr als 80%% (%d/%d) der IK-Versuche für beide Implementierungen nicht erfolgreich', num_niO, max_single_points);
      end
      fprintf(['Klassen- und Template-Methode für Einzelpunkt-IK stimmen überein (Fall %d, %d Punkte).\n', ...
        'Zeiten: invkin_ser: mean %1.1fms (std %1.1fms), invkin2: mean %1.1fms (std %1.1fms)\n'], ...
        jj, max_single_points, 1e3*mean(calctimes(:,1)), ...
        1e3*std(calctimes(:,1)), 1e3*mean(calctimes(:,2)), 1e3*std(calctimes(:,2)));
    end
    fprintf('Vergleich invkin_ser vs invkin2 für 3 Fälle erfolgreich\n');
    %% Teste inverse Kinematik für Einzelpunkte nach zweiter Methode
    fprintf('Vergleiche Klassenmethode (invkin3) und Template-Methode (invkin4) für Einzelpunkt-IK Variante 2\n');
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
          s_jj.I_EE_Task = logical([1 1 1 1 1 0]);
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
    fprintf('Vergleich invkin3 vs invkin4 für 2 Fälle erfolgreich\n');
    %% Teste inverse Kinematik für Trajektorie
    fprintf('Vergleiche Klassenmethode (invkin_traj) und Template-Methode (invkin2_traj) für Trajektorien-IK\n');
    for jj = 1:5
      calctimes = NaN(1,2);
      % Zurücksetzen der Aufgaben-FG auf Standard-Wert
      RP.update_EE_FG(RP.I_EE, RP.I_EE);
      s_jj = s;
      s_jj.debug = true; % Abbruch, wenn interne Fehler in Traj.-IK-Fkt.
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
      % Berechne IK des ersten Traj.-Punktes, damit beide Verfahren
      % gleich anfangen (sonst mehr Möglichkeit für Abweichungen).
      [q0_traj, Phi_q0] = RP.invkin_ser(Traj_0.X(1,:)', q0);
      if any(abs(Phi_q0) > 1e-8) || any(isnan(Phi_q0))
        error('IK für Anfangspunkt nicht berechenbar. Muss funktionieren, da Parameter aus Maßsynthese');
      end
      t1=tic();
      [Q_t_kls, QD_t_kls, QDD_t_kls, Phi1_t_kls,~,~, JointPos_all_kls] = ...
        RP.invkin_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0_traj, s_jj);% Klassen
      calctimes(1)=toc(t1);
      t1=tic();
      [Q_t_tpl, QD_t_tpl, QDD_t_tpl, Phi1_t_tpl,~,~, JointPos_all_tpl] = ...
        RP.invkin2_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0_traj, s_jj);% Template
      calctimes(2)=toc(t1);
      trajik_iO = [any(abs(Phi1_t_kls(:)) > 1e-6), any(abs(Phi1_t_tpl(:)) > 1e-6)];
      if any(~trajik_iO)
        error('Eine Trajektorien-IK funktioniert nicht. Darf nicht sein, da Parameter aus Maßsynthese kommen');
      end
      % Prüfe das Ergebnis der Traj.-IK in sich. Wenn die IK erfolgreich
      % ist, muss die Plattform-Bewegung von allen Beinketten aus gezählt
      % übereinstimmen
      for iktyp = 1:2
        if ~trajik_iO(iktyp), continue; end % Nicht prüfbar, ob IK stimmt.
        if iktyp == 1
          Q = Q_t_kls; QD = QD_t_kls; QDD = QDD_t_kls;
        else
          Q = Q_t_tpl; QD = QD_t_tpl; QDD = QDD_t_tpl;
        end
        for j = 1:RP.NLEG
          [X3,XD3,XDD3] = RP.fkineEE_traj(Q, QD, QDD, j);
          if j == 1 % Speichere die erste Beinkette als Referenz
            X = X3; XD = XD3; XDD = XDD3;
          end
          test_X = X(:,1:6) - X3(:,1:6);
          test_X([false(size(test_X,1),3),abs(abs(test_X(:,4:6))-2*pi)<1e-3]) = 0; % 2pi-Fehler entfernen
          test_XD = XD(:,1:6) - XD3(:,1:6);
          test_XDD = XDD(:,1:6) - XDD3(:,1:6);
          if max(abs(test_X(:)))>1e-6
            Ifirst = find(any(abs(test_X)>1e-6,2), 1, 'first');
            error('Die Endeffektor-Trajektorie X aus Beinkette %d stimmt nicht gegen Beinkette 1. Erstes Vorkommnis: Zeitschritt %d', j, Ifirst);
          end
          if max(abs(test_XD(:)))>1e-6
            Ifirst = find(any(abs(test_XD)>1e-6,2), 1, 'first');
            error('Die Endeffektor-Trajektorie XD aus Beinkette %d stimmt nicht gegen Beinkette 1. Erstes Vorkommnis: Zeitschritt %d', j, Ifirst);
          end
          test_XDD(abs(test_XDD)<1e-3) = 0;
          if max(abs(test_XDD(:)))>1e-3
            Ifirst = find(any(abs(test_XDD)>1e-3,2), 1, 'first');
            error('Die Endeffektor-Trajektorie XDD aus Beinkette %d stimmt nicht gegen Beinkette 1. Erstes Vorkommnis: Zeitschritt %d', j, Ifirst);
          end
        end
      end
      test_JP = JointPos_all_kls - JointPos_all_tpl;
      test_Q = Q_t_kls - Q_t_tpl;
      if any(abs(test_Q(:))>1e-3)
        Ifirst = find(any(abs(test_Q)>1e-3,2), 1, 'first');
        warning('invkin_traj vs invkin2_traj: Q1 aus kls und tpl stimmen nicht überein. Max. Abweichung %1.1e. Erstes Vorkommnis: Zeitschritt %d', ...
          max(abs(test_Q(:))), Ifirst);
      end
      if any(abs(test_JP(:))>1e-6)
        warning('invkin_traj vs invkin2_traj: JointPos aus kls und tpl stimmen nicht überein. Max. Abweichung %1.1e', ...
          max(abs(test_JP(:))));
      end
      fprintf(['Klassen- und Template-Methode für Traj.-IK stimmen überein (Fall %d).\n', ...
        'Zeiten: invkin_traj: %1.2fs, invkin2_traj: %1.2fs\n'], ...
        jj, calctimes(1), calctimes(2));
    end
  end % for ii (PKM)
end % for i_FG
