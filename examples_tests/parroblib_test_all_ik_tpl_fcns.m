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
resdir = fullfile(rob_path, 'examples_tests', 'results');
% Alle FG-Kombinationen durchgehen
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 0; ...
            1 1 1 1 1 1];
ResStat = table();
for i_FG = 1:size(EEFG_Ges,1)
  % Aufgabenredundanz ist bei 2T1R, 3T1R und 3T3R möglich.
  if any(i_FG == [1 3 5]), taskred_possible = true;
  else,                    taskred_possible = false; end
  EE_FG = logical(EEFG_Ges(i_FG,:));
  EE_FG_red = EE_FG; EE_FG_red(6) = 0;
  [PNames_Kin, ~] = parroblib_filter_robots(EE_FG, 6);
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
%   elseif i_FG == 3
%     III = find(strcmp(PNames_Kin, 'P4PRRRR3V1G1P2'));
%   elseif i_FG == 4
%     III = find(strcmp(PNames_Kin, 'P5RPRRR8V1G9P8'));
%   elseif i_FG == 5
%     III = find(strcmp(PNames_Kin, 'P6RRPRRR14V3G1P4'));
%   end
  for ii = III(1:min(max_num_pkm, length(III))) % Debug: find(strcmp(PNames_Kin, 'P6RRPRRR14V3G1P4'));
    PName = [PNames_Kin{ii},'A1']; % Nehme nur die erste Aktuierung (ist egal)
    if size(ResStat,1)>0 && any(strcmp(ResStat.Name, PName))
      warning('PKM %s steht schon in der Ergebnis-Tabelle. Überspringe', PName);
      continue
    end
    row_ii = {PName, 0, NaN(1,3), NaN(1,3), NaN(1,4), NaN(1,4), false, NaN(1,8)};
    ResStat = [ResStat; row_ii]; %#ok<AGROW>
    ResStat.Properties.VariableNames = {'Name', 'DimSynth_Erfolg', ...
      'IK_M1_Anteil_Identisch', 'IK_M1_Anteil_Erfolg', 'IK_M2_Anteil_Identisch', ...
      'IK_M2_Anteil_Erfolg', 'IK_Traj_Erfolg', 'IK_Traj_ErfolgDetail'};
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    paramfile_robot = fullfile(tmpdir_params, sprintf('%s_params.mat', PName));
    fprintf('Untersuche PKM %s\n', PName);
    RP = parroblib_create_robot_class(PName, 1, 0.3);
    % Erzeuge IK-Funktionen der Beinkette neu (automatische Auswahl des
    % Varianten-Modells, falls notwendig)
    serroblib_create_template_functions({LEG_Names{1}},  ~tpl_fcn_neu,false); %#ok<CCAT1>
    if recompile_mex
      matlabfcn2mex({[RP.Leg(1).mdlname,'_invkin_eulangresidual'], [RP.Leg(1).mdlname,'_invkin_traj']});
    end
    parroblib_create_template_functions({PNames_Kin{ii}},~tpl_fcn_neu,false); %#ok<CCAT1>
    % Nur die benötigten Funktionen neu kompilieren
    if recompile_mex
      matlabfcn2mex({[PName_Legs,'_invkin'], [PName_Legs,'_invkin3'], ...
        [PName_Legs,'_invkin_traj']});
    end
    % Initialisierung der Funktionen: Kompilierte Funktionen nehmen
    % (notwendig für Struktursynthese)
    RP.fill_fcn_handles(true, true); % zur Prüfung, ob kompilierte Funktionen vorhanden sind
    RP.fill_fcn_handles(test_mex, true); % zur Einstellung mex ja/nein
    
    %% Kinematikparameter durch Optimierung erzeugen (oder gespeichert laden)
    Set = cds_settings_defaults(struct('DoF', EE_FG));
    % Trajektorie definieren, die auch in der Struktursynthese eingesetzt
    % wird. Siehe parroblib_add_robots_symact.m
    % TODO: Trajektorien zentral ablegen. Spezialfälle zentral definieren.
    Set.task.Ts = 1e-2;
    Set.task.Tv = 1e-1;
    Set.task.profile = 1; % Zeitverlauf mit Geschwindigkeit
    Set.task.maxangle = 5*pi/180;
    if i_FG == 1,    trajno = 2; % 2T1R-Traj. muss auch Drehung haben. Sonst bei 3RPP Gelenkwinkel konstant.
    elseif i_FG < 4, trajno = 1;
    else,            trajno = 3; end
    if all(EE_FG==[1 1 1 1 1 0])
      Set.task.maxangle = 3*pi/180;
    end
    Traj_W = cds_gen_traj(EE_FG, trajno, Set.task);
    if all(EE_FG==[1 1 1 1 1 0])
      Traj_W.X(:,4:5) = Traj_W.X(:,4:5) + 4*pi/180;
      Traj_W.XE(:,4:5) = Traj_W.XE(:,4:5) + 4*pi/180;
    end
    % Reduziere Punkte (geht dann schneller, aber auch schlechtere KinPar.)
    % Traj = timestruct_select(Traj, [1, 2]);
    % Lade die bestehenden Parameter und prüfe, ob sie gültig sind oder mit
    % einer veralteten Version der Maßsynthese erstellt wurden.
    params_valid = false;
    if exist(paramfile_robot, 'file')
      params = load(paramfile_robot);
      if length(intersect(fieldnames(params), {'pkin','qlim','qDlim','qDDlim',...
          'r_W_0','phi_W_0','r_P_E','phi_P_E','DesPar_ParRob','q0'}))==10
        params_valid = true;
      end
    end
    params_success = false;
    % Prüfe, ob die Kinematikparameter schon optimiert wurden.
    if params_valid
      q0 = params.q0;
      for il = 1:RP.NLEG
        RP.Leg(il).update_mdh(params.pkin); 
        RP.Leg(il).qlim = params.qlim(RP.I1J_LEG(il):RP.I2J_LEG(il),:);
        RP.Leg(il).qDlim = params.qDlim(RP.I1J_LEG(il):RP.I2J_LEG(il),:);
        RP.Leg(il).qDDlim = params.qDDlim(RP.I1J_LEG(il):RP.I2J_LEG(il),:);
      end
      RP.update_base(params.r_W_0, params.phi_W_0);
      RP.update_EE(params.r_P_E, params.phi_P_E);
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
      Set.optimization.obj_limit = 1e3; % Sofort abbrechen, falls Ergebnis irgendwie funktionierend (hinsichtlich IK der Beingelenke)
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
        resfile1 = fullfile(resmaindir, sprintf('Rob%d_%s_Details.mat', Structures{i}.Number, PName));
        tmp1 = load(resfile1, 'RobotOptDetails');
        resfile1 = fullfile(resmaindir, sprintf('Rob%d_%s_Endergebnis.mat', Structures{i}.Number, PName));
        tmp2 = load(resfile1, 'RobotOptRes');
        if isfield(tmp1, 'RobotOptDetails') && tmp2.RobotOptRes.fval <= 1000 % Singuläre Konditionszahl der aktiven Gelenke heißt exakt 1e3 als Fitness-Wert. Reicht hier.
          i_select = i;
          RobotOptDetails = tmp1.RobotOptDetails;
          RobotOptRes = tmp2.RobotOptRes;
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
      RP = RobotOptDetails.R;
      r_W_0 = RP.r_W_0;
      phi_W_0 = RP.phi_W_0;
      phi_P_E = RP.phi_P_E;
      r_P_E = RP.r_P_E;
      pkin = RP.Leg(1).pkin;
      DesPar_ParRob = RP.DesPar;
      q0 = RobotOptRes.q0;
      qlim = cat(1, RP.Leg.qlim); % Wichtig für Mehrfach-Versuche der IK
      qDlim = cat(1, RP.Leg.qDlim); qDDlim = cat(1, RP.Leg.qDDlim);
      save(paramfile_robot, 'pkin', 'DesPar_ParRob', 'q0', 'r_W_0', ...
        'phi_W_0', 'qlim', 'qDlim', 'qDDlim', 'r_P_E', 'phi_P_E');
      fprintf('Maßsynthese beendet\n');
      Traj_0 = cds_transform_traj(RP, Traj_W);
    end
    RP.fill_fcn_handles(test_mex, true); % mex wird in Maßsynthese aktiviert
    % IK nochmal durchführen, falls eine Maßsynthese gemacht wurde und die
    % Parameter aktualisiert wurden.
    s_ik = struct('normalize', false);
    [q_test,Phi]=RP.invkin_ser(Traj_0.X(1,:)', q0, s_ik);
    assert(all(abs(Phi)<1e-6) && all(~isnan(Phi)), 'Phi stimmt nicht mit Startwert-IK');
    % Extrahiere Gelenk-Grenzen und passe
    qlim_pkm = cat(1,RP.Leg(:).qlim);
    I_inf = isinf(qlim_pkm);
    qlim_pkm(I_inf) = sign(qlim_pkm(I_inf))*pi;
    % Erzeuge 90%-Bereich innerhalb der Grenzen
    qlim_90 = repmat(mean(qlim_pkm,2),1,2) + repmat(diff(qlim_pkm')',1,2) .* ...
      repmat([-0.45, 0.45], size(qlim_pkm,1), 1);
    % Weite die Gelenkwinkelgrenzen auf. In Maßsynthese werden sie auf die
    % Ergebnis-Trajektorie eingestellt.
    if any(q0<qlim_90(:,1)) || any(q0 > qlim_90(:,2))
      q0_norm1 = (q0-qlim_pkm(:,1))./(qlim_pkm(:,2)-qlim_pkm(:,1));
      warning('IK-Anfangswerte q0 sind außerhalb der 90%-Grenzen. Erweitere die Grenzen.');
      % Generiere neue Grenzen, die mindestens 10% der Spannweite von dem
      % Startpunkt weg sind. Dadurch werden die ursprünglichen Grenzen
      % aufgeweitet, wenn die Anfangswinkel nahe daran sind.
      qlim_safe1 = repmat(q0,1,2) + repmat(diff(qlim_pkm')',1,2) .* ...
        repmat([-0.1, 0.1], size(qlim_pkm,1), 1);
      qlim_neu = minmax2([qlim_safe1, qlim_pkm]);
      for i = 1:RP.NLEG
        RP.Leg(i).qlim = qlim_neu(RP.I1J_LEG(i):RP.I2J_LEG(i),:);
      end
      q0_norm2 = (q0-qlim_neu(:,1))./(qlim_neu(:,2)-qlim_neu(:,1)); % muss zwischen 0 und 1 liegen
      fprintf(['Gelenkwinkel-Grenzen aktualisiert. Vorher q0 zwischen ', ...
        '%1.1f%% und %1.1f%% der alten Grenzen, jetzt zwischen %1.1f%% und %1.1f%% der neuen Grenzen.\n'], ...
        100*min(q0_norm1), 100*max(q0_norm1), 100*min(q0_norm2), 100*max(q0_norm2));
      assert(all(q0>=qlim_neu(:,1) & q0 <= qlim_neu(:,2)), ...
        'q0 außerhalb der Grenzen, obwohl gerade erst eingestellt');
    end
    qlim = cat(1,RP.Leg(:).qlim);
    
    ResStat.DimSynth_Erfolg(strcmp(ResStat.Name, PName)) = 1;
    % Allgemeine Einstellungen der IK:
    s = struct('Phit_tol', 1e-9, 'Phir_tol', 1e-9, 'retry_limit', 0, ...
      'normalize', false, 'n_max', 5000);
    % Prüfe die Konditionierung der PKM
    Phi_q = RP.constr3grad_q(q_test, Traj_0.X(1,:)');
    if cond(Phi_q) > 1e6
      warning('PKM ist bereits in Startpose stark singulär. cond(Phi_q)=%1.1e', cond(Phi_q));
      continue
    end
    %% Teste inverse Kinematik für Einzelpunkte nach erster Methode
    fprintf('Vergleiche Klassenmethode (invkin_ser) und Template-Methode (invkin2) für Einzelpunkt-IK Variante 1\n');
    % Wähle nur wenige Punkte der Trajektorie zum Testen der Einzelpunkt-IK
    II_traj = 1:size(Traj_0.X,1);
    II_traj = II_traj(randperm(length(II_traj)));
    max_single_points = min(max_single_points, length(II_traj));
    ikstat_neq_m1 = NaN(1,3);
    ikstat_niO_m1 = NaN(1,3);
    for jj = 1:3
      % Zurücksetzen der EE-FG (und Aufgaben-FG) des Roboters
      RP.update_EE_FG(EE_FG, EE_FG);
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
          if taskred_possible
            % s_jj.I_EE = logical([1 1 1 1 1 0]);
            % s_jj.scale_lim = 0;
            % TODO: Die IK bei 3T2R scheint schlechter lösbar zu sein.
            % Vermutlich führt das Ergebnis von Beinkette 1 zu einer
            % Lösung, die für eine der restlichen Beinketten besonders
            % ungünstig ist. Diese Methode ist für 3T2R aber sowieso nicht
            % wirklich zu empfehlen. Test wird nur gemacht, um die Gleich-
            % heit der Implementierungen zu prüfen.
            s_jj.retry_limit = 0; % mehr würde helfen, ist aber egal.
            RP.update_EE_FG(EE_FG, EE_FG_red);
          else
            % Aufruf mit Aufgabenredundanz nur für 2T1R/3T1R/3T3R sinnvoll
            continue
          end
      end
      calctimes = NaN(max_single_points,2); % zum Abspeichern der Rechenzeit
      num_niO = 0;
      num_neq = 0;
      for i = 1:max_single_points
        k = II_traj(i);
        t1=tic();
        [q_kls, Phi_kls, Tc_stack_kls, Stats_kls]=RP.invkin_ser(Traj_0.X(k,:)', q0_jj, s_jj); % Klassen
        calctimes(i,1)=toc(t1);
        t1=tic();
        [q_tpl, Phi_tpl, Tc_stack_tpl, Stats_tpl]=RP.invkin2(Traj_0.X(k,:)', q0_jj, s_jj); % Template
        calctimes(i,2)=toc(t1);
        ik_res_iks = (all(abs(Phi_kls(RP.I_constr_t_red))<s.Phit_tol) && ... 
                      all(abs(Phi_kls(RP.I_constr_r_red))<s.Phir_tol)); % IK-Status Klassenmethode
        ik_res_ik2 = (all(abs(Phi_tpl(RP.I_constr_t_red))<s.Phit_tol) && ...
                      all(abs(Phi_tpl(RP.I_constr_r_red))<s.Phir_tol));% IK-Status Funktionsdatei
        % Ergebnis der IK (Gültigkeit der Lösung) neu nachrechnen
        Phi_kls_v2 = RP.constr3(q_kls, Traj_0.X(k,:)');
        Phi_tpl_v2 = RP.constr3(q_tpl, Traj_0.X(k,:)');
        % Die Toleranz wird gröber gewählt, da in der invkin_ser die
        % Zwangsbedingungen nur nach der Beinkette (und nicht nach PKM)
        % berechnet werden. Die genaue Höhe der möglichen Abweichung ist
        % noch unklar. Daher kein Fehler.
        ik_res2_ik2_v2 = (all(abs(Phi_tpl_v2(RP.I_constr_t_red))<10*s.Phit_tol) && ...
                          all(abs(Phi_tpl_v2(RP.I_constr_r_red))<10*s.Phir_tol));
        ik_res2_iks_v2 = (all(abs(Phi_kls_v2(RP.I_constr_t_red))<10*s.Phit_tol) && ...
                          all(abs(Phi_kls_v2(RP.I_constr_r_red))<10*s.Phir_tol));
        if ik_res_ik2 ~= ik_res2_ik2_v2
          warning(['%d/%d: Ausgabe Phi von invkin2 stimmt nicht gegen ', ...
            'Prüfung mit constr3'], i, max_single_points);
        end
        if ik_res_iks ~= ik_res2_iks_v2
          warning(['%d/%d: Ausgabe Phi von invkin_ser stimmt nicht gegen ', ...
            'Prüfung mit constr3'], i, max_single_points);
        end
        if ik_res_ik2 ~= ik_res_iks % Vergleiche IK-Status (Erfolg / kein Erfolg)
          % Debuggen:
          change_current_figure(7);clf;
          RPstr = ['R', 'P'];
          for kkk = 1:RP.NJ
            legnum = find(kkk>=RP.I1J_LEG, 1, 'last');
            legjointnum = kkk-(RP.I1J_LEG(legnum)-1);
            subplot(ceil(sqrt(RP.NJ)), ceil(RP.NJ/ceil(sqrt(RP.NJ))), kkk);
            hold on; grid on;
            plot(Stats_kls.Q(1,kkk), 'ko');
            hdl1=plot(Stats_kls.Q(:,kkk), '-');
            hdl2=plot(Stats_tpl.Q(:,kkk), '--');
            plot([1;max(Stats_tpl.iter)], repmat(RP.Leg(legnum).qlim(legjointnum,:),2,1), 'r-');
            title(sprintf('q %d (%s), L%d,J%d', kkk, RPstr(RP.MDH.sigma(kkk)+1), legnum, legjointnum));
            grid on;
          end
          legend([hdl1;hdl2], {'invkin_ser', 'invkin2'}, 'interpreter', 'none');
          linkxaxes
          change_current_figure(8);clf;
          for kkk = 1:RP.NLEG
            subplot(2,3,kkk); hold on; grid on;
            hdl1=plot(Stats_kls.PHI(:,(kkk-1)*6+1:kkk*6), '-');  set(gca, 'ColorOrderIndex', 1);
            hdl2=plot(Stats_tpl.PHI(:,(kkk-1)*6+1:kkk*6), '--');
            title(sprintf('Phi Leg %d', kkk));
          end
          legend([hdl1(1);hdl2(1)], {'invkin_ser', 'invkin2'}, 'interpreter', 'none');
          linkxaxes
          change_current_figure(9);clf;
          for kkk = 1:RP.NLEG
            subplot(2,RP.NLEG,sprc2no(2,RP.NLEG,1,kkk)); hold on; grid on;
            plot(Stats_kls.condJ(:,kkk), '-');
            plot(Stats_tpl.condJ(:,kkk), '-');
            ylabel('Cond. Legs');
            title(sprintf('Leg %d', kkk));
            subplot(2,RP.NLEG,sprc2no(2,RP.NLEG,2,kkk)); hold on; grid on;
            plot(Stats_kls.lambda(:,kkk), '-');
            plot(Stats_tpl.lambda(:,kkk), '-');
            ylabel('Lambda. Legs');
          end
          linkxaxes
          error(['invkin_ser vs invkin2 (Punkt %d/%d): IK Status nicht ', ...
            'gleich. Klassenmethode: %d, Vorlagen-Funktion: %d'], i, ...
            max_single_points, ik_res_iks, ik_res_ik2);
        elseif ~ik_res_iks  % beide IK-Status false
          num_niO = num_niO + 1; % Darf passieren, aber nicht zu oft (s.u.)
          warning('invkin_ser vs invkin2 (Punkt %d/%d): IK Status beide nicht erfolgreich', i, max_single_points);
        end
        test_q = q_kls-q_tpl;
        if any(abs(test_q)>1e-6)
          num_neq = num_neq + 1;
          % Ein Umklappen der IK-Konfiguration ist möglich, wenn z.B.
          % aufgrund von Euler-Winkel-Singularitäten sich die Winkel ändern
%           figure(1);clf;
%           subplot(1,2,1);hold on;grid on;view(3); title('Klasse');
%           RP.plot( q_kls, Traj_0.X(k,:)');
%           subplot(1,2,2);hold on;grid on;view(3); title('Vorlage');
%           RP.plot( q_tpl, Traj_0.X(k,:)');
          warning('invkin_ser vs invkin2: Gelenkwinkel aus kls und tpl stimmen nicht überein');
        end
        ik_test_Tc_stack = Tc_stack_kls - Tc_stack_tpl;
        if all(abs(test_q)<1e-6) && max(abs(ik_test_Tc_stack(:))) > 1e-3
          % Die Ausgabe muss übereinstimmen, wenn die Winkel identisch sind
          error('invkin_ser vs invkin2: Ausgabe Tc_stack stimmt nicht überein');
        end
        
        % Teste IK bezogen auf Plattform-KS
        if jj == 3
          continue
        end
        % Nehme Ergebnis der EE-IK als Startwert. Die IK sollte sofort
        % konvergieren, ohne dass sich q ändert.
        [q_plf, Phi_plf, Tc_stack_plf, Stats_plf]=RP.invkin_ser( ...
          RP.xE2xP(Traj_0.X(k,:)'), q_kls, s_jj, struct('platform_frame', true));
        if any(abs(q_plf-q_kls)>1e-3)
          error('invkin_ser mit Plattform-KS stimmt nicht mit IK bezogen auf EE-KS überein');
        end
      end
      fprintf(['Klassen- und Template-Methode für Einzelpunkt-IK stimmen überein (Fall %d, %d Punkte).\n', ...
        'Zeiten: invkin_ser: mean %1.1fms (std %1.1fms), invkin2: mean %1.1fms (std %1.1fms)\n'], ...
        jj, max_single_points, 1e3*mean(calctimes(:,1)), ...
        1e3*std(calctimes(:,1)), 1e3*mean(calctimes(:,2)), 1e3*std(calctimes(:,2)));
      ikstat_neq_m1(jj) = num_niO / max_single_points;
      ikstat_niO_m1(jj) = num_neq / max_single_points;
      if jj == 1 % Nur bei Nr. 1 bei Fehler abbrechen. Nr. 2 hat evtl durch die anderen Startwerte eine schlechte Konvergenz. Bei Nr. 3 (3T2R) noch nicht sehr robust
        if ikstat_niO_m1(jj) > 0.80 % bei ungünstigen Kinematikparametern evtl schlechte Konvergenz
          error(['invkin_ser vs invkin2: Mehr als 80%% (%d/%d) der IK-Versuche ', ...
            'für beide Implementierungen nicht erfolgreich'], num_niO, max_single_points);
        end
        if ikstat_neq_m1(jj) > 0.50 % bei ungünstigen Kinematikparametern evtl schlechte Konvergenz
          error(['invkin_ser vs invkin2: Mehr als 50%% (%d/%d) der IK-Versuche ', ...
            'für beide Implementierungen nicht gleich'], num_niO, max_single_points);
        end
      end
    end
    ResStat.IK_M1_Anteil_Identisch(strcmp(ResStat.Name, PName),:) = 1-ikstat_neq_m1;
    ResStat.IK_M1_Anteil_Erfolg(strcmp(ResStat.Name, PName),:) = 1-ikstat_niO_m1;
    fprintf('Vergleich invkin_ser vs invkin2 für 3 Fälle erfolgreich\n');
    %% Teste inverse Kinematik für Einzelpunkte nach zweiter Methode
    fprintf('Vergleiche Klassenmethode (invkin3) und Template-Methode (invkin4) für Einzelpunkt-IK Variante 2\n');
    ikstat_neq_m2 = NaN(1,4);
    ikstat_niO_m2 = NaN(1,4);
    for jj = 1:4
      calctimes = NaN(max_single_points,2);
      % Geht vorerst nur für 3T3R
      if ~taskred_possible
        break % Die folgenden Tests nur für 2T1R/3T1R/3T3R. TODO: Wieso eigentlich?
      end
      % Standard-FG-Auswahl
      RP.update_EE_FG(EE_FG, EE_FG);
      s_jj = s;
      s_jj.n_max = 800; % damit schneller fertig. Muss reichen.
      s_jj.Kn = ones(RP.NJ,1);
      s_jj.maxstep_ns = 1e-6; % damit schneller aufgehört wird (Nullraumbewegung erst fast abgeschlossen)
      switch jj
        case 1
          % Teste IK mit vollständigen Freiheitsgraden
          % Standard-Einstellungen: 3T3R-IK
          s_jj.scale_lim = 0; % keine Skalierung an Grenze sinnvoll
        case 2
          % Teste IK mit Aufgabenredundanz
          % 3T2R-IK
          s_jj.scale_lim = 0.0; % Keine Skalierung an Grenze
          s_jj.wn = [0;0]; % Ohne Optimierung
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 3
          % Teste IK mit Aufgabenredundanz
          % 3T2R-IK
          s_jj.scale_lim = 0.0;
          s_jj.wn = [1;0;0]; % Mit Optimierung der Gelenkgrenzen
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 4
          % Teste 3T2R-IK mit Aufgabenredundanz (und Optimierung)
          s_jj.scale_lim = 0.0;
          s_jj.wn = [0;0;1]; % Mit Optimierung der Konditionszahl (ohne Betrachtung der Grenzen)
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 5
          % Teste 3T2R-IK mit Aufgabenredundanz (und Optimierung)
          % TODO: Das konvergiert im allgemeinen nicht wirklich gut.
          % Ohne scale_lim werden die Grenzen verletzt, mit gibt es keine
          % Lösung.
          s_jj.scale_lim = 0.8; % Keine Überschreitung der Grenzen erlauben (wegen hyperb. Bestrafung)
          s_jj.wn = [1;1;1]; % Mit Optimierung der Konditionszahl (mit Grenzen)
          RP.update_EE_FG(EE_FG, EE_FG_red);
      end
      num_niO = 0;
      num_neq = 0;
      for i = 1:max_single_points
        k = II_traj(i);
        t1=tic();
        [q_kls_3, Phi_kls_3, Tc_stack_kls_3, Stats_kls_3]=RP.invkin3(Traj_0.X(k,:)', q0, s_jj);
        calctimes(i,1)=toc(t1);
%         if any(q_kls_3<qlim(:,1)) || any(q_kls_3>qlim(:,2))
%           warning('Fall %d: %d/%d. Ergebnis überschreitet Gelenkwinkel-Grenzen', ...
%             jj, i, max_single_points);
%         end
        t1=tic();
        [q_tpl_3, Phi_tpl_3, Tc_stack_tpl_3, Stats_tpl_3]=RP.invkin4(Traj_0.X(k,:)', q0, s_jj);
        calctimes(i,2)=toc(t1);
        if i == 1
          fprintf(['Fall %d: Zeiten für eine Auswertung: %1.1fms für ', ...
            'invkin3 und %1.1fms für invkin4\n'], jj, 1e3*calctimes(i,1), 1e3*calctimes(i,2));
        end
        ik_res2_ik2 = (all(abs(Phi_tpl_3(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_tpl_3(RP.I_constr_r_red))<s.Phir_tol));% IK-Status Funktionsdatei
        ik_res2_iks = (all(abs(Phi_kls_3(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_kls_3(RP.I_constr_r_red))<s.Phir_tol)); % IK-Status Klassenmethode
        % Ergebnis der IK (Gültigkeit der Lösung) neu nachrechnen
        Phi_kls_3_v2 = RP.constr3(q_kls_3, Traj_0.X(k,:)');
        Phi_tpl_3_v2 = RP.constr3(q_tpl_3, Traj_0.X(k,:)');
        ik_res2_ik2_v2 = (all(abs(Phi_tpl_3_v2(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_tpl_3_v2(RP.I_constr_r_red))<s.Phir_tol));
        ik_res2_iks_v2 = (all(abs(Phi_kls_3_v2(RP.I_constr_t_red))<s.Phit_tol) && ...
          all(abs(Phi_kls_3_v2(RP.I_constr_r_red))<s.Phir_tol));
        if ik_res2_ik2 ~= ik_res2_ik2_v2
          error('Ausgabe Phi von invkin4 stimmt nicht bei erneuter Berechnung');
        end
        if ik_res2_iks ~= ik_res2_iks_v2
          error('Ausgabe Phi von invkin3 stimmt nicht bei erneuter Berechnung');
        end
        if ik_res2_ik2 ~= ik_res2_iks % Vergleiche IK-Status (Erfolg / kein Erfolg)
          error('invkin3 vs invkin4: IK Status nicht gleich');
        elseif ~ik_res2_iks  % beide IK-Status false
          % Gebe keinen Fehler aus. In der Struktursynthese wird die
          % Einzelpunkt-IK nur für die Eckpunkte benutzt.
          num_niO = num_niO + 1;
          warning(['invkin3 vs invkin4 %d/%d: Pose %d/%d: IK Status beide nicht ', ...
            'erfolgreich (Fall %d)'], i, max_single_points, k, size(Traj_0.X,1), jj);
          continue % TODO: Entscheiden, was passieren soll.
          % Debuggen:
          change_current_figure(10);clf;
          RPstr = ['R', 'P'];
          for kkk = 1:RP.NJ
            legnum = find(kkk>=RP.I1J_LEG, 1, 'last');
            legjointnum = kkk-(RP.I1J_LEG(legnum)-1);
            subplot(ceil(sqrt(RP.NJ)), ceil(RP.NJ/ceil(sqrt(RP.NJ))), kkk);
            hold on; grid on;
            plot(Stats_kls_3.Q(1,kkk), 'ko');
            plot(Stats_kls_3.Q(:,kkk), '-');
            hdl1=plot(Stats_tpl_3.Q(:,kkk), '-');
            hdl2=plot([1,Stats_kls_3.iter], repmat(RP.Leg(legnum).qlim(legjointnum,:),2,1), 'r:');
            title(sprintf('q %d (%s), L%d,J%d', kkk, RPstr(RP.MDH.sigma(kkk)+1), legnum, legjointnum));
            grid on;
          end
          linkxaxes
          legend([hdl1(1);hdl2(1)], {'invkin3', 'invkin4'}, 'interpreter', 'none');
          change_current_figure(11);clf;
          for kkk = 1:RP.NLEG
            Ic_Leg = RP.I_constr_red(RP.I_constr_red>(kkk-1)*6 & RP.I_constr_red<(kkk*6+1));
            subplot(3,3,kkk); hold on; grid on;
            plot(Stats_kls_3.PHI(:,Ic_Leg), '-');
            plot(Stats_tpl_3.PHI(:,Ic_Leg), '--');
            title(sprintf('Phi Leg %d', kkk));
            subplot(3,3,RP.NLEG+1); hold on; grid on;
            plot(sum(Stats_kls_3.PHI(:,Ic_Leg).^2,2), '-');
            plot(sum(Stats_tpl_3.PHI(:,Ic_Leg).^2,2), '--');
            title('Norm Phi Leg');
          end
          subplot(3,3,RP.NLEG+2); hold on; grid on;
          hdl1=plot(sum(Stats_kls_3.PHI(:,RP.I_constr_red).^2,2), '-');
          hdl2=plot(sum(Stats_tpl_3.PHI(:,RP.I_constr_red).^2,2), '--');
          legend([hdl1(1);hdl2(1)], {'invkin3', 'invkin4'}, 'interpreter', 'none');
          title('Norm Phi');
          linkxaxes
          change_current_figure(12);clf;
          for kkk = 1:4
            if kkk == 1, hlabel = 'hges';
            else, hlabel = sprintf('h %d', kkk-1); end
            subplot(2,4,sprc2no(2,4,1,kkk)); hold on; grid on;
            plot(Stats_kls_3.h(:,kkk), '-');
            plot(Stats_tpl_3.h(:,kkk), '--');
            ylabel(hlabel);
            subplot(2,4,sprc2no(2,4,2,kkk)); hold on; grid on;
            plot(diff(Stats_kls_3.h(:,kkk)), '-');
            plot(diff(Stats_tpl_3.h(:,kkk)), '--');
            ylabel(sprintf('diff %s', hlabel));
          end
          sgtitle(sprintf('Zielfunktion Nebenoptimierung'));
          linkxaxes
          error('Halt');
        end
        test_q_2 = q_kls_3-q_tpl_3;
        if any(abs(test_q_2)>5e-2)
          num_neq = num_neq + 1;
          warning(['invkin3 vs invkin4 %d/%d: Gelenkwinkel aus kls und tpl stimmen ', ...
            'nicht überein. Max. Abweichung %1.1e'], i, max_single_points, max(abs(test_q_2)));
        end 
        ik_test_Tc_stack_2 = Tc_stack_kls_3 - Tc_stack_tpl_3;
        if all(abs(test_q_2)<5e-2) && max(abs(ik_test_Tc_stack_2(:))) > 5e-2
          % Muss identisch sein, wenn die Winkel gleich sind.
          error(['invkin3 vs invkin4: Ausgabe Tc_stack stimmt nicht ', ...
            'überein. Max. Abweichung %1.1e'], max(abs(ik_test_Tc_stack_2(:))));
        end
      end
      fprintf(['Klassen- und Template-Methode für Einzelpunkt-IK Variante 2 stimmen überein (Fall %d, %d Punkte).\n', ...
        'Zeiten: invkin3: mean %1.1fms (std %1.1fms), invkin4: mean %1.1fms (std %1.1fms)\n'], ...
        jj, max_single_points, 1e3*mean(calctimes(:,1)), 1e3*std(calctimes(:,1)), 1e3*mean(calctimes(:,2)), 1e3*std(calctimes(:,2)));
      ikstat_neq_m2(jj) = num_niO / max_single_points;
      ikstat_niO_m2(jj) = num_neq / max_single_points;
      % Prüfe, ob IK-Ergebnis zufriedenstellend.
      if jj ~= 4 && ikstat_niO_m2(jj) > 0.80 % bei ungünstigen Kinematikparametern evtl schlechte Konvergenz
        error(['invkin3 vs invkin4: Mehr als 80%% (%d/%d) der IK-Versuche ', ...
          'für beide Implementierungen nicht erfolgreich'], num_niO, max_single_points);
      end
      if jj ~= 4 && ikstat_neq_m2(jj) > 0.50 % bei ungünstigen Kinematikparametern evtl schlechte Konvergenz
        error(['invkin3 vs invkin4: Mehr als 50%% (%d/%d) der IK-Versuche ', ...
          'für beide Implementierungen nicht gleich'], num_niO, max_single_points);
      end
    end
    ResStat.IK_M2_Anteil_Identisch(strcmp(ResStat.Name, PName),:) = 1-ikstat_neq_m2;
    ResStat.IK_M2_Anteil_Erfolg(strcmp(ResStat.Name, PName),:) = 1-ikstat_niO_m2;
    fprintf('Vergleich invkin3 vs invkin4 für 2 Fälle erfolgreich\n');
    %% Teste inverse Kinematik für Trajektorie
    % Debug:
%     RP.fill_fcn_handles(false);
%     parroblib_create_template_functions({PName}, false);
%     matlabfcn2mex({[PName(1:end-6), '_invkin_traj']});
    fprintf('Vergleiche Klassenmethode (invkin_traj) und Template-Methode (invkin2_traj) für Trajektorien-IK\n');
    trajikstat = NaN(1,8);
    for jj = 1:8
      calctimes = NaN(1,2);
      % Zurücksetzen der Aufgaben-FG auf Standard-Wert
      RP.update_EE_FG(RP.I_EE, RP.I_EE);
      s_jj = rmfield(s, 'retry_limit'); % für Traj.-IK Optionen anpassen.
      s_jj.debug = true; % Abbruch, wenn interne Fehler in Traj.-IK-Fkt.
      s_jj.Phir_tol = 1e-12; % Sehr feine Toleranz, damit Fehler-Schwellwerte ...
      s_jj.Phit_tol = 1e-12; % ... nicht dadurch überschritten werden
      s_jj.simplify_acc = false;
      switch jj
        case 1
          % Teste IK mit vollständigen Freiheitsgraden
          % Standard-Einstellungen: 3T3R-IK
          s_jj.mode_IK = 1; % benutze invkin_ser für Positionskorrektur
        case 2
          % Standard-Einstellungen: 3T3R-IK
          s_jj.mode_IK = 2; % benutze invkin3 für Positionskorrektur
        case 3
          % Vereinfachte Berechnung der Beschleunigung, vollständige FG
          s_jj.simplify_acc = true;
        case 4
          % Teste IK mit Aufgabenredundanz (ohne Nebenbedingung)
          if ~taskred_possible, continue; end
          s_jj.mode_IK = 3;
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 5
          % Teste IK mit Aufgabenredundanz ohne Nebenbedingungen
          if ~taskred_possible, continue; end
          s_jj.mode_IK = 3;
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 6
          % Vereinfachte Berechnung der Beschleunigung, Aufgabenredundanz
          % ohne Nebenbedingungen
          if ~taskred_possible, continue; end
          s_jj.simplify_acc = true;
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 7
          % Teste IK mit Aufgabenredundanz (mit Nebenbedingung)
          if ~taskred_possible, continue; end
          s_jj.mode_IK = 3;
          s_jj.wn = [1;0;1;0;0]; % quadratische Abstandsfunktion von Grenzen (Pos./Geschw.)
          RP.update_EE_FG(EE_FG, EE_FG_red);
        case 8
          % Teste IK mit Aufgabenredundanz (mit Nebenbedingung)
          if ~taskred_possible, continue; end
          s_jj.mode_IK = 3; % Benutze beide IK-Verfahren für Positions-Korrektur
          s_jj.wn = [0;0;0;0;1]; % Konditionszahl
          RP.update_EE_FG(EE_FG, EE_FG_red);
      end
      fprintf('Fall %d: Starte Trajektorien-IK (%d Bahnpunkte)\n', jj, size(Traj_0.X,1));
      % Berechne IK des ersten Traj.-Punktes, damit beide Verfahren
      % gleich anfangen (sonst mehr Möglichkeit für Abweichungen).
      [q0_traj, Phi_q0, ~, Stats] = RP.invkin_ser(Traj_0.X(1,:)', q0, s);
      Stats.PHI(isnan(Stats.PHI))=0; % Für Summenbildung unten
      [q0_traj2, Phi_q02, ~, Stats2] = RP.invkin3(Traj_0.X(1,:)', q0, s);
      [~,Phivoll_q0_test] = RP.constr3(q0_traj2, Traj_0.X(1,:)');
      Phi_q0_test = Phivoll_q0_test(RP.I_constr_red);
      if any(abs(Phi_q0) > 1e-6) || any(isnan(Phi_q0)) % invkin_ser geht nicht
        error('IK für Anfangspunkt nicht berechenbar. Muss funktionieren, da Parameter aus Maßsynthese');
      end
      if any(abs(Phi_q02) > 1e-6) ~= any(abs(Phi_q0_test) > 1e-6)
        error('Berechnung der kinematischen ZB am Anfang mit constr3 ist anders als aus invkin3.');
      end
      % invkin3 geht nicht, aber invkin_ser
      if all(abs(Phi_q0) < 1e-6) && any(abs(Phi_q02) > 1e-6)
        warning('IK zum Startpunkt der Trajektorie stimmt nicht mit invkin3, aber mit invkin_ser');
        if false % Debuggen
          figure(1); clf; hold on; grid on;
          xlabel('x in m');ylabel('y in m');zlabel('z in m'); view(3);
          s_plot = struct( 'ks_legs', [RP.I1L_LEG; RP.I1L_LEG+1; RP.I2L_LEG], 'straight', 0);
          RP.plot( q0_traj, Traj_0.X(1,:)', s_plot );
          change_current_figure(13);clf;
          RPstr = ['R', 'P'];
          for kkk = 1:RP.NJ
            legnum = find(kkk>=RP.I1J_LEG, 1, 'last');
            legjointnum = kkk-(RP.I1J_LEG(legnum)-1);
            subplot(ceil(sqrt(RP.NJ)), ceil(RP.NJ/ceil(sqrt(RP.NJ))), kkk);
            hold on; grid on;
            plot(Stats.Q(1,kkk), 'ko');
            hdl1=plot(Stats.Q(:,kkk), '-');
            hdl2=plot(Stats2.Q(:,kkk), '--');
            plot([1,max([Stats.iter,Stats2.iter])], repmat(RP.Leg(legnum).qlim(legjointnum,:),2,1), 'r-');
            title(sprintf('q %d (%s), L%d,J%d', kkk, RPstr(RP.MDH.sigma(kkk)+1), legnum, legjointnum));
            grid on;
          end
          legend([hdl1;hdl2], {'invkin_ser', 'invkin3'}, 'interpreter', 'none');
          linkxaxes
          change_current_figure(14);clf;hold on; grid on;
          plot(sum(Stats.PHI.^2,2));
          plot(sum(Stats2.PHI.^2,2));
          legend({'invkin_ser', 'invkin3'}, 'interpreter', 'none');
        end
      end
      t1=tic();
      [Q_t_kls, QD_t_kls, QDD_t_kls, Phi1_t_kls,~,~, JointPos_all_kls] = ...
        RP.invkin_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0_traj, s_jj);% Klassen
      calctimes(1)=toc(t1);
      t1=tic();
      [Q_t_tpl, QD_t_tpl, QDD_t_tpl, Phi1_t_tpl,~,~, JointPos_all_tpl] = ...
        RP.invkin2_traj(Traj_0.X, Traj_0.XD, Traj_0.XDD, Traj_0.t, q0_traj, s_jj);% Template
      calctimes(2)=toc(t1);
      trajik_iO = [all(abs(Phi1_t_kls(:)) < 1e-6), all(abs(Phi1_t_tpl(:)) < 1e-6)];
      trajikstat(jj) = trajik_iO(1);
      if any(~trajik_iO)
        if any(jj == [4 5 6 7 8]) % TODO: Das ist noch unbefriedigend. Sollte eigentlich immer gehen.
          Ifirst = find(any(abs(Phi1_t_kls)>1e-6|isnan(Phi1_t_kls),2), 1, 'first');
          warning(['Traj.-IK funktioniert nicht (Abbruch bei %d/%d). Fall ', ...
            '%d allerdings ungünstig und kann zu Singularitäten führen. ', ...
            'Also kein Fehler.'], Ifirst, size(Phi1_t_kls,1), jj);
          continue
        end
        if ~all(~trajik_iO)
          % Bei numerischen Problemen kann es rundungsbedingt vorkommen,
          % dass eine Implementierung funktioniert, und die andere nicht.
          error(['Ergebnis der Traj.-IK zwischen beiden Implementierungen ', ...
            'unterschiedlich. Klasse %d vs Vorlage %d'], trajik_iO(1), trajik_iO(2));
        end
        test_QDD = QDD_t_kls - QDD_t_tpl;
        test_QDD(abs(test_QDD)<1e-10) = 0;
        test_QD = QD_t_kls - QD_t_tpl;
        test_QD(abs(test_QD)<1e-10) = 0;
        test_Q = Q_t_kls - Q_t_tpl;
        test_Q(abs(test_Q)<1e-10) = 0;
        RPstr = ['R', 'P'];
        change_current_figure(13);clf;
        for kkk = 1:RP.NJ
          legnum = find(kkk>=RP.I1J_LEG, 1, 'last');
          legjointnum = kkk-(RP.I1J_LEG(legnum)-1);
          subplot(ceil(sqrt(RP.NJ)), ceil(RP.NJ/ceil(sqrt(RP.NJ))), kkk);
          hold on; grid on;
          hdl1=plot(Q_t_kls(:,kkk), '-');
          hdl2=plot(Q_t_tpl(:,kkk), '--');
          title(sprintf('q %d (%s), L%d,J%d', kkk, RPstr(RP.MDH.sigma(kkk)+1), legnum, legjointnum));
          grid on;
        end
        linkxaxes
        change_current_figure(14);clf;
        for kkk = 1:RP.NJ
          legnum = find(kkk>=RP.I1J_LEG, 1, 'last');
          legjointnum = kkk-(RP.I1J_LEG(legnum)-1);
          subplot(ceil(sqrt(RP.NJ)), ceil(RP.NJ/ceil(sqrt(RP.NJ))), kkk);
          hold on; grid on;
          hdl1=plot(QD_t_kls(:,kkk), '-');
          hdl2=plot(QD_t_tpl(:,kkk), '--');
          title(sprintf('qD %d (%s), L%d,J%d', kkk, RPstr(RP.MDH.sigma(kkk)+1), legnum, legjointnum));
          grid on;
        end
        linkxaxes
        change_current_figure(15);clf;
        for kkk = 1:RP.NJ
          legnum = find(kkk>=RP.I1J_LEG, 1, 'last');
          legjointnum = kkk-(RP.I1J_LEG(legnum)-1);
          subplot(ceil(sqrt(RP.NJ)), ceil(RP.NJ/ceil(sqrt(RP.NJ))), kkk);
          hold on; grid on;
          hdl1=plot(QDD_t_kls(:,kkk), '-');
          hdl2=plot(QDD_t_tpl(:,kkk), '--');
          title(sprintf('qDD %d (%s), L%d,J%d', kkk, RPstr(RP.MDH.sigma(kkk)+1), legnum, legjointnum));
          grid on;
        end
        linkxaxes
        error('Eine Trajektorien-IK funktioniert nicht. Darf nicht sein, da Parameter aus Maßsynthese kommen');
      end
      % Prüfe das Ergebnis der Traj.-IK in sich. Wenn die IK erfolgreich
      % ist, muss die Plattform-Bewegung von allen Beinketten aus gezählt
      % übereinstimmen. Vergleiche nur die Beinketten gegen Soll-EE-Traj.
      % Dadurch auch Fehler der IK erkennbar.
      for iktyp = 1:2
        if ~trajik_iO(iktyp), continue; end % Nicht prüfbar, ob IK stimmt.
        if iktyp == 1
          Q = Q_t_kls; QD = QD_t_kls; QDD = QDD_t_kls; PHI = Phi1_t_kls;
        else
          Q = Q_t_tpl; QD = QD_t_tpl; QDD = QDD_t_tpl; PHI = Phi1_t_tpl;
        end
        Phi3_Traj_red_pre = NaN(size(Phi1_t_kls));
        for j = 1:size(Q,1)
          Phi3_Traj_red_pre(j,:) = RP.constr3(Q(j,:)', Traj_0.X(j,:)');
        end
        if~(all(abs(Phi3_Traj_red_pre(:)) < 1e-6))
          figure(20);clf;
          subplot(2,1,1); hold on; grid on;
          plot(Traj_0.t, PHI);
          subplot(2,1,2); hold on; grid on;
          plot(Traj_0.t, Phi3_Traj_red_pre);
          ylabel('Phi');
          linkxaxes
          error(['Widersprüchliche Bestimmung des Residuums in der ', ...
            'Trajektorie (constr3)']);
        end
        % Plattform-Bewegung neu für 3T2R-Roboter berechnen
        Xsoll = Traj_0.X;
        XDsoll = Traj_0.XD;
        XDDsoll = Traj_0.XDD;
        if ~RP.I_EE_Task(6) && RP.I_EE(6) || ... % Fall der Aufgabenredundanz
            all(RP.I_EE==[1 1 1 1 1 0])  % oder 3T2R-PKM
          [X2,XD2,XDD2] = RP.fkineEE2_traj(Q, QD, QDD);
          Xsoll(:,6) = X2(:,6);
          XDsoll(:,6) = XD2(:,6);
          XDDsoll(:,6) = XDD2(:,6);
        end
        Phi1_Traj_full = NaN(size(Q,1),RP.NLEG*6);
        Phi3_Traj_full = NaN(size(Q,1),RP.NLEG*6);
        Phi3_Traj_red = NaN(size(Phi1_t_kls));
        for j = 1:size(Q,1)
          [~,Phi1_Traj_full(j,:)] = RP.constr1(Q(j,:)', Xsoll(j,:)');
          [Phi3_Traj_red(j,:),Phi3_Traj_full(j,:)] = RP.constr3(Q(j,:)', Xsoll(j,:)');
        end
        assert(all(abs(Phi3_Traj_red(:)) < 1e-6), ['Widersprüchliche ', ...
          'Bestimmung des Residuums in der Trajektorie (constr3)']);
        assert(all(abs(Phi1_Traj_full(:)) < 1e-6), ['Widersprüchliche ', ...
          'Bestimmung des Residuums in der Trajektorie (constr1)']);
        assert(all(abs(Phi3_Traj_full(:)) < 1e-6), ['Widersprüchliche ', ...
          'Bestimmung des Residuums in der Trajektorie (constr3)']);
        for j = 1:RP.NLEG
          [X3,XD3,XDD3] = RP.fkineEE_traj(Q, QD, QDD, j);
          test_X = Xsoll(:,1:6) - X3(:,1:6);
          test_X([false(size(test_X,1),3),abs(abs(test_X(:,4:6))-2*pi)<1e-3]) = 0; % 2pi-Fehler entfernen
          test_XD_abs = XDsoll(:,1:6) - XD3(:,1:6);
          test_XDD_abs = XDDsoll(:,1:6) - XDD3(:,1:6);
          test_XD_rel = test_XD_abs ./ XDsoll(:,1:6);
          test_XDD_rel = test_XDD_abs ./ XDDsoll(:,1:6);
          if max(abs(test_X(:)))>1e-6
            Ifirst = find(any(abs(test_X)>1e-6,2), 1, 'first');
            Icols = find(any(abs(test_X)>1e-6,1),1);
            error(['Die Endeffektor-Trajektorie X aus Beinkette %d stimmt ', ...
              'nicht gegen Soll-Trajektorie. Erstes Vorkommnis: Zeitschritt %d. Max Fehler %1.2e. Komponenten: [%s]'], ...
              j, Ifirst, max(abs(test_X(:))), disp_array(Icols, '%d'));
          end
          if max(abs(test_XD_abs(:)))>1e-6
            Ifirst = find(any(abs(test_XD_abs)>1e-6,2), 1, 'first');
            Icols = find(any(abs(test_XD_abs)>1e-6,1),1);
            error(['Die Endeffektor-Trajektorie XD aus Beinkette %d stimmt ', ...
              'nicht gegen Soll-Trajektorie. Erstes Vorkommnis: Zeitschritt %d. Max Fehler %1.2e. Komponenten: [%s]'], ...
              j, Ifirst, max(abs(test_XD_abs(:))), disp_array(Icols, '%d'));
          end
          % Vereinfachte Beschleunigungsberechung führt zu höherem Fehler.
          % Ausmaß abhängig von Schrittweite und Höhe der Geschwindigkeit.
          if s_jj.simplify_acc, tolscal = 1e3; else, tolscal = 1; end
          % Teste nicht die letzten Zeitschritte, falls die Trajektorie
          % vorzeitig abgebrochen wird. Ursache: Fehler steigt vor Abbruch
          % der Trajektorie bereits an.
          I_DDabserr = abs(test_XDD_abs) > tolscal*1e-3;
          I_DDrelerr = abs(test_XDD_rel) > 1e-3;
          I_firstnan = find(any(isnan(test_XDD_abs),2),1,'first');
          if isempty(I_firstnan)
            I_firstnan = size(test_XDD_abs,1)+1;
          end
          I_DDabserr = I_DDabserr(1:end-5,:);
          I_DDrelerr = I_DDrelerr(1:end-5,:);
          test_XDD_abs(abs(test_XDD_abs)<tolscal*1e-3) = 0; % Zur Übersichtlichkeit
          if any(I_DDrelerr(:) & I_DDabserr(:))
            Ifirst = find(any(I_DDrelerr & I_DDabserr,2), 1, 'first');
            errtxt = sprintf(['Die Endeffektor-Trajektorie XDD aus Beinkette ', ...
              '%d stimmt nicht gegen Soll-Trajektorie. Erstes Vorkommnis: ', ...
              'Zeitschritt %d. Max Fehler %1.2e'], j, Ifirst, max(abs(test_XDD_abs(:))));
            if s_jj.simplify_acc
              warning(['Da vereinfachte Beschleunigung berechnet wurde, ', ...
                'ist das folgende kein (unerwarteter) Fehler:\n%s'], errtxt); %#ok<SPWRN>
            else
              error(errtxt); %#ok<SPERR>
            end
          end
        end
      end
      test_JP = JointPos_all_kls - JointPos_all_tpl;
      test_Q = Q_t_kls - Q_t_tpl;
      if any(abs(test_Q(:))>1e-4)
        Ifirst = find(any(abs(test_Q)>1e-3,2), 1, 'first');
        warning(['Fall %d: invkin_traj vs invkin2_traj: Q1 aus kls und tpl ', ...
          'stimmen nicht überein. Max. Abweichung %1.1e. Erstes Vorkommnis: ', ...
          'Zeitschritt %d'], jj, max(abs(test_Q(:))), Ifirst);
      end
      if any(abs(test_JP(:))>1e-3)
        warning(['Fall %d: invkin_traj vs invkin2_traj: JointPos aus kls ', ...
          'und tpl stimmen nicht überein. Max. Abweichung %1.1e'], jj, ...
          max(abs(test_JP(:))));
      end
      fprintf(['Fall %d: Klassen- und Template-Methode für Traj.-IK stimmen überein.\n', ...
        'Gesamt %1.1fs. Zeiten: invkin_traj: %1.2fs (%1.1fms pro Bahnpunkt),', ...
        'invkin2_traj: %1.2fs (%1.1fms pro Bahnpunkt).\n'], ...
        jj, sum(calctimes(:)), 1000*calctimes(1)/size(Q,1), calctimes(1), ...
        calctimes(2), 1000*calctimes(2)/size(Q,1));
    end
    ResStat.IK_Traj_Erfolg(strcmp(ResStat.Name, PName)) = all(trajikstat(1:3)); % Marker ob alle erfolgreich. Zählt auch NaN als erfolgreich (passt). Nur die betrachten, die funktionieren müssen.
    ResStat.IK_Traj_ErfolgDetail(strcmp(ResStat.Name, PName),:) = trajikstat;
    fprintf('Untersuchung für PKM %s abgeschlossen.\n', PName);
  end % for ii (PKM)
  fprintf('Untersuchung für EE-FG %dT%dR abgeschlossen.\n', sum(EE_FG(1:3)), sum(EE_FG(4:6)));
end % for i_FG
save(fullfile(resdir, 'parroblib_ik_tpl_test.mat'), 'ResStat');
writetable(ResStat, fullfile(resdir, 'parroblib_ik_tpl_test.csv'), 'Delimiter', ';');
fprintf('Gesamtergebnis-Tabelle nach %s gespeichert.\n', resdir);
