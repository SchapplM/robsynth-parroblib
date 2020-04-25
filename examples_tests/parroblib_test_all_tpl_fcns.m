% Rufe alle Funktionen aller in der Bibliothek enthaltenen PKM auf
% Dadurch wird sichergestellt, dass alles vorhanden ist und funktioniert.
% Fokus liegt hier auf den Template-Funktionen und deren Korrektheit

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc
repopath=fileparts(which('serroblib_path_init.m'));

%% Alle PKM durchgehen
% Alle bisher implementierten EE-FG
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

  for ii =  1:length(PNames_Kin) % Debug: find(strcmp(PNames_Kin, 'P6RRPRRR14V3G1P4'));
    PName = [PNames_Kin{ii},'A1']; % Nehme nur die erste Aktuierung (ist egal)
    [~, ~, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    fprintf('Untersuche PKM %s\n', PName);
    
    % Vorlagen-Funktionen aktualisieren
    parroblib_create_template_functions({PNames_Kin{ii}}, false, true);
    %% Roboter Initialisieren
    R = parroblib_create_robot_class(PName, 1, 0.3);
    % Initialisierung der Funktionen: Kompilierte Funktionen nehmen
    R.fill_fcn_handles(false, false); % true, true
    
    %% Kinematikparameter zufällig setzen (ist egal)
    pkin_rand = rand(length(R.Leg(1).pkin),1); % Symmetrische PKM mit zufälligen Parametern
    for il = 1:R.NLEG
      R.Leg(il).update_mdh(pkin_rand);
      R.Leg(il).qlim = repmat([-pi,pi], R.Leg(il).NJ,1);
    end
    R.update_base(rand(3,1), rand(3,1));
    R.update_EE(rand(3,1), rand(3,1));
    %% Test-Szenario definieren
    qlim_pkm = cat(1,R.Leg(:).qlim);
    q = qlim_pkm(:,1) + (qlim_pkm(:,2)-qlim_pkm(:,1)).*rand(R.NJ,1);
    
    n = 20;
    Q_test = qlim_pkm(:,1) + (qlim_pkm(:,2)-qlim_pkm(:,1)).*rand(R.NJ,n);
    QD_test = rand(R.NJ,n);
    X_test = rand(6, n);
    XD_test = rand(6, n);
    X_test(:,~R.I_EE) = 0;
    XD_test(:,~R.I_EE) = 0;
    %% Rufe Funktionen auf
    % Struktur aus ParRob/invkin2_traj
    Leg_I_EE_Task = true(R.NLEG,6);
    Leg_pkin_gen = zeros(R.NLEG,length(R.Leg(1).pkin_gen));
    Leg_T_N_E_vec = zeros(6,R.NLEG);% 1:3 eularwinkel 4:6 Position
    Leg_T_0_W_vec = zeros(6,R.NLEG);% 1:3 eularwinkel 4:6 Position
    Leg_I_EElink = zeros(R.NLEG,1);
    Leg_phi_W_0 = zeros(3,R.NLEG);
    Leg_phiconv_W_0 = uint8(zeros(R.NLEG,1));
    Leg_NQJ = zeros(R.NLEG,1);
    phiconv_W_E = uint8(R.phiconv_W_E);
    Leg_sigmaJ = zeros(R.Leg(1).NJ,R.NLEG);
    Leg_qlim = zeros(6,2*R.NLEG);
    Leg_phiconv_W_E = uint8(zeros(R.NLEG,1));
    for i = 1:R.NLEG
      Leg_I_EE_Task(i,:) = R.Leg(i).I_EE_Task;
      Leg_pkin_gen(i,:) = R.Leg(i).pkin_gen;
      Leg_I_EElink(i,:) = R.Leg(i).I_EElink;
      T_N_E = R.Leg(i).T_N_E;
      Leg_T_N_E_vec(1:3,i) = r2eulxyz(T_N_E(1:3,1:3));
      Leg_T_N_E_vec(4:6,i) = T_N_E(1:3,4);
      T_0_W = R.Leg(i).T_0_W;
      Leg_T_0_W_vec(1:3,i) = r2eulxyz(T_0_W(1:3,1:3));
      Leg_T_0_W_vec(4:6,i) = T_0_W(1:3,4);
      Leg_phi_W_0(:,i) = R.Leg(i).phi_W_0;
      Leg_phiconv_W_0(i) = R.Leg(i).phiconv_W_0;
      Leg_NQJ(i) = R.Leg(i).NJ;
      Leg_sigmaJ(:,i) = R.Leg(i).MDH.sigma(R.Leg(i).MDH.mu>=1);
      Leg_qlim(1:R.Leg(i).NJ,(1+2*(i-1)):(2+2*(i-1))) = R.Leg(i).qlim;
      Leg_phiconv_W_E(i) = R.Leg(i).phiconv_W_E;
    end
    % Struktur 1 aus pkm_invkin_traj
    s = struct('I_EE', R.I_EE,...
      'I_EE_Task', R.I_EE_Task,...
      'simplify_acc', false,...
      'mode_IK', 2,...
      'debug', false,...
      'I_constr_t_red', R.I_constr_t_red,...
      'I_constr_r_red', R.I_constr_r_red,...
      'I1constr_red', R.I1constr_red,...
      'I2constr_red', R.I2constr_red,...
      'I_constr_t', R.I_constr_t,...
      'I_constr_r', R.I_constr_r,...
      'I_qa', R.I_qa,...
      'r_P_B_all', R.r_P_B_all, ...
      'phi_P_B_all', R.phi_P_B_all, ...
      'NLEG', R.NLEG, ...
      'NJ', R.NJ, ...
      'phiconv_W_E', phiconv_W_E, ...
      'T_P_E', R.T_P_E, ...
      'I1J_LEG', R.I1J_LEG, ...
      'I2J_LEG', R.I2J_LEG, ...
      'Leg_I_EE_Task', Leg_I_EE_Task, ...
      'Leg_pkin_gen', Leg_pkin_gen, ...
      'Leg_T_N_E_vec', Leg_T_N_E_vec, ...
      'Leg_T_0_W_vec', Leg_T_0_W_vec, ...
      'Leg_I_EElink', Leg_I_EElink, ...
      'Leg_phi_W_0', Leg_phi_W_0, ...
      'Leg_phiconv_W_0', Leg_phiconv_W_0, ...
      'Leg_NQJ', Leg_NQJ,...
      'Leg_sigmaJ', Leg_sigmaJ, ...
      'Leg_qlim', Leg_qlim, ...
      'Leg_phiconv_W_E', Leg_phiconv_W_E);
    K_def = 0.5*ones(R.Leg(1).NQJ,1);
    % Struktur 2 aus pkm_invkin_traj
    s_ser = struct('reci', true, ...
      'K', K_def, ... % Verstärkung
      'Kn', 1e-2*ones(R.Leg(1).NQJ,1), ... % Verstärkung
      'wn', zeros(2,1), ... % Gewichtung der Nebenbedingung
      'scale_lim', 0.0, ... % Herunterskalierung bei Grenzüberschreitung
      'maxrelstep', 0.05, ... % Maximale auf Grenzen bezogene Schrittweite
      'normalize', true, ... % Normalisieren auf +/- 180°
      'n_min', 0, ... % Minimale Anzahl Iterationen
      'n_max', 1000, ... % Maximale Anzahl Iterationen
      'rng_seed', NaN, ... Initialwert für Zufallszahlengenerierung
      'Phit_tol', 1e-8, ... % Toleranz für translatorischen Fehler
      'Phir_tol', 1e-8, ... % Toleranz für rotatorischen Fehler
      'retry_limit', 100);
    % Struktur aus constr4grad_q und constr4gradD_q
    s_q = struct(   'I_EE', s.I_EE,...
      'I_constr_t_red', s.I_constr_t_red,...
      'I_constr_r_red', s.I_constr_r_red,...     
      'I_constr_t', s.I_constr_t,...
      'I_constr_r', s.I_constr_r,...
      'NLEG', s.NLEG,...
      'NJ', s.NJ,...
      'I1J_LEG', s.I1J_LEG,...
      'I2J_LEG', s.I2J_LEG,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EElink', s.Leg_I_EElink,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0,...
      'Leg_NQJ', s.Leg_NQJ);
    % Struktur aus constr4grad_x und constr4gradD_x
    s_x = struct(   'I_EE', s.I_EE,...
      'I_constr_t_red', s.I_constr_t_red,...
      'I_constr_r_red', s.I_constr_r_red,...
      'I_constr_t', s.I_constr_t,...
      'I_constr_r', s.I_constr_r,...
      'NLEG', s.NLEG,...
      'r_P_B_all', s.r_P_B_all,...
      'phiconv_W_E', s.phiconv_W_E,...
      'T_P_E', s.T_P_E,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task);
    for mextest = 1:2 % mextest = 2: mex datei aufrufen
      for jj = 1:n
        q = Q_test(:,jj);
        qd = QD_test(:,jj);
        x = X_test(:,jj);
        xd = XD_test(:,jj);
        if mextest == 1
          eval(sprintf('Phi4_dq_file=%s_constr4grad_q(q, s_q);', PName_Legs));
        else
          eval(sprintf('Phi4_dq_file=%s_constr4grad_q_mex(q, s_q);', PName_Legs));
        end
        Phi4_dq_class = R.constr4grad_q(q);
        test_Phi4dq = Phi4_dq_class - Phi4_dq_file;
        if any(abs(test_Phi4dq(:)) > 1e-10)
          error('Berechnung von constr4grad_q stimmt nicht zwischen Template-Funktion und Klasse');
        end

        if mextest == 1
          eval(sprintf('Phi4_dqD_file=%s_constr4gradD_q(q, qd, s_q);', PName_Legs));
        else
          eval(sprintf('Phi4_dqD_file=%s_constr4gradD_q_mex(q, qd, s_q);', PName_Legs));
        end
        Phi4_dqD_class = R.constr4gradD_q(q, qd);
        test_Phi4dqD = Phi4_dqD_class - Phi4_dqD_file;
        if any(abs(test_Phi4dqD(:)) > 1e-10)
          error('Berechnung von constr4gradD_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        
        if mextest == 1
          eval(sprintf('Phi4_dx_file=%s_constr4grad_x(x, s_x);', PName_Legs));
        else
          eval(sprintf('Phi4_dx_file=%s_constr4grad_x_mex(x, s_x);', PName_Legs));
        end
        Phi4_dx_class = R.constr4grad_x(x);
        test_Phi4dx = Phi4_dx_class - Phi4_dx_file;
        if any(abs(test_Phi4dx(:)) > 1e-10)
          error('Berechnung von constr4grad_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
        
        if mextest == 1
          eval(sprintf('Phi4_dxd_file=%s_constr4gradD_x(x, xd, s_x);', PName_Legs));
        else
          eval(sprintf('Phi4_dxd_file=%s_constr4gradD_x_mex(x, xd, s_x);', PName_Legs));
        end
        Phi4_dxd_class = R.constr4gradD_x(x, xd);
        test_Phi4dxd = Phi4_dxd_class - Phi4_dxd_file;
        if any(abs(test_Phi4dxd(:)) > 1e-10)
          error('Berechnung von constr4gradD_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
      end    
    end
    fprintf('PKM %s erfolgreich getestet\n', PName);
  end
end
