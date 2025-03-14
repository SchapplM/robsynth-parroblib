% Rufe alle Funktionen aller in der Bibliothek enthaltenen PKM auf
% Dadurch wird sichergestellt, dass alles vorhanden ist und funktioniert.
% Fokus liegt hier auf den Template-Funktionen und deren Korrektheit
% (es erfolgt keine inhaltliche Prüfung der Ergebnisse).
% Geprüfte Szenarien:
% * kompilierte und nicht kompilierte Version
% * Kinematik mit und ohne Aufgabenredundanz (nicht definierte z-Drehung am EE)

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc
repopath=fileparts(which('serroblib_path_init.m'));

%% Initialisierung
% Schalter zur Neugenerierung der Template-Funktionen
tpl_fcn_neu = true;
test_mex = true;
test_taskred = true;
recompile_mex = true; % nur notwendig, falls Template-Dateien geändert wurden.
max_num_pkm = 2; % Reduziere die Anzahl der geprüften PKM pro FG
shuffle_pkm_selection = true; % Zufällige Auswahl der PKM. Zum Debuggen deaktivieren.

%% Alle PKM durchgehen
% Alle bisher implementierten EE-FG
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 0; ...
            1 1 1 1 1 1];

for i_FG = 1:size(EEFG_Ges,1)
  EE_FG = EEFG_Ges(i_FG,:);
  [PNames_Kin, ~] = parroblib_filter_robots(EE_FG, 6);
  if isempty(PNames_Kin)
    continue % Es gibt keine PKM mit diesen FG.
  end
  III = 1:length(PNames_Kin);
  if shuffle_pkm_selection
    III = III(randperm(length(III)));
  end
  % Debug:
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
  for ii =  III(1:min(max_num_pkm, length(III))) % Debug: find(strcmp(PNames_Kin, 'P6RRPRRR14V3G1P4'));
    PName = [PNames_Kin{ii},'A1']; % Nehme nur die erste Aktuierung (ist egal)
    [~, LEG_Names, ~, ~, ~, ~, ~, ~, PName_Legs, ~] = parroblib_load_robot(PName);
    fprintf('Untersuche PKM %s\n', PName);
    
    % Vorlagen-Funktionen aktualisieren. Die Funktionen müssen auf jeden
    % Fall kompiliert werden, da das unten getestet wird.
    serroblib_create_template_functions({LEG_Names{1}},  ~tpl_fcn_neu,recompile_mex); %#ok<CCAT1>
    parroblib_create_template_functions({PNames_Kin{ii}},~tpl_fcn_neu,recompile_mex); %#ok<CCAT1>

    %% Roboter Initialisieren
    R = parroblib_create_robot_class(PName, '', 1, 0.3);
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
    n = 20;
    Q_test = qlim_pkm(:,1)' + (qlim_pkm(:,2)-qlim_pkm(:,1))'.*rand(n, R.NJ);
    QD_test = rand(n, R.NJ);
    QDD_test = rand(n, R.NJ);
    X_test = rand(n, 6);
    XD_test = rand(n, 6);
    X_test(:,~R.I_EE) = 0;
    XD_test(:,~R.I_EE) = 0;
    for trtest = 1:(1+test_taskred) % trtest = 2: mit Aufgabenredundanz
    if trtest == 2
      if all(EE_FG == [1 1 1 1 1 1])
        R.I_EE_Task = logical([1 1 1 1 1 0]);
      elseif all(EE_FG == [1 1 0 0 0 1])
        R.I_EE_Task = logical([1 1 0 0 0 0]);
      elseif all(EE_FG == [1 1 1 0 0 1])
        R.I_EE_Task = logical([1 1 1 0 0 0]);
      else
        % Aufgabenredundanz für diese Struktur-FG noch nicht definiert.
        continue
      end
      R.update_EE_FG(R.I_EE, R.I_EE_Task);
    end
    %% Initialisierung der Parameter-Strukturen
    % Struktur aus ParRob/invkin2_traj
    Leg_I_EE_Task = true(R.NLEG,6);
    Leg_pkin_gen = zeros(R.NLEG,length(R.Leg(1).pkin_gen));
    Leg_T_N_E_vec = zeros(6,R.NLEG);% 1:3 Euler-Winkel 4:6 Position
    Leg_T_0_W_vec = zeros(6,R.NLEG);% 1:3 Euler-Winkel 4:6 Position
    Leg_phi_W_0 = zeros(3,R.NLEG);
    Leg_phiconv_W_0 = uint8(zeros(R.NLEG,1));
    phiconv_W_E = uint8(R.phiconv_W_E);
    Leg_sigmaJ = zeros(R.Leg(1).NJ,R.NLEG);
    Leg_qlim = zeros(6,2*R.NLEG);
    Leg_phiconv_W_E = uint8(zeros(R.NLEG,1));
    for i = 1:R.NLEG
      Leg_I_EE_Task(i,:) = R.Leg(i).I_EE_Task;
      Leg_pkin_gen(i,:) = R.Leg(i).pkin_gen;
      T_N_E = R.Leg(i).T_N_E;
      Leg_T_N_E_vec(1:3,i) = r2eulxyz(T_N_E(1:3,1:3));
      Leg_T_N_E_vec(4:6,i) = T_N_E(1:3,4);
      T_0_W = R.Leg(i).T_0_W;
      Leg_T_0_W_vec(1:3,i) = r2eulxyz(T_0_W(1:3,1:3));
      Leg_T_0_W_vec(4:6,i) = T_0_W(1:3,4);
      Leg_phi_W_0(:,i) = R.Leg(i).phi_W_0;
      Leg_phiconv_W_0(i) = R.Leg(i).phiconv_W_0;
      Leg_sigmaJ(:,i) = R.Leg(i).MDH.sigma(R.Leg(i).MDH.mu>=1);
      Leg_qlim(1:R.Leg(i).NJ,(1+2*(i-1)):(2+2*(i-1))) = R.Leg(i).qlim;
      Leg_phiconv_W_E(i) = R.Leg(i).phiconv_W_E;
    end
    % Struktur 1 aus pkm_invkin_traj
    s = struct( ...
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
      'phiconv_W_E', phiconv_W_E, ...
      'T_P_E', R.T_P_E, ...
      'Leg_I_EE_Task', Leg_I_EE_Task, ...
      'Leg_pkin_gen', Leg_pkin_gen, ...
      'Leg_T_N_E_vec', Leg_T_N_E_vec, ...
      'Leg_T_0_W_vec', Leg_T_0_W_vec, ...
      'Leg_phi_W_0', Leg_phi_W_0, ...
      'Leg_phiconv_W_0', Leg_phiconv_W_0, ...
      'Leg_sigmaJ', Leg_sigmaJ, ...
      'Leg_qlim', Leg_qlim, ...
      'Leg_phiconv_W_E', Leg_phiconv_W_E);
    % Struktur 2 aus pkm_invkin_traj
    s_ser = struct( ...
      'reci', true, ...
      'K', ones(R.Leg(1).NQJ,1), ... % Verstärkung
      'Kn', ones(R.Leg(1).NQJ,1), ... % Verstärkung
      'wn', zeros(3,1), ... % Gewichtung der Nebenbedingung
      'scale_lim', 0.0, ... % Herunterskalierung bei Grenzüberschreitung
      'maxrelstep', 0.05, ... % Maximale auf Grenzen bezogene Schrittweite
      'normalize', true, ... % Normalisieren auf +/- 180°
      'n_min', 0, ... % Minimale Anzahl Iterationen
      'n_max', 1000, ... % Maximale Anzahl Iterationen
      'rng_seed', NaN, ... Initialwert für Zufallszahlengenerierung
      'Phit_tol', 1e-8, ... % Toleranz für translatorischen Fehler
      'Phir_tol', 1e-8, ... % Toleranz für rotatorischen Fehler
      'retry_limit', 100);
    % Struktur aus fkineEE_traj
    s_fkine = struct( ...
      'T_P_E', s.T_P_E, ...
      'r_P_B_all', s.r_P_B_all,...
      'phi_P_B_all', R.phi_P_B_all,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_T_0_W_vec', Leg_T_0_W_vec, ...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus fkine_coll
    s_fkine_coll = struct( ...
      'r_P_B_all', s.r_P_B_all,...
      'phi_P_B_all', R.phi_P_B_all,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'T_P_E', s.T_P_E, ...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_T_0_W_vec', Leg_T_0_W_vec);
    % Struktur aus constr1_trans
    s_ct_1 = struct( ...
      'r_P_B_all', s.r_P_B_all,...
      'phi_P_B_all', R.phi_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_T_0_W_vec', Leg_T_0_W_vec, ...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr1grad_tq
    s_q_1 = struct( ...
      'I_constr_t_red', s.I_constr_t_red,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr1D_trans
    s_qD_1 = struct( ...
      'r_P_B_all', s.r_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr4D_rot
    s_qD_4r = struct( ...
      'I_constr_r_red', s.I_constr_r_red, ...
      'r_P_B_all', s.r_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr4D_traj
    s_qD_4 = struct( ...
      'I_constr_t', s.I_constr_t, ...
      'I_constr_r', s.I_constr_r, ...
      'I_constr_r_red', s.I_constr_r_red,...
      'r_P_B_all', s.r_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr4grad_q und constr4gradD_q
    s_q = struct( ...
      'I_constr_t_red', s.I_constr_t_red,...
      'I_constr_r_red', s.I_constr_r_red,...     
      'I_constr_t', s.I_constr_t,...
      'I_constr_r', s.I_constr_r,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr2grad_q
    s_q_2 = struct( ...
      'r_P_B_all', s.r_P_B_all,...
      'phi_P_B_all', R.phi_P_B_all,...
      'phiconv_W_E', s.phiconv_W_E, ...
      'T_P_E', s.T_P_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr2gradD_q
    s_qD_2 = struct( ...
      'I_EE_Task', R.I_EE_Task,...
      'I_constr_t_red', s.I_constr_t_red,...
      'phi_P_B_all', R.phi_P_B_all,...
      'r_P_B_all', s.r_P_B_all,...
      'phiconv_W_E', s.phiconv_W_E, ...
      'T_P_E', s.T_P_E, ...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr2grad_x
    s_x_2 = struct( ...
      'I_EE_Task', R.I_EE_Task,...
      'phi_P_B_all', R.phi_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'phiconv_W_E', s.phiconv_W_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_T_0_W_vec', Leg_T_0_W_vec, ...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr2gradD_x
    s_xD_2 = struct( ...
      'I_EE_Task', R.I_EE_Task,...
      'I_constr_t_red', s.I_constr_t_red,...
      'phi_P_B_all', R.phi_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'phiconv_W_E', s.phiconv_W_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr2grad_rr
    s_x_4 = rmfield(s_x_2, 'I_EE_Task');
    % Struktur aus constr2gradD_rr
    s_xD_4 = rmfield(s_xD_2, {'I_EE_Task', 'I_constr_t_red'});
    % Struktur aus constr3
    s3 = struct( ...
      'I_constr_t_red', s.I_constr_t_red,...
      'I_constr_r_red', s.I_constr_r_red,...
      'r_P_B_all', s.r_P_B_all,...
      'phi_P_B_all', R.phi_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'phiconv_W_E', s.phiconv_W_E, ...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_T_0_W_vec', s.Leg_T_0_W_vec,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr3grad_q
    s_q_3 = struct( ...
      'I_EE_Task', R.I_EE_Task,...
      'r_P_B_all', s.r_P_B_all,...
      'phi_P_B_all', R.phi_P_B_all,...
      'phiconv_W_E', s.phiconv_W_E, ...
      'T_P_E', s.T_P_E, ...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr3gradD_q
    s_qD_3 = struct( ...
      'I_EE_Task', R.I_EE_Task,...
      'I_constr_t_red', s.I_constr_t_red,...
      'phi_P_B_all', R.phi_P_B_all,...
      'r_P_B_all', s.r_P_B_all,...
      'phiconv_W_E', s.phiconv_W_E, ...
      'T_P_E', s.T_P_E, ...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr4grad_x und constr4gradD_x
    s_x = struct( ...
      'I_constr_t_red', s.I_constr_t_red,...
      'I_constr_r_red', s.I_constr_r_red,...
      'I_constr_t', s.I_constr_t,...
      'I_constr_r', s.I_constr_r,...
      'r_P_B_all', s.r_P_B_all,...
      'phiconv_W_E', s.phiconv_W_E,...
      'T_P_E', s.T_P_E,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task);
    % Struktur aus constr3grad_x
    s_x_3 = struct( ...
      'I_constr_red', R.I_constr_red, ...
      'phi_P_B_all', R.phi_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'phiconv_W_E', s.phiconv_W_E,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_T_0_W_vec', Leg_T_0_W_vec, ...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    % Struktur aus constr3gradD_x
    s_xD_3 = struct( ...
      'I_constr_red', R.I_constr_red,...
      'I_constr_t_red', s.I_constr_t_red,...
      'phi_P_B_all', R.phi_P_B_all,...
      'T_P_E', s.T_P_E, ...
      'phiconv_W_E', s.phiconv_W_E,...
      'Leg_pkin_gen', s.Leg_pkin_gen,...
      'Leg_T_N_E_vec', s.Leg_T_N_E_vec,...
      'Leg_I_EE_Task', s.Leg_I_EE_Task,...
      'Leg_phi_W_0', s.Leg_phi_W_0,...
      'Leg_phiconv_W_0', s.Leg_phiconv_W_0);
    %% Rufe Funktionen auf und vergleiche Implementierungen
    for mextest = 1:(1+test_mex) % mextest = 2: mex-Datei aufrufen
      if mextest == 2
        % Benutze auch in der Matlab-Klasse kompilierte Funktionen
        R.fill_fcn_handles(true, true);
      end
      % Teste Trajektorien-Funktionen
      %% Teste fkineEE_traj
      for ll = 1:R.NLEG
        if mextest == 1
          eval(sprintf('[X_file, XD_file, XDD_file]=%s_fkineEE_traj(Q_test, QD_test, QDD_test, uint8(ll), s_fkine);', PName_Legs));
        else
          % if ll == 1, matlabfcn2mex({sprintf('%s_fkineEE_traj',PName_Legs)}); end
          eval(sprintf('[X_file, XD_file, XDD_file]=%s_fkineEE_traj_mex(Q_test, QD_test, QDD_test, uint8(ll), s_fkine);', PName_Legs));
        end
        [X_class, XD_class, XDD_class] = R.fkineEE_traj(Q_test, QD_test, QDD_test, uint8(ll));
        test_X = X_class - X_file;
        test_XD = XD_class - XD_file;
        test_XDD = XDD_class - XDD_file;
        if any(abs([test_X(:);test_XD(:);test_XDD(:)]) > 1e-8) || ...
           any(isnan([test_X(:);test_XD(:);test_XDD(:)]))
          error('Berechnung von fkineEE_traj stimmt nicht zwischen Template-Funktion und Klasse');
        end
        [X_class2, XD_class2, XDD_class2] = R.fkineEE2_traj(Q_test, QD_test, QDD_test, uint8(ll));
        test_X2 = X_class - X_class2;
        test_XD2 = XD_class - XD_class2;
        test_XDD2 = XDD_class - XDD_class2;
        if any(abs([test_X2(:);test_XD2(:);test_XDD2(:)]) > 1e-8) || ...
           any(isnan([test_X2(:);test_XD2(:);test_XDD2(:)]))
          error('Berechnung von fkineEE_traj stimmt nicht als Aufruf der Template-Funktion in der Klasse');
        end
      end
      %% Teste fkine_coll
      for jj = 1:n
        q = Q_test(jj,:)';
        if mextest == 1
          eval(sprintf('[Tc_file, JP_file]=%s_fkine_coll(q, s_fkine_coll);', PName_Legs));
        else
          % if jj == 1, matlabfcn2mex({sprintf('%s_fkine_coll', PName_Legs)}); end
          eval(sprintf('[Tc_file, JP_file]=%s_fkine_coll_mex(q, s_fkine_coll);', PName_Legs));
        end
        [Tc_class, JP_class] = R.fkine_coll(q);
        test_Tc = Tc_class - Tc_file;
        test_JP = JP_class - JP_file;
        if any(abs([test_Tc(:);test_JP(:)]) > 1e-8) || ...
           any(isnan([test_Tc(:);test_JP(:)]))
          error('Berechnung von fkine_coll stimmt nicht zwischen Template-Funktion und Klasse');
        end
      end
      %% Teste constr4D_traj
      if mextest == 1
        eval(sprintf('PHID_file=%s_constr4D_traj(Q_test, QD_test, X_test, XD_test, s_qD_4);', PName_Legs));
      else
        % matlabfcn2mex({sprintf('%s_constr4D_traj',PName_Legs)});
        eval(sprintf('PHID_file=%s_constr4D_traj_mex(Q_test, QD_test, X_test, XD_test, s_qD_4);', PName_Legs));
      end
      PHID_class = NaN(n,6*R.NLEG);
      for i = 1:n
        [~,Phi_i] = R.constr4D(Q_test(i,:)', QD_test(i,:)', X_test(i,:)', XD_test(i,:)');
        PHID_class(i,:) = Phi_i;
      end
      test_PHID = PHID_class - PHID_file;
      if any(abs(test_PHID(:)) > 1e-8) || any(isnan(test_PHID(:)))
        error('Berechnung von constr4D_traj stimmt nicht zwischen Template-Funktion und Klasse');
      end
      PHID_class2 = R.constr4D2_traj(Q_test, QD_test, X_test, XD_test);
      test_PHID2 = PHID_class2 - PHID_class;
      if any(abs(test_PHID2(:)) > 1e-8) || any(isnan(test_PHID2(:)))
        error('Berechnung von constr4D_traj stimmt nicht als Aufruf der Template-Funktion in der Klasse');
      end
      %% Teste Einzelpunkt-Funktionen
      for jj = 1:n
        q = Q_test(jj,:)';
        qD = QD_test(jj,:)';
        x = X_test(jj,:)';
        xD = XD_test(jj,:)';

        %% Teste constr1_trans
        if mextest == 1
          eval(sprintf('[Phi1t_file,Phi1t_file_full]=%s_constr1_trans(q, x, s_ct_1);', PName_Legs));
        else
          % if jj == 1, matlabfcn2mex({sprintf('%s_constr1_trans', PName_Legs)}); end
          eval(sprintf('[Phi1t_file,Phi1t_file_full]=%s_constr1_trans_mex(q, x, s_ct_1);', PName_Legs));
        end
        [Phi1t_class,Phi1t_class_full] = R.constr1_trans(q, x);
        test_Phi1t = Phi1t_class - Phi1t_file;
        test_Phi1t_full = Phi1t_class_full - Phi1t_file_full;
        if any(abs([test_Phi1t(:); test_Phi1t_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi1t(:); test_Phi1t_full(:)]))
          error('Berechnung von constr1_trans stimmt nicht zwischen Template-Funktion und Klasse');
        end
        continue
        %% Teste constr1D_trans
        if mextest == 1
          eval(sprintf('[Phi1Dt_file, Phi1Dt_file_full]=%s_constr1D_trans(q, qD, x, xD, s_qD_1);', PName_Legs));
        else
          % if jj == 1, matlabfcn2mex({sprintf('%s_constr1D_trans',PName_Legs)}); end
          eval(sprintf('[Phi1Dt_file, Phi1Dt_file_full]=%s_constr1D_trans_mex(q, qD, x, xD, s_qD_1);', PName_Legs));
        end
        [Phi1Dt_class, Phi1Dt_class_full] = R.constr1D_trans(q, qD, x, xD);
        test_Phi1Dt = Phi1Dt_class - Phi1Dt_file;
        test_Phi1Dt_full = Phi1Dt_class_full - Phi1Dt_file_full;
        if any(abs([test_Phi1Dt(:); test_Phi1Dt_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi1Dt(:); test_Phi1Dt_full(:)]))
          error('Berechnung von constr1D_trans stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr4D_rot
        if mextest == 1
          eval(sprintf('[Phi4Dr_file,Phi4Dr_file_full]=%s_constr4D_rot(q, qD, x, xD, s_qD_4r);', PName_Legs));
        else
          % if jj == 1, matlabfcn2mex({sprintf('%s_constr4D_rot',PName_Legs)}); end
          eval(sprintf('[Phi4Dr_file,Phi4Dr_file_full]=%s_constr4D_rot_mex(q, qD, x, xD, s_qD_4r);', PName_Legs));
        end
        [Phi4Dr_class, Phi4Dr_class_full] = R.constr4D_rot(q, qD, x, xD);
        if i_FG == 4 || trtest == 2 % nicht für 3T2R definiert
          test_Phi4Dr = 0; % für 3T2R Platzhalter für reduzierten Ausdruck
        else % 
          test_Phi4Dr = Phi4Dr_class - Phi4Dr_file;
        end
        test_Phi4Dr_full = Phi4Dr_class_full - Phi4Dr_file_full;
        if any(abs([test_Phi4Dr(:); test_Phi4Dr_full(:)]) > 1e-10)
          error('Berechnung von constr4D_rot stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr1grad_tq
        if mextest == 1
          eval(sprintf('[Phi1t_dq_file, Phi1t_dq_file_full]=%s_constr1grad_tq(q, s_q_1);', PName_Legs));
        else
          eval(sprintf('[Phi1t_dq_file, Phi1t_dq_file_full]=%s_constr1grad_tq_mex(q, s_q_1);', PName_Legs));
        end
        [Phi1t_dq_class, Phi1t_dq_class_full] = R.constr1grad_tq(q);
        if i_FG == 4 || trtest == 2 % nicht für 3T2R definiert
          test_Phi1tdq = 0;
        else
          test_Phi1tdq = Phi1t_dq_class - Phi1t_dq_file;
        end
        test_Phi1tdq_full = Phi1t_dq_class_full - Phi1t_dq_file_full;
        if any(abs([test_Phi1tdq(:);test_Phi1tdq_full(:)]) > 1e-10) || ...
           any(isnan([test_Phi1tdq(:);test_Phi1tdq_full(:)]))
          error('Berechnung von constr1grad_tq stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr4grad_q
        if mextest == 1
          eval(sprintf('[Phi4_dq_file,Phi4_dq_file_full]=%s_constr4grad_q(q, s_q);', PName_Legs));
        else
          eval(sprintf('[Phi4_dq_file,Phi4_dq_file_full]=%s_constr4grad_q_mex(q, s_q);', PName_Legs));
        end
        [Phi4_dq_class, Phi4_dq_class_full] = R.constr4grad_q(q);
        if i_FG == 4 || trtest == 2 % nicht für 3T2R definiert
          test_Phi4dq = 0; % Reduzierung von Zeilen nicht sinnvoll. Platzhalter.
        else
          test_Phi4dq = Phi4_dq_class - Phi4_dq_file;
        end
        test_Phi4dq_full = Phi4_dq_class_full - Phi4_dq_file_full;
        if any(abs([test_Phi4dq(:);test_Phi4dq_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi4dq(:);test_Phi4dq_full(:)]))
          error('Berechnung von constr4grad_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr4gradD_q
        if mextest == 1
          eval(sprintf('[Phi4_dqD_file,Phi4_dqD_file_full]=%s_constr4gradD_q(q, qD, s_q);', PName_Legs));
        else
          eval(sprintf('[Phi4_dqD_file,Phi4_dqD_file_full]=%s_constr4gradD_q_mex(q, qD, s_q);', PName_Legs));
        end
        [Phi4_dqD_class,Phi4_dqD_class_full] = R.constr4gradD_q(q, qD);
        if i_FG == 4 || trtest == 2 % nicht für 3T2R definiert
          test_Phi4dqD = 0;
        else
          test_Phi4dqD = Phi4_dqD_class - Phi4_dqD_file;
        end
        test_Phi4dqD_full = Phi4_dqD_class_full - Phi4_dqD_file_full;
        if any(abs([test_Phi4dqD(:);test_Phi4dqD_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi4dqD(:);test_Phi4dqD_full(:)]))
          error('Berechnung von constr4gradD_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr4grad_x
        if mextest == 1
          eval(sprintf('[Phi4_dx_file,Phi4_dx_file_full]=%s_constr4grad_x(x, s_x);', PName_Legs));
        else
          eval(sprintf('[Phi4_dx_file,Phi4_dx_file_full]=%s_constr4grad_x_mex(x, s_x);', PName_Legs));
        end
        [Phi4_dx_class,Phi4_dx_class_full] = R.constr4grad_x(x);
        if i_FG == 4 || trtest == 2 % nicht für 3T2R definiert
          test_Phi4dx = 0;
        else
          test_Phi4dx = Phi4_dx_class - Phi4_dx_file;
        end
        test_Phi4dx_full = Phi4_dx_class_full - Phi4_dx_file_full;
        if any(abs([test_Phi4dx(:);test_Phi4dx_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi4dx(:);test_Phi4dx_full(:)]))
          error('Berechnung von constr4grad_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr4gradD_x
        if mextest == 1
          eval(sprintf('[Phi4_dxD_file,Phi4_dxD_file_full]=%s_constr4gradD_x(x, xD, s_x);', PName_Legs));
        else
          eval(sprintf('[Phi4_dxD_file,Phi4_dxD_file_full]=%s_constr4gradD_x_mex(x, xD, s_x);', PName_Legs));
        end
        [Phi4_dxD_class,Phi4_dxD_class_full] = R.constr4gradD_x(x, xD);
        if i_FG == 4 || trtest == 2 % nicht für 3T2R definiert
          test_Phi4dxD = 0;
        else
          test_Phi4dxD = Phi4_dxD_class - Phi4_dxD_file;
        end
        test_Phi4dxD_full = Phi4_dxD_class_full - Phi4_dxD_file_full;
        if any(abs([test_Phi4dxD(:);test_Phi4dxD(:)]) > 1e-10) || ...
            any(isnan([test_Phi4dxD(:);test_Phi4dxD(:)]))
          error('Berechnung von constr4gradD_x stimmt nicht zwischen Template-Funktion und Klasse');
        end

        %% Teste constr3
        if mextest == 1
          eval(sprintf('[Phi3_file,Phi3_file_full]=%s_constr3(q, x, s3);', PName_Legs));
        else
          eval(sprintf('[Phi3_file,Phi3_file_full]=%s_constr3_mex(q, x, s3);', PName_Legs));
        end
        [Phi3_class,Phi3_class_full] = R.constr3(q, x);
        test_Phi3 = Phi3_class - Phi3_file;
        test_Phi3_full = Phi3_class_full - Phi3_file_full;
        if any(abs([test_Phi3(:); test_Phi3_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi3(:); test_Phi3_full(:)]))
          error('Berechnung von constr3 stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr3grad_q
        if mextest == 1
          eval(sprintf('[Phi3_dq_file,Phi3_dq_file_full]=%s_constr3grad_q(q, x, s_q_3);', PName_Legs));
        else
          eval(sprintf('[Phi3_dq_file,Phi3_dq_file_full]=%s_constr3grad_q_mex(q, x, s_q_3);', PName_Legs));
        end
        [Phi3_dq_class,Phi3_dq_class_full] = R.constr3grad_q(q, x);
        test_Phi3dq = Phi3_dq_class - Phi3_dq_file;
        test_Phi3dq_full = Phi3_dq_class_full - Phi3_dq_file_full;
        if any(abs([test_Phi3dq(:); test_Phi3dq_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi3dq(:); test_Phi3dq_full(:)]))
          error('Berechnung von constr3grad_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr3gradD_q
        if mextest == 1
          eval(sprintf('[Phi3_dqD_file, Phi3_dqD_file_full]=%s_constr3gradD_q(q, qD, x, xD, s_qD_3);', PName_Legs));
        else
          eval(sprintf('[Phi3_dqD_file, Phi3_dqD_file_full]=%s_constr3gradD_q_mex(q, qD, x, xD, s_qD_3);', PName_Legs));
        end
        [Phi3_dqD_class, Phi3_dqD_class_full] = R.constr3gradD_q(q, qD, x, xD);
        test_Phi3dqD = Phi3_dqD_class - Phi3_dqD_file;
        test_Phi3dqD_rel = test_Phi3dqD ./ Phi3_dqD_class;
        test_Phi3dqD_full = Phi3_dqD_class_full - Phi3_dqD_file_full;
        test_Phi3dqD_full_rel = test_Phi3dqD_full ./ Phi3_dqD_class_full;
        if any(abs([test_Phi3dqD(:);test_Phi3dqD_full(:)]) > 1e-10 & ...
            abs([test_Phi3dqD_rel(:);test_Phi3dqD_full_rel(:)]) > 1e-6) || ...
            any(isnan([test_Phi3dqD(:);test_Phi3dqD_full(:)]))
          error('Berechnung von constr3gradD_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr3grad_x
        if mextest == 1
          eval(sprintf('[Phi3_dx_file,Phi3_dx_file_full]=%s_constr3grad_x(q, x, s_x_3);', PName_Legs));
        else
          eval(sprintf('[Phi3_dx_file,Phi3_dx_file_full]=%s_constr3grad_x_mex(q, x, s_x_3);', PName_Legs));
        end
        [Phi3_dx_class,Phi3_dx_class_full] = R.constr3grad_x(q, x);
        test_Phi3dx = Phi3_dx_class - Phi3_dx_file;
        test_Phi3dx_full = Phi3_dx_class_full - Phi3_dx_file_full;
        if any(abs([test_Phi3dx(:);test_Phi3dx_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi3dx(:);test_Phi3dx_full(:)]))
          error('Berechnung von constr3grad_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr3gradD_x
        if mextest == 1
          eval(sprintf('[Phi3_dxD_file,Phi3_dxD_file_full]=%s_constr3gradD_x(q, qD, x, xD, s_xD_3);', PName_Legs));
        else
          eval(sprintf('[Phi3_dxD_file,Phi3_dxD_file_full]=%s_constr3gradD_x_mex(q, qD, x, xD, s_xD_3);', PName_Legs));
        end
        [Phi3_dxD_class,Phi3_dxD_class_full] = R.constr3gradD_x(q, qD, x, xD);
        test_Phi3dxD = Phi3_dxD_class - Phi3_dxD_file;
        test_Phi3dxD_rel = test_Phi3dxD ./ Phi3_dxD_class;
        test_Phi3dxD_full = Phi3_dxD_class_full - Phi3_dxD_file_full;
        test_Phi3dxD_full_rel = test_Phi3dxD_full ./ Phi3_dxD_class_full;
        if any(abs([test_Phi3dxD(:);test_Phi3dxD_full(:)]) > 1e-10 & ...
            abs([test_Phi3dxD_rel(:);test_Phi3dxD_full_rel(:)]) > 1e-6) || ...
            any(isnan([test_Phi3dxD(:);test_Phi3dxD_full(:)]))
          error('Berechnung von constr3gradD_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
        if i_FG == 2
          continue % Für 3T0R nicht definiert.
        end
        %% Teste constr2grad_q
        if mextest == 1
          eval(sprintf('[Phi2_dq_file,Phi2_dq_file_full]=%s_constr2grad_q(q, x, s_q_2);', PName_Legs));
        else
          eval(sprintf('[Phi2_dq_file,Phi2_dq_file_full]=%s_constr2grad_q_mex(q, x, s_q_2);', PName_Legs));
        end
        [Phi2_dq_class, Phi2_dq_class_full] = R.constr2grad_q(q, x);
        test_Phi2dq = Phi2_dq_class - Phi2_dq_file;
        test_Phi2dq_full = Phi2_dq_class_full - Phi2_dq_file_full;
        if any(abs([test_Phi2dq(:); test_Phi2dq_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi2dq(:); test_Phi2dq_full(:)]))
          error('Berechnung von constr2grad_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr2gradD_q
        if mextest == 1
          eval(sprintf('[Phi2D_dq_file,Phi2D_dq_file_full]=%s_constr2gradD_q(q, qD, x, xD, s_qD_2);', PName_Legs));
        else
          eval(sprintf('[Phi2D_dq_file,Phi2D_dq_file_full]=%s_constr2gradD_q_mex(q, qD, x, xD, s_qD_2);', PName_Legs));
        end
        [Phi2D_dq_class,Phi2D_dq_class_full] = R.constr2gradD_q(q, qD, x, xD);
        % Absoluten und relativen Fehler gleichzeitig testen, da 1e-10 auch
        % mal knapp überschritten werden kann.
        test_Phi2Ddq = Phi2D_dq_class - Phi2D_dq_file;
        test_Phi2Ddq_rel = test_Phi2Ddq ./ Phi2D_dq_class;
        test_Phi2Ddq_full = Phi2D_dq_class_full - Phi2D_dq_file_full;
        test_Phi2Ddq_full_rel = test_Phi2Ddq_full ./ Phi2D_dq_class_full;
        if any(abs([test_Phi2Ddq(:);test_Phi2Ddq_full(:)])         > 1e-10 & ...
               abs([test_Phi2Ddq_rel(:);test_Phi2Ddq_full_rel(:)]) > 1e-6) || ...
            any(isnan([test_Phi2Ddq(:);test_Phi2Ddq_full(:)]))
          error('Berechnung von constr2gradD_q stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr2grad_rr
        if mextest == 1
          eval(sprintf('[Phi2_rr_file,Phi2_rr_full_file]=%s_constr2grad_rr(q, x, s_x_4);', PName_Legs));
        else
          eval(sprintf('[Phi2_rr_file,Phi2_rr_full_file]=%s_constr2grad_rr_mex(q, x, s_x_4);', PName_Legs));
        end
        [Phi2_rr_class, Phi2_rr_full_class] = R.constr2grad_rr(q, x);
        test_Phi2rr = Phi2_rr_class - Phi2_rr_file;
        test_Phi2rr_full = Phi2_rr_full_class - Phi2_rr_full_file;
        if any(abs([test_Phi2rr(:);test_Phi2rr_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi2rr(:);test_Phi2rr_full(:)]))
          error('Berechnung von constr2grad_rr stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr2grad_x
        if mextest == 1
          eval(sprintf('[Phi2_dx_file,Phi2_dx_full_file]=%s_constr2grad_x(q, x, s_x_2);', PName_Legs));
        else
          eval(sprintf('[Phi2_dx_file,Phi2_dx_full_file]=%s_constr2grad_x_mex(q, x, s_x_2);', PName_Legs));
        end
        [Phi2_dx_class, Phi2_dx_full_class] = R.constr2grad_x(q, x);
        test_Phi2dx = Phi2_dx_class - Phi2_dx_file;
        test_Phi2dx_full = Phi2_dx_full_class - Phi2_dx_full_file;
        if any(abs([test_Phi2dx(:);test_Phi2dx_full(:)]) > 1e-10) || ...
            any(isnan([test_Phi2dx(:);test_Phi2dx_full(:)]))
          error('Berechnung von constr2grad_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
        %% Teste constr2gradD_rr
        if mextest == 1
          eval(sprintf('[Phi2D_rr_file,Phi2D_rr_file_full]=%s_constr2gradD_rr(q, qD, x, xD, s_xD_4);', PName_Legs));
        else
          eval(sprintf('[Phi2D_rr_file,Phi2D_rr_file_full]=%s_constr2gradD_rr_mex(q, qD, x, xD, s_xD_4);', PName_Legs));
        end
        [Phi2D_rr_class, Phi2D_rr_class_full] = R.constr2gradD_rr(q, qD, x, xD);
        test_Phi2Drr = Phi2D_rr_class - Phi2D_rr_file;
        test_Phi2Drr_rel = test_Phi2Drr ./ Phi2D_rr_class;
        test_Phi2Drr_full = Phi2D_rr_class_full - Phi2D_rr_file_full;
        test_Phi2Drr_full_rel = test_Phi2Drr_full ./ Phi2D_rr_class_full;
        if any(abs([test_Phi2Drr(:);test_Phi2Drr_full(:)])         > 1e-10 & ...
               abs([test_Phi2Drr_rel(:);test_Phi2Drr_full_rel(:)]) > 1e-6) || ...
            any(isnan([test_Phi2Drr(:);test_Phi2Drr_full(:)]))
          error(['Berechnung von constr2gradD_rr stimmt nicht zwischen ', ...
            'Template-Funktion und Klasse. Fehler %1.1e bzw. %1.1e'], ...
            max(abs(test_Phi2Drr(:))), max(abs(test_Phi2Drr_full(:))));
        end
        %% Teste constr2gradD_x
        if mextest == 1
          eval(sprintf('[Phi2D_dx_file,Phi2D_dx_file_full]=%s_constr2gradD_x(q, qD, x, xD, s_xD_2);', PName_Legs));
        else
          eval(sprintf('[Phi2D_dx_file,Phi2D_dx_file_full]=%s_constr2gradD_x_mex(q, qD, x, xD, s_xD_2);', PName_Legs));
        end
        [Phi2D_dx_class, Phi2D_dx_class_full] = R.constr2gradD_x(q, qD, x, xD);
        test_Phi2Ddx = Phi2D_dx_class - Phi2D_dx_file;
        test_Phi2Ddx_rel = test_Phi2Ddx ./ Phi2D_dx_class;
        test_Phi2Ddx_full = Phi2D_dx_class_full - Phi2D_dx_file_full;
        test_Phi2Ddx_full_rel = test_Phi2Ddx_full ./ Phi2D_dx_class_full;
        if any(abs([test_Phi2Ddx(:);test_Phi2Ddx_full(:)]) > 1e-10 & ...
               abs([test_Phi2Ddx_rel(:);test_Phi2Ddx_full_rel(:)]) > 1e-10) || ...
            any(isnan([test_Phi2Ddx(:);test_Phi2Ddx_full(:)]))
          error('Berechnung von constr2gradD_x stimmt nicht zwischen Template-Funktion und Klasse');
        end
      end
      fprintf('PKM %s getestet: mex=%d, taskred=%d\n', PName, mextest-1, trtest-1);
    end % for mextest
    end % for trtest
    fprintf('PKM %s erfolgreich getestet\n', PName);
  end % for ii (PKM)
end % for EE_FG
