% Aufruf und Test aller symbolisch generierten Funktionen der Roboter
% Dadurch wird sichergestellt, dass keine Fehler in den Funktionen vorliegen.
% 
% Funktionsweise:
% * In Ordner wechseln, in dem die Ergebnisse gespeichert werden sollen
% * Durchlauf mit usr_testfaillist=0 zur Erkennung mit usr_codegen=0
% * Durchlauf zur Neu-Generierung mit usr_testfaillist=1 und usr_codegen=1

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-11
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc
repopath=fileparts(which('serroblib_path_init.m'));
usr_testselection = false;
usr_shuffle_pkm_selection = true;
usr_testfaillist = true; % Zuerst Durchlauf notwendig um Listendatei zu erstellen
usr_codegen = true;
%% Alle PKM durchgehen
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 1];
robot_list_succ = {}; % Liste erfolgreicher PKM
robot_list_notest = {}; % Liste ohne Test-Möglichkeit
robot_list_fail = {}; % Liste fehlgeschlagener
if usr_testfaillist
  d = load('parroblib_call_all_fcns_result.mat');
  FiltList = d.robot_list_fail;
else
  FiltList = {};
end
if usr_testselection
  FiltList = {'P3RRR1G1P1', 'P3PRRRR8V2G1P2', 'P4PRRRR8V1G3P1', 'P6RRPRRR14V4G7P4', 'P6RRRRRR10V6G6P1A1'};
end
for i_FG = 1:size(EEFG_Ges,1)
  EE_FG = EEFG_Ges(i_FG,:);
  [PNames_Kin, PNames_Akt] = parroblib_filter_robots(EE_FG, 6);
  III = 1:length(PNames_Kin);
  if usr_shuffle_pkm_selection
    III = III(randperm(length(III)));
  end
  if isempty(III), warning('Keine PKM mit FG %s gefunden', EE_FGstr); continue; end
  for ii = III(:)' % 1:length(PNames_Kin)
    % Suche erste gültige Aktuierung in der Datenbank (ist für Dynamik
    % egal; es wird nur die Plattform-Dynamik betrachtet ohne Projektion in
    % die Antriebskoordinaten).
    IIa = find(contains(PNames_Akt, [PNames_Kin{ii},'A']),1,'first');
    PName = PNames_Akt{IIa};
    if ~isempty(FiltList) && ~any(contains(PName, FiltList))
      continue;
    end
    fprintf('Untersuche PKM %s\n', PName);
    %% Roboter Initialisieren
    R = parroblib_create_robot_class(PName, '', 1, 0.3);
    % Initialisierung der Funktionen: Kompilierte Funktionen nehmen
    files_missing = R.fill_fcn_handles(false, false);
    % Prüfe, ob die Dynamik-Implementierung in symbolischer Form vorliegt
    if length(files_missing) >= 6
      fprintf('Keine symbolisch generierten Dynamik-Funktionen verfügbar.\n');
      robot_list_notest = [robot_list_notest(:)', {PName}];
      continue
    end
    if any(contains(files_missing, 'Jinv'))
      fprintf('Inverse Jacobi-Matrix ist nicht in symbolischer Form vorhanden\n');
      % Ist nicht so schlimm. Eigentlich geht es ja um die Dynamik.
    end
    % Initialisierung der Funktionen: Keine kompilierte Funktionen nehmen
    R.fill_fcn_handles(false, false);
    %% Kinematikparameter zufällig setzen (ist egal)
    pkin_rand = rand(length(R.Leg(1).pkin),1); % Symmetrische PKM mit zufälligen Parametern
    for il = 1:R.NLEG
      R.Leg(il).update_mdh(pkin_rand);
      R.Leg(il).qlim = repmat([-pi,pi], R.Leg(il).NJ,1);
    end
    R.update_base(rand(3,1), rand(3,1));
    R.update_EE(rand(3,1), rand(3,1));
    %% Parameter für Dynamik-Test
    mges_PKM = rand(size(R.DynPar.mges));
    rSges_PKM = rand(size(R.DynPar.rSges));
    ISges_PKM = rand(size(R.DynPar.Icges));
    mges_PKM(R.NQJ_LEG_bc+1:end-1) = 0;
    rSges_PKM(R.NQJ_LEG_bc+1:end-1,:) = 0;
    ISges_PKM(R.NQJ_LEG_bc+1:end-1,:) = 0;
    R.update_dynpar1 (mges_PKM, rSges_PKM, ISges_PKM);
    %% Test-Szenario definieren
    qlim_pkm = cat(1,R.Leg(:).qlim);
    n = 5; % Verschiedene Konfigurationen
    Q_test = qlim_pkm(:,1)' + (qlim_pkm(:,2)-qlim_pkm(:,1))'.*rand(n, R.NJ);
    QD_test = rand(n, R.NJ);
    QDD_test = rand(n, R.NJ);
    X_test = rand(n, 6);
    XD_test = rand(n, 6);
    XDD_test = rand(n, 6);
    X_test(:,~R.I_EE) = 0;
    XD_test(:,~R.I_EE) = 0;
    XDD_test(:,~R.I_EE) = 0;
    %% Funktionen aufrufen
    valid = true;
    for i = 1:n
      R.DynPar.mode = 1;
      Mx1 = R.inertia_platform(Q_test(i,:)', X_test(i,:)');
      Cx1 = R.coriolisvec_platform(Q_test(i,:)', X_test(i,:)', XD_test(i,:)');
      Gx1 = R.gravload_platform(Q_test(i,:)', X_test(i,:)');
      Fx1 = R.invdyn_platform(Q_test(i,:)', X_test(i,:)', XD_test(i,:)', XDD_test(i,:)');
      R.DynPar.mode = 2;
      Mx2 = R.inertia_platform(Q_test(i,:)', X_test(i,:)');
      Cx2 = R.coriolisvec_platform(Q_test(i,:)', X_test(i,:)', XD_test(i,:)');
      Gx2 = R.gravload_platform(Q_test(i,:)', X_test(i,:)');
      Fx2 = R.invdyn_platform(Q_test(i,:)', X_test(i,:)', XD_test(i,:)', XDD_test(i,:)');
      % Prüfe Dynamikparameter 1 gegen Nummer 2 (wird eigentlich schon bei
      % Generierung gemacht, aber hier nochmal, falls dort nicht gemacht).
      if any(abs([Cx1;Cx2]) > 1e8), continue; end  % Singularität
      test_G_abs = Gx1-Gx2;
      test_C_abs = Cx1-Cx2;
      test_F_abs = Fx1-Fx2;
      test_M_abs = Mx1-Mx2;
      test_C_rel = test_C_abs ./ Cx2;
      test_G_rel = test_G_abs ./ Gx2;
      test_F_rel = test_F_abs ./ Fx2;
      % Prüfe Dynamikparameter mit Regressorform
      if ~any(contains(files_missing, 'para_pf_mdp'))
        R.DynPar.mode = 4;
        [Mx4, Mx4reg] = R.inertia_platform(Q_test(i,:)', X_test(i,:)');
        [Cx4, Cx4reg] = R.coriolisvec_platform(Q_test(i,:)', X_test(i,:)', XD_test(i,:)');
        [Gx4, Gx4reg] = R.gravload_platform(Q_test(i,:)', X_test(i,:)');
        [Fx4, Fx4reg] = R.invdyn_platform(Q_test(i,:)', X_test(i,:)', XD_test(i,:)', XDD_test(i,:)');
        test_G4_abs = Gx1-Gx4;
        test_C4_abs = Cx1-Cx4;
        test_F4_abs = Fx1-Fx4;
        test_M4_abs = Mx1-Mx4;
        test_C4_rel = test_C4_abs ./ Cx1;
        test_G4_rel = test_G4_abs ./ Gx1;
        test_F4_rel = test_F4_abs ./ Fx1;
        test_G4_reg = Gx4reg*R.DynPar.mpv_sym-Gx4;
        test_C4_reg = Cx4reg*R.DynPar.mpv_sym-Cx4;
        test_F4_reg = Fx4reg*R.DynPar.mpv_sym-Fx4;
        test_M4_reg = reshape(Mx4reg*R.DynPar.mpv_sym, R.NLEG, R.NLEG)' - Mx4;
      end
      try
        assert(all(abs(test_G_abs)<1e-8) | all(abs(test_G_rel)<1e-8), ...
          'Gravitationsmoment stimmt nicht zwischen DynPar Methode 1 und 2');
        assert(all(abs(test_C_abs)<1e-8) | all(abs(test_C_rel)<1e-8), ...
          'Coriolismoment stimmt nicht zwischen DynPar Methode 1 und 2');
        assert(all(abs(test_F_abs)<1e-8) | all(abs(test_F_rel)<1e-8), ...
          'Inversynamik-Kraft stimmt nicht zwischen DynPar Methode 1 und 2');
        assert(all(abs(test_M_abs(:))<1e-8), 'Massenmatrix stimmt nicht zwischen DynPar Methode 1 und 2');
%         assert(~all(test_G==0), 'Alle Einträge in G-Termen exakt identisch zwischen DynPar Methode 1 und 2. Unlogisch, da Rundungsfehler vorliegen müssen.')
%         assert(~all(test_C_abs==0), 'Alle Einträge in C-Termen exakt identisch zwischen DynPar Methode 1 und 2. Unlogisch, da Rundungsfehler vorliegen müssen.')
%         assert(~all(test_F_abs==0), 'Alle Einträge in F-Termen exakt identisch zwischen DynPar Methode 1 und 2. Unlogisch, da Rundungsfehler vorliegen müssen.')
%         assert(~all(test_M_abs(:)==0), 'Alle Einträge in M-Termen exakt identisch zwischen DynPar Methode 1 und 2. Unlogisch, da Rundungsfehler vorliegen müssen.')
        if ~any(contains(files_missing, 'para_pf_mdp'))
          assert(all(abs(test_G4_abs)<1e-8) | all(abs(test_G4_rel)<1e-8), ...
            'Gravitationsmoment stimmt nicht zwischen DynPar Methode 1 und 4');
          assert(all(abs(test_C4_abs)<1e-8) | all(abs(test_C4_rel)<1e-8), ...
            'Coriolismoment stimmt nicht zwischen DynPar Methode 1 und 4');
          assert(all(abs(test_F4_abs)<1e-8) | all(abs(test_F4_rel)<1e-8), ...
            'Inversynamik-Kraft stimmt nicht zwischen DynPar Methode 1 und 4');
          assert(all(abs(test_M4_abs(:))<1e-8), 'Massenmatrix stimmt nicht zwischen DynPar Methode 1 und 4');
          assert(all(abs(test_F4_reg)<1e-8), 'Inversynamik-Regressor stimmt nicht');
          assert(all(abs(test_C4_reg)<1e-8), 'Coriolis-Regressor stimmt nicht');
          assert(all(abs(test_G4_reg)<1e-8), 'Gravitationslast-Regressor stimmt nicht');
          assert(all(abs(test_M4_reg(:))<1e-8), 'Massenmatrix-Regressor stimmt nicht');
        end
      catch err
        warning('parrob:fcnfail', 'Fehler in generiertem Code für %s: %s. Neu-Generierung notwendig.', PName, err.message);
        valid = false;
        break;
      end
    end
    fprintf('Dynamik-Terme für %s getestet\n', PName);
    if valid
      robot_list_succ = [robot_list_succ(:)', {PName}];
    else
      robot_list_fail = [robot_list_fail(:)', {PName}];
    end
    if valid == false && usr_codegen
      fprintf('Dynamik für %s neu generieren\n', PName)
      parroblib_generate_mapleinput({PName});
      parroblib_generate_code({PName}, true);
    end
  end
end
if ~usr_testfaillist
  save('parroblib_call_all_fcns_result.mat', 'robot_list_succ', 'robot_list_fail', 'robot_list_notest');
end
