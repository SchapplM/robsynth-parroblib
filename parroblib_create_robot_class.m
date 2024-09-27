% Instanz der Roboterklasse für gegebenen Roboter initialisieren
% 
% Eingabe:
% Name
%   Name des PKM-Robotermodells nach dem Schema "PxRRRyyGuuPvvAzz"
% RobName
%   Name der Roboterparameter entsprechend der Tabelle models.csv
% p_Base
%   Parameter der Gestell-Koppelpunkte. Z.B. Radius des Kreises.
%   Siehe ParRob/align_base_coupling
% p_platform
%   Parameter der Plattform-Koppelpunkte. Z.B. Radius des Kreises der
%   Plattform-Koppelpunkte
%   Siehe ParRob/align_platform_coupling
% phi_RS_EE (optional)
%   Orientierung des Beinketten-Koppel-KS. Bei einigen planaren Systemen
%   muss eine Drehung um Pi erfolgen, damit das Bein-Koppel-KS mit dem
%   Plattform-KS übereinstimmen kann.
%   Diese Drehung ist eigentlich schon in der SerRobLib gespeichert
% 
% Ausgabe:
% RP [ParRob]
%   Instanz der ParRob-Klasse mit den Eigenschaften und Funktionen des zu
%   ladenden Roboters
% 
% TODO:
% * EE-FG für 5FG-Roboter anders definieren
% * Nicht-symmetrische PKM berücksichtigen
% 
% Siehe auch serroblib_create_robot_class

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function RP = parroblib_create_robot_class(Name, RobName, p_Base, p_platform, phi_RS_EE)

if nargin < 5
  phi_RS_EE = [];
end
assert(isa(Name,'char'), 'Eingabe Name muss Name als char sein');
assert(isa(RobName,'char') || isempty(RobName), ...
  'Eingabe Roboter-Modellname muss Name als char sein');

%% Daten laden

[NLEG,LEG_Names,Actuation,Coupling,~,~,EE_dof0,~,PName_Legs]=parroblib_load_robot(Name);

%% Instanz der seriellen Roboterklasse erstellen
% TODO: Nicht-symmetrische PKM fehlen noch
RS = serroblib_create_robot_class(LEG_Names{1});
RS.fill_fcn_handles(false);
RS.update_pkin();
% Entferne die gespeicherte EE-Transformation der seriellen Kette.
% Diese Transformation ist nur für serielle Roboter relevant. Für PKM wird
% die Transformation durch die virtuellen KS "P->Bi" erledigt.
% Nur eine Drehung um 180° ist teilweise noch notwendig. Das betrifft 
% Beinketten mit nur einem rotatorischen FG.
if sum(RS.I_EE(4:6)) > 1
  RS.update_EE(zeros(3,1),zeros(3,1));
end

% Manuelle Eingabe verarbeiten
if ~isempty(phi_RS_EE)
  RS.update_EE([], phi_RS_EE);
end

%% Technische Gelenke in Beinkette eintragen
% Das wird nicht direkt beim Erstellen der seriellen Beinkette gemacht, da
% dort noch nicht bekannt ist, ob die serielle Kette ein serieller Roboter
% oder eine PKM-Beinkette ist. Siehe serroblib_gen_bitarrays.m.
% Laden der Seriellroboter-Datenbank und extrahieren der Information:
N = str2double(LEG_Names{1}(2));
mdllistfile_Ndof = fullfile(fileparts(which('serroblib_path_init.m')), ...
  sprintf('mdl_%ddof', N), sprintf('S%d_list.mat',N));
l = load(mdllistfile_Ndof, 'Names_Ndof', 'AdditionalInfo');
I_robot = strcmp(l.Names_Ndof,LEG_Names{1});
SName_TechJoint = fliplr(regexprep(num2str(l.AdditionalInfo(I_robot,7)), ...
    {'1','2','3','4','5'}, {'R','P','C','U','S'}));
RS.set_techjoints(SName_TechJoint);

%% Parameter-Zahlenwerte setzen (falls Robotermodell geladen wird)
if ~isempty(RobName)
  pardat = fullfile(fileparts(which('parroblib_path_init.m')), ...
    sprintf('sym_%dT%dR', sum(EE_dof0(1:3)), sum(EE_dof0(4:6))), PName_Legs, 'models.csv');
  assert(exist(pardat, 'file'), sprintf('Parameter-Datei existiert nicht: %s', pardat));
  ParTable = readtable(pardat, 'Delimiter', ';');
  unitmult_angle = pi/180; % Angaben in Tabelle in Grad. Umrechnung in Radiant
  unitmult_dist = 1/1000; % Angaben in Millimeter
  % Überschreibe die PKM-Parameter aus der Eingabe
  I = strcmp(ParTable.Name, RobName);
  assert(sum(I)==1, sprintf(['Name %s wurde nicht eindeutig in Tabelle ' ...
    'gefunden: %s'], RobName, pardat));
  p_Base_tmp = [unitmult_dist*ParTable.base_radius(I); ...
                unitmult_dist*ParTable.base_morph_pairdist(I)];
  p_Base = p_Base_tmp(~isnan(p_Base_tmp));
  p_platform_tmp = [unitmult_dist*ParTable.platform_radius(I); ...
                    unitmult_dist*ParTable.platform_morph_pairdist(I)];
  p_platform = p_platform_tmp(~isnan(p_platform_tmp));
  % Setze die DH-Parameter für die Beinkette
  pkin_neu = RS.pkin;
  for i = 1:length(RS.pkin_names)
    VarColName = sprintf('pkin_%s', RS.pkin_names{i});
    if ~any(strcmp(ParTable.Properties.VariableNames, VarColName))
      TabCol_i = 0; % Parameter ist in Tabelle nicht definiert. Daher wohl Null.
    else
      if any(contains(VarColName, {'_alpha', '_theta'}))
        unitmult_i = unitmult_angle;
      else
        unitmult_i = unitmult_dist;
      end
      TabCol_i = unitmult_i * ParTable.(sprintf('pkin_%s', RS.pkin_names{i}));
    end
    pkin_neu(i) = TabCol_i(I);
  end
  % Eintragen der neuen Kinematikparameter
  RS.update_mdh(pkin_neu);
  % Setze die Gelenkwinkelgrenzen
  qlim = RS.qlim;
  for i = 1:RS.NQJ
    if any(RS.MDH.sigma(i) == 0)
      unitmult_i = unitmult_angle;
    else
      unitmult_i = unitmult_dist;
    end
    qmin_tab = ParTable.(sprintf('q%dmin',i));
    qmax_tab = ParTable.(sprintf('q%dmax',i));
    qlim(i,:) = unitmult_i * [qmin_tab(I), qmax_tab(I)];
  end
  RS.update_qlim(qlim);
end
%% Instanz der parallelen Roboterklasse erstellen
parroblib_addtopath({Name})

RP = ParRob(Name);
RP.create_symmetric_robot(NLEG, RS);
RP.I_EE = logical(EE_dof0); % Schon hier belegnDient zu Logik-Prüfungen in Methode initialize
RP.initialize();
% Vervollständige Koppelpunkt-Parameter mit Standard-Einstellungen
if length(p_Base) > 1
  % Annahme: Bei Vorgabe mehrere Parameter hat der Benutzer alle
  % notwendigen Parameter angegeben und weiß was er tut.
  p_Base_all = p_Base;
elseif any(Coupling(1) == [1,2,3])
  % Einziger Parameter ist der Gestell-Radius
  p_Base_all = p_Base;
elseif Coupling(1) == 4
  % Pyramide symmetrische Anordnung. Nehme standardmäßig halben
  % Punktradius als Punktabstand und 30 Grad Steigung
  p_Base_all = [p_Base; 30*pi/180];
elseif any(Coupling(1) == [5,6,7])
  % Paarweiser Anordnung. Nehme standardmäßig halben Radius als
  % Paar-Abstand
  p_Base_all = [p_Base; p_Base/2];
elseif Coupling(1) == 8
  % Pyramide mit paarweiser Anordnung. Nehme standardmäßig halben
  % Punktradius als Punktabstand und 30 Grad Steigung
  p_Base_all = [p_Base; p_Base/2; 30*pi/180];
elseif any(Coupling(1) == [9 10])
  % Einziger Parameter ist der Gestell-Radius
  p_Base_all = p_Base;
else
  error('Gestell-Methode %d nicht definiert', Coupling(1));
end
RP.align_base_coupling(Coupling(1), p_Base_all);
if length(p_platform) > 1
  % Annahme: Bei Vorgabe mehrere Parameter hat der Benutzer alle
  % notwendigen Parameter angegeben und weiß was er tut.
  p_platform_all = p_platform;
elseif any(Coupling(2) == [1,2,3])
  p_platform_all = p_platform;
elseif any(Coupling(2) == [4,5,6])
  p_platform_all = [p_platform; p_platform/2];
elseif Coupling(2) == 7
  p_platform_all = p_platform;
elseif Coupling(2) == 8
  p_platform_all = [p_platform; 0];
elseif any(Coupling(2) == [9 10]) % konisch (=vertikal/geneigt) und tangential/geneigt
  p_platform_all = [p_platform; 0];
else
  error('Plattform-Methode %d nicht definiert', Coupling(2));
end
RP.align_platform_coupling(Coupling(2), p_platform_all);

% EE-FG eintragen. Die Logik für die Zuordnung der Beinketten-FG liegt in
% der Klassenmethode. Sonderfall von PKM mit Beinketten, deren FG identisch
% zur Plattform sind, noch nicht implementiert (für 3T0R, 3T1R)
RP.update_EE_FG(logical(EE_dof0), logical(EE_dof0));

% Aktuierung eintragen
I_qa = false(RP.NJ,1);
for i = 1:NLEG
  if isempty(Actuation{i})
    continue
  end
  I_Legi = RP.I1J_LEG(i):RP.I2J_LEG(i);
  I_qa(I_Legi(Actuation{i})) = true;
end
RP.update_actuation(I_qa);

RP.fill_fcn_handles(false);
