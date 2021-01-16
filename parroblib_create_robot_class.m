% Instanz der Roboterklasse für gegebenen Roboter initialisieren
% 
% Eingabe:
% Name
%   Name des PKM-Robotermodells nach dem Schema "PxRRRyyGuuPvvAzz"
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

function RP = parroblib_create_robot_class(Name, p_Base, p_platform, phi_RS_EE)

if nargin < 4
  phi_RS_EE = [];
end
assert(isa(Name,'char'), 'Eingabe Name muss Name als char sein');

%% Daten laden

[NLEG,LEG_Names,Actuation,Coupling,~,~, EE_dof0]=parroblib_load_robot(Name);

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
I_robot = find(strcmp(l.Names_Ndof,LEG_Names{1}));
joints_number_str = sprintf('%d', l.AdditionalInfo(I_robot,7));
numtechjoints = l.AdditionalInfo(I_robot,6);
i_coord = 0; % Zähler über Gelenke im Modell
for i = 1:numtechjoints % über technische Gelenke (können mehrwertig sein)
  % Der Gelenktyp ist rückwärts gelesen kodiert (einer-Stelle ist erstes
  % Gelenk).
  joint_i = str2double(joints_number_str(end-i+1));
  switch joint_i % 1=R,2=P,3=C,4=U,5=S
    case 1 % R% 0=Drehgelenk, 1=Schub- (allgemein), 2=Kardan-, 3=Kugel
      i_coord = i_coord + 1;
      RS.DesPar.joint_type(i_coord) = 0; % Drehgelenk
    case 2 % P
      i_coord = i_coord + 1;
      RS.DesPar.joint_type(i_coord) = 1; % Schubgelenk (allgemein)
    case 3 % C
      i_coord = i_coord + 2;
      % Dreh-Schubgelenk noch nicht implementiert. Setze einfach direkt die
      % MDH-Parameter. Keine Auswirkung.
      RS.DesPar.joint_type(i_coord-1:i_coord) = RS.MDH.sigma(i_coord-1:i_coord);
      warning('Gelenk Typ C noch nicht implementiert');
    case 4 % U
      i_coord = i_coord + 2;
      RS.DesPar.joint_type(i_coord-1:i_coord) = 2; % Kardan
    case 5 % S
      i_coord = i_coord + 3;
      RS.DesPar.joint_type(i_coord-2:i_coord) = 3; % Kugel
  end
end
%% Instanz der parallelen Roboterklasse erstellen

parroblib_addtopath({Name})

RP = ParRob(Name);
RP.create_symmetric_robot(NLEG, RS);
RP.initialize();
% Vervollständige Koppelpunkt-Parameter mit Standard-Einstellungen
if length(p_Base) > 1
  % Annahme: Bei Vorgabe mehrere Parameter hat der Benutzer alle
  % notwendigen Parameter angegeben und weiß was er tut.
  p_Base_all = p_Base;
elseif any(Coupling(1) == [1,2,3])
  p_Base_all = p_Base;
elseif Coupling(1) == 4
  % Pyramide symmetrische Anordnung. Nehme standardmäßig halben
  % Punktradius als Punktabstand und 30 Grad Steigung
  p_Base_all = [p_Base; 30*pi/180];
elseif any(Coupling(1) == [5,6,7])
  % Paarweiser Anordnung. Nehme standardmäßig halben
  p_Base_all = [p_Base; p_Base/2];
elseif Coupling(1) == 8
  % Pyramide mit paarweiser Anordnung. Nehme standardmäßig halben
  % Punktradius als Punktabstand und 30 Grad Steigung
  p_Base_all = [p_Base; p_Base/2; 30*pi/180];
elseif Coupling(1) == 9
  p_Base_all = p_Base;
else
  error('Gestell-Methode %d nicht definiert', Coupling(1));
end
RP.align_base_coupling(Coupling(1), p_Base_all);
if length(p_platform) > 1
  % Annahme: Bei Vorgabe mehrere Parameter hat der Benutzer alle
  % notwendigen Parameter angegeben und weiß was er tut.
  p_platform_all = p_platform;
elseif any(Coupling(2) == [1,2,3 8])
  p_platform_all = p_platform;
elseif any(Coupling(2) == [4,5,6])
  p_platform_all = [p_platform; p_platform/2];
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
