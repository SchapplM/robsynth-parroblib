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
else
  error('Plattform-Methode %d nicht definiert', Coupling(1));
end
RP.align_platform_coupling(Coupling(2), p_platform_all);

% EE-FG eintragen
if all(EE_dof0 == [1 1 0 0 0 1]) % 2T1R (planar)
  RP.update_EE_FG(logical(EE_dof0), logical(EE_dof0), logical(repmat(logical(EE_dof0),RP.NLEG,1)));
elseif all(EE_dof0 == [1 1 1 1 1 0])
  error('3T2R-Aufgaben für PKM noch nicht implementiert');
else
  RP.update_EE_FG(logical(EE_dof0), logical(EE_dof0), true(RP.NLEG,6));
end

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
