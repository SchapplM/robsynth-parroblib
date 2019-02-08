% Instanz der Roboterklasse für gegebenen Roboter initialisieren
% 
% Eingabe:
% Name
%   Name des PKM-Robotermodells nach dem Schema "PxRRRyyAzz"
% d_0A
%   Durchmesser des Kreises, auf dem die Basis-Koppelpunkte liegen
% d_PB
%   Durchmesser des Kreises der Plattform-Koppelpunkte
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

function RP = parroblib_create_robot_class(Name, d_0A, d_PB, phi_RS_EE)

if nargin < 4
  phi_RS_EE = [];
end

%% Daten laden

[NLEG, LEG_Names, Actuation, ActNr, ~, EE_dof0, PName_Kin] = parroblib_load_robot(Name);

parroblibpath=fileparts(which('parroblib_path_init.m'));

%% Instanz der seriellen Roboterklasse erstellen
% TODO: Nicht-symmetrische PKM fehlen noch
RS = serroblib_create_robot_class(LEG_Names{1});
RS.fill_fcn_handles(false);
RS.update_pkin();

RS.I_EE = logical(EE_dof0); % Für IK der Beinketten mit invkin_ser

% TODO: Das ist keine automatische Lösung
% EE anpassen für 2T1R-PKM, bei denen die Plattform-KoppelKS falsch gedreht
% sind.
if ~isempty(phi_RS_EE)
  RS.update_EE([], phi_RS_EE);
end
%% Instanz der parallelen Roboterklasse erstellen

% Pfade für Matlab-Funktionen hinzufügen
addpath(fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Kin, 'hd_A0'));
addpath(fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Kin, sprintf('hd_A%d', ActNr)));

RP = ParRob(Name);
RP = RP.create_symmetric_robot(NLEG, RS, d_0A, d_PB);
RP = RP.initialize();

% EE-FG eintragen
RP.update_EE_FG(logical(EE_dof0)); % Für IK der PKM

% Aktuierung eintragen
I_qa = false(RP.NJ,1);
for i = 1:NLEG
  I_Legi = RP.I1J_LEG(i):RP.I2J_LEG(i);
  I_qa(I_Legi(Actuation{i})) = true;
end
RP.update_actuation(I_qa);

RP.fill_fcn_handles(false);
