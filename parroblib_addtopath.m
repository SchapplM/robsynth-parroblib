% Füge die Modelle gegebener Roboter zum Pfad hinzu
% 
% Eingabe:
% Names
%   Cell-Array mit Namen der Robotermodelle, deren Funktionen zum
%   Matlab-Pfad hinzugefügt werden sollen.
% 
% Siehe auch: serroblib_addtopath.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_addtopath(Names)
if ~iscell(Names)
  error('parroblib_addtopath: Eingegebene Roboternamen müssen ein Cell-Array sein');
end
parroblibpath=fileparts(which('parroblib_path_init.m'));
for i = 1:length(Names)
  Name = Names{i};
  [NLEG, LEG_Names, ~, Coupling, ActNr, ~, ~, ~, PName_Legs] = parroblib_load_robot(Name);

  % Pfade für Matlab-Funktionen der PKM hinzufügen
  fcn_dir1 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Legs, ...
    sprintf('hd_G%dP%dA0', Coupling(1), Coupling(2)));
  fcn_dir2 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Legs, ...
    sprintf('hd_G%dP%dA%d', Coupling(1), Coupling(2), ActNr));
  if exist(fcn_dir1, 'file')
    addpath(fcn_dir1);
  end
  if exist(fcn_dir2, 'file')
    addpath(fcn_dir2);
  end
  % Zusätzlich Pfade für Funktionen der seriellen Beinketten hinzufügen.
  % Damit wird das Laden eins
  LEG_Names_unique = unique(LEG_Names);
  for j = 1:length(LEG_Names_unique)
    serroblib_addtopath({LEG_Names_unique{j}});
  end
end