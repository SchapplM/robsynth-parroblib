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
  [NLEG, LEG_Names, ~, Coupling, ActNr, ~, ~, PName_Kin, PName_Legs] = parroblib_load_robot(Name);

  % Pfade für Matlab-Funktionen der PKM hinzufügen
  fcn_dir1 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Legs, ...
    sprintf('hd_G%dP%dA0', Coupling(1), Coupling(2)));
  fcn_dir2 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Legs, ...
    sprintf('hd_G%dP%dA%d', Coupling(1), Coupling(2), ActNr));
  fcn_dir3 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Legs, ...
    'tpl');
  if exist(fcn_dir1, 'file') % existiert nur, wenn symbolischer Code generiert
    addpath(fcn_dir1);
  end
  if exist(fcn_dir2, 'file') % existiert nur, wenn für diese Aktuierung symbolisch generiert
    addpath(fcn_dir2);
  end
  if ~exist(fcn_dir3, 'file')
    fprintf('Vorlagen-Funktion existieren nicht für %s. Erstelle.\n', PName_Kin)
    parroblib_create_template_functions({PName_Kin});
  end
  addpath(fcn_dir3);
  % Zusätzlich Pfade für Funktionen der seriellen Beinketten hinzufügen.
  % Damit wird das Laden eins
  LEG_Names_unique = unique(LEG_Names);
  for j = 1:length(LEG_Names_unique)
    serroblib_addtopath({LEG_Names_unique{j}});
  end
end