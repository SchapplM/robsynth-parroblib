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
  switch Coupling(1) % siehe align_base_coupling.m
    case 1
      basecoupling_equiv = [5,4,8];
    case 2
      basecoupling_equiv = [6,4,8];
    case 3
      basecoupling_equiv = [7,4,8];
    case 4
      basecoupling_equiv = 8;
    case 5
      basecoupling_equiv = [1 4 8];
    case 6
      basecoupling_equiv = [2 4 8];
    case 7
      basecoupling_equiv = [3 4 8];
    case 8
      basecoupling_equiv = 4;
    case 9
      basecoupling_equiv = [4,8];
  end
  % Code für die eigentlich gesuchte Darstellung vorne anstellen (zuerst
  % suchen)
  basecoupling_equiv = [Coupling(1), basecoupling_equiv]; %#ok<AGROW>
  for k = basecoupling_equiv
    fcn_dir1 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Legs, ...
      sprintf('hd_G%dA0', k)); % keine P-Nummer für A0-Code
    if exist(fcn_dir1, 'file'), break; end % nehme den existierenden Ordner
  end
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