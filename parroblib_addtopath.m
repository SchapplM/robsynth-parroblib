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
  [NLEG, ~, ~, ActNr, ~, ~, PName_Kin] = parroblib_load_robot(Name);

  % Pfade für Matlab-Funktionen hinzufügen
  fcn_dir1 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Kin, 'hd_A0');
  fcn_dir2 = fullfile(parroblibpath, sprintf('sym%dleg', NLEG), PName_Kin, sprintf('hd_A%d', ActNr));
  if exist(fcn_dir1, 'file')
    addpath(fcn_dir1);
  end
  if exist(fcn_dir2, 'file')
    addpath(fcn_dir2);
  end
end