% Generiere Matlab-Code mit Maple-Dynamik-Toolbox für PKM-Roboterstrukturen
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen
% force_par
%   Erzwinge Neu-Kompilierung, egal ob bereits generierter Code (für die
%   PKM) vorliegt 
% force_ser
%   Erzwinge Neu-Kompilierung der seriellen Kette
% 
% Vorher: 
% * Funktion maplerepo_path.m muss vorliegen mit Rückgabe des
%   Repo-Pfades der Maple-Dynamik-Toolbox ("HybrDyn")
% * Maple-Eingabedaten müssen für die Roboterstruktur mit
%   parroblib_generate_mapleinput.m erzeugt werden
% 
% Siehe auch: serroblib_generate_code.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_generate_code(Names, force_par, force_ser)
%% Init
if nargin < 2
  force_par = false;
end
if nargin < 3
  force_ser = false;
end

repopath=fileparts(which('parroblib_path_init.m'));
%% Roboterstrukturen durchgehen
for i = 1:length(Names)
  n = Names{i};

  % Daten über den Roboter zusammenstellen
  [NLEG, LEG_Names, ~, ActNr] = parroblib_load_robot(n);
  % Robotereigenschaften aus dem Namen auslesen.
  % TODO: Einbindung nicht-symmetrischer PKM
  expression = 'P(\d)([RP]+)(\d+)A(\d+)'; % Format "P3RRR1A1"
  [tokens, ~] = regexp(n,expression,'tokens','match');
  res = tokens{1};
  PName_Kin = ['P', res{1}, res{2}, res{3}];
  
  % Pfad zur Maple-Dynamik-Toolbox (muss im Repo abgelegt werden; s.o.)
  mrp = maplerepo_path();
  
  % Maple-Toolbox-Eingabe laden (wurde an anderer Stelle erzeugt)
  % (durch parroblib_generate_mapleinput.m)
  mapleinputfile=fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, ...
    sprintf('hd_A%d',ActNr), sprintf('robot_env_par_%s', n));
  if ~exist(mapleinputfile, 'file')
    error('Datei %s existiert nicht. Wurde `parroblib_generate_mapleinput.m` ausgeführt?', fileparts(mapleinputfile) );
  end
  % Verzeichnisse für die zu erzeugenden Matlab-Funktionen
  outputdir_tb_par = fullfile(mrp, 'codeexport', n, 'matlabfcn'); % Verzeichnis in der Maple-Toolbox
  outputdir_local = fileparts(mapleinputfile); % Verzeichnis in der Bibliothek
  
  % Prüfe, ob Code schon einmal generiert wurde 
  % (und im Zielverzeichnis vorliegt)
  if ~force_par && length(dir(fullfile(outputdir_local, '*.m'))) > 5
    % das werden wohl schon genug .m-Dateien sein.
    continue
  end
  
  % Code-Erstellung für serielle Beinkette starten, falls diese nicht
  % vorliegt
  for k = 1:length(LEG_Names)
    % Suche nach temporären Dateien im Arbeitsverzeichnis von HybrDyn
    tmpdir_tb_ser = fullfile(mrp, 'codeexport', LEG_Names{k}, 'tmp'); % tmp-Verzeichnis in der Maple-Toolbox
    if ~force_ser && length(dir(fullfile(tmpdir_tb_ser, '*_maple.m'))) > 20
      % Die serielle Beinkette wurde wahrscheinlich schon generiert.
      continue
    end
    % Generiere diese Beinkette neu (ohne Rückkopieren der
    % Matlab-Funktionen ins SerRobLib-Repo)
    serroblib_generate_code({LEG_Names{k}}, true, true)  
  end
  
  % Eingabedatei für parallelen Roboter kopieren
  copyfile( mapleinputfile, fullfile(mrp, 'robot_codegen_definitions', 'robot_env_par') );
  
  % Code-Erstellung für parallelen Roboter starten (ohne Generierung der
  % Beinketten; das wurde oben schon gemacht).
  fprintf('Starte Code-Generierung %d/%d für %s\n', i, length(Names), n);
  system( sprintf('cd %s && ./robot_codegen_start.sh --fixb_only --notest --parrob --not_gen_serial', mrp) ); %  > /dev/null
  
  % generierten Code zurückkopieren (alle .m-Dateien)
  for f = dir(fullfile(outputdir_tb_par, '*.m'))'
    copyfile(fullfile(outputdir_tb_par, f.name), fullfile(outputdir_local, f.name));
  end
end