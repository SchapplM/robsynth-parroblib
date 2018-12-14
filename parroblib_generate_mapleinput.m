% Generiere die Eingabedateien zur Code-Generierung für gegebene
% Roboterstrukturen
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen, für die der Code
%   erzeugt werden soll. Der Name entspricht dem Schema "PxRRRyyyAzzz" mit
%   x=Anzahl Beine, yyy laufende Nummer für "RRR" (Benennung der seriellen
%   Beinkette) und zzz Nummer der Aktuierung.
% 
% Siehe auch: serroblib_generate_mapleinput.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_generate_mapleinput(Names)
%% Init
repopath=fileparts(which('parroblib_path_init.m'));
%% Roboterstrukturen durchgehen
for i = 1:length(Names)
  n = Names{i};
  
  %% Daten für diesen Roboter laden
  [NLEG, LEG_Names, Actuation, ActNr, ~, EE_FG0] = parroblib_load_robot(n);
  % Robotereigenschaften aus dem Namen auslesen.
  % TODO: Einbindung nicht-symmetrischer PKM
  expression = 'P(\d)([RP]+)(\d+)A(\d+)'; % Format "P3RRR1A1"
  [tokens, ~] = regexp(n,expression,'tokens','match');
  res = tokens{1};
  PName_Kin = ['P', res{1}, res{2}, res{3}];
  %% Maple-Toolbox-Eingabe erzeugen
  % Zur Definition des Formats der Eingabedatei; siehe HybrDyn-Repo
  mapleinputfile=fullfile(repopath, sprintf('sym%dleg', NLEG), PName_Kin, ...
    sprintf('hd_A%d',ActNr), sprintf('robot_env_par_%s', n));
  mkdirs(fileparts(mapleinputfile));
  fid = fopen(mapleinputfile, 'w');
  
  % Allgemeine Definitionen
  fprintf(fid, 'robot_name := "%s":\n', n);
  fprintf(fid, 'leg_name := "%s":\n', LEG_Names{1});
  fprintf(fid, 'AKTIV := Matrix(%d,1,[',NLEG);
  for j = 1:NLEG
    if length(Actuation{j}) > 1
      error('Beinketten mit mehrfacher Aktuierung noch nicht unterstützt');
    end
    if j > 1, fprintf(fid, ', '); end
    fprintf(fid, '%d', Actuation{j});
  end
  fprintf(fid, ']);\n');
  fprintf(fid, 'N_LEGS := %d:\n', NLEG);

  % TODO: Anzahl der Parameter der Basis-Koppelpunkte abhängig von den
  % EE-FG. nur die FG bearbeiten, die für EE-FG relevant sind.
  % Aktuelle Annahme: Gestellt ist eine Ebene, alle Beinketten liegen mit
  % Basis in dieser Ebene
  % TODO: xA und yA scheinen für die Dynamik gar nicht relevant zu sein.
  fprintf(fid, 'leg_frame := Matrix(7,1,[xA, yA, 0, 0, 0, gammaLeg, xyz]):\n');
  
  % EE-FG
  fprintf(fid, 'xE_s := Matrix(7,1,[');
  k=0;
  for j = 1:6
    if j > 1, fprintf(fid, ', '); end
    if EE_FG0(j) == 0
      fprintf(fid, '0');
    else
      k=k+1;
      fprintf(fid, 'x_all[%d]', k);
    end
  end
  fprintf(fid, ', xyz]):\n');
  fclose(fid);
end