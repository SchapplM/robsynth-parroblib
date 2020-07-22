% Erstelle die Vorlagen-Funktionen aller Robotermodelle in der Bibliothek
% 
% Eingabe:
% Names
%   Cell array mit Namen der zu erstellenden Robotermodelle
%   Optional: Wenn nicht angegeben, werden alle erstellt.
% skip_existing
%   true: Bereits existierende tpl-Ordner überspringen (Standard)
%   false: Alle Dateien immer neu erstellen
% mex_results
%   true: Die aus Vorlagen generierten Funktionen werden zum testen direkt
%   kompiliert.
% 
% Siehe auch: serroblib_create_template_functions
% 
% TODO: Wenn es mehrere PKM mit gleichen Beinketten, aber verschiedenen
% Koppelpunkten in der Datenbank gibt, werden die selben
% Vorlagen-Funktionen (für IK) mehrfach direkt nacheinander erstellt.

% Junnan Li, Hiwi bei Moritz Schappler, 2020-02
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-04
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_create_template_functions(Names, skip_existing, mex_results)

%% Initialisierung
old_dir = pwd();
% Prüfe Eingabeargument
repopath=fileparts(which('parroblib_path_init.m'));
if nargin < 1 || isempty(Names) % keine Eingabe einer Liste. Nehme alle.
  % Stelle Liste aller Roboter zusammen
  Names = {};  
  for N = 3:6
    [PNames_Kin, ~, ~] = parroblib_filter_robots(N, ones(1,6), zeros(1,6), 6); 
    Names = {Names{:}, PNames_Kin{:}}; %#ok<CCAT>
  end
end
if nargin < 2
  skip_existing = true;
end
if nargin < 3
  mex_results = false;
end

%% Alle Vorlagen-Funktionen aus HybrDyn-Repo kopieren
% Alle Funktionen dieser Liste werden roboterspezifisch aus der Liste
% erstellt. Die Anpassung sind nur geringfügig und ermöglichen Kompilierung
function_list_copy_robotics = {...
  {'kinematics', 'fkineEE_traj.m'},...
  {'kinematics', 'invkin.m'},...
  {'kinematics', 'invkin_traj.m'},...
  {'kinematics', 'constr1D_trans.m'},...
  {'kinematics', 'constr4D_rot.m'},...
  {'kinematics', 'constr4D_traj.m'},...
  {'kinematics', 'constr1grad_rt.m'},...
  {'kinematics', 'constr1grad_tq.m'},...
  {'kinematics', 'constr1grad_tr.m'},...
  {'kinematics', 'constr1grad_tt.m'},...
  {'kinematics', 'constr1gradD_tq.m'},...
  {'kinematics', 'constr1gradD_tr.m'},...
  {'kinematics', 'constr1gradD_tt.m'},...
  {'kinematics', 'constr2_trans.m'},...
  {'kinematics', 'constr2grad_q.m'},...
  {'kinematics', 'constr2grad_rq.m'},...
  {'kinematics', 'constr2grad_rr.m'},...
  {'kinematics', 'constr2grad_rt.m'},...
  {'kinematics', 'constr2grad_tq.m'},...
  {'kinematics', 'constr2grad_tr.m'},...
  {'kinematics', 'constr2grad_tt.m'},...
  {'kinematics', 'constr2gradD_tq.m'},...
  {'kinematics', 'constr2gradD_tr.m'},...
  {'kinematics', 'constr2gradD_tt.m'},...
  {'kinematics', 'constr2grad_x.m'},...
  {'kinematics', 'constr3_rot.m'},...
  {'kinematics', 'constr3grad_q.m'},...
  {'kinematics', 'constr3.m'},...
  {'kinematics', 'constr3grad_rq.m'},...
  {'kinematics', 'constr3grad_rr.m'},...
  {'kinematics', 'constr3grad_rt.m'},...
  {'kinematics', 'constr3grad_x.m'},...
  {'kinematics', 'constr3gradD_q.m'},...
  {'kinematics', 'constr3gradD_rq.m'},...
  {'kinematics', 'constr3gradD_rr.m'},...
  {'kinematics', 'constr3gradD_x.m'},...
  {'kinematics', 'invkin3.m'}, ...
  {'kinematics', 'constr4grad_q.m'},...
  {'kinematics', 'constr4grad_rq.m'},...
  {'kinematics', 'constr4grad_rr.m'},...
  {'kinematics', 'constr4grad_x.m'},...
  {'kinematics', 'constr4gradD_q.m'},...
  {'kinematics', 'constr4gradD_rq.m'},...
  {'kinematics', 'constr4gradD_rr.m'},...
  {'kinematics', 'constr4gradD_x.m'}};

mkdirs(fullfile(repopath, sprintf('template_functions')))
rtp = fileparts(which('robotics_toolbox_path_init.m'));
if isempty(rtp)
  warning('Die Robotik-Toolbox muss im Pfad sein (siehe README.MD)');
  return
end
for tmp = function_list_copy_robotics
  tplf = tmp{1};
  file1 = fullfile(rtp, tplf{1}, ['pkm_',tplf{2},'.template']);
  if ~exist(file1, 'file')
    warning('Die Vorlagen-Datei %s existiert nicht. Vermutlich Inkonsistenz der Repo-Versionen', tplf{2});
    continue
  end
  file2 = fullfile(repopath, 'template_functions');
  % Prüfe nur lesend, ob die Dateien identisch sind. Nur dann in
  % Sammel-Ordner kopieren. Das verringert die Gefahr von Dateikonflikten.
  data1 = dir(file1);
  data2 = dir(fullfile(file2, ['pkm_',tplf{2},'.template']));
  if isempty(data2) || ... % Die Zieldatei existiert noch nicht
      (data1.datenum ~= data2.datenum || data1.bytes ~= data2.bytes) % beide Dateien sind nicht identisch
    copyfile(file1, file2); % Nur in diesem Fall die Vorlagen im Ordner anpassen
  end
end

% Generiere die Liste der Vorlagen-Funktionen aus den tatsächlich
% existierenden Dateien. Dadurch können durch Benutzer hinzugefügte
% Funktionen auch generiert werden
fl_tmp = dir(fullfile(repopath, 'template_functions', '*.template'));
function_list = cell(1,length(fl_tmp));
for i = 1:length(fl_tmp)
  % Präfix "robot_" und Suffix ".template" entfernen
  function_list{i} = strrep(fl_tmp(i).name(5:end), '.template', '');
end
%% Gehe alle Modellnamen durch und erstelle alle Vorlagen-Funktionen
for i = 1:length(Names)
  Name_i = Names{i};
  % Daten zur PKM laden
  [N, LEG_Names,~,~,~,~,~,~,PName_Legs] = parroblib_load_robot([Name_i,'A0']);
  RS = serroblib_create_robot_class(LEG_Names{1});
  fcn_dir = fullfile(repopath, sprintf('sym%dleg', N), PName_Legs, 'tpl');
  % Platzhalter-Ausdrücke für diese PKM. Werden aus ParRob-Klassen
  % generiert. Ersetze Ausdruck später von Spalte 1 mit Ausdruck in Spalte 2
  subsexp_array = { ...
     'PN', PName_Legs; ...
     'SN', RS.mdlname;...
     'NJ', sprintf('%d',N * RS.NQJ);... % bezogen auf PKM
     'NLEG', sprintf('%d',N);...
     'NKP', sprintf('%d',length(RS.pkin_gen));...
     'NQJ', sprintf('%d',RS.NQJ);... % bezogen auf serielle Beinkette
     'VERSIONINFO', sprintf(['Generated in ParRobLib from ', ...
                  'template on %s'], datestr(now,'yyyy-mm-dd HH:MM:SS'))};

  % Kopiere alle Vorlagen-Funktionen an die Ziel-Orte und Ersetze die
  % Platzhalter-Ausdrücke.
  % Prüfe, ob alle Matlab-Funktionen oder Mex-Dateien da sind
  mkdirs(fcn_dir);  % Ordner existiert vielleicht noch nicht. Neu erstellen.
  cd(fcn_dir); % In Ordner wechseln für kürzeren sed-Befehl (und zum Finden der Dateien)
  function_list_mex = {};
  num_files_written = 0;
  for tmp = function_list
    tplf = tmp{1};
    file1=fullfile(repopath, 'template_functions', ['pkm_',tplf,'.template']);
    file2=fullfile(fcn_dir, [PName_Legs,'_',tplf]);
    function_list_mex = [function_list_mex(:); [PName_Legs,'_',tplf(1:end-2)]];
    if skip_existing && exist(file2, 'file')
      continue % Keine Datei erzeugen. Nur mex-Liste erstellen.
    end
    fid1 = fopen(file1);
    fid2 = fopen(file2, 'w');
    tline = fgetl(fid1);
    while ischar(tline)
      % Zeile weiterverarbeiten: Platzhalter-Ausdrücke ersetzen
      for ii = 1:size(subsexp_array,1)
        tline = strrep(tline, ['%',subsexp_array{ii,1},'%'], subsexp_array{ii,2});
      end
      fwrite(fid2, [tline, newline()]); % Zeile in Zieldatei schreiben
      tline = fgetl(fid1); % nächste Zeile
    end
    fclose(fid1);fclose(fid2);
    
    % Prüfe, ob die Erstellung erfolgreich war (falls die Code-Generierung
    % nicht vollständig durchgeführt wurde, gehen nicht alle Funktionen)
    fid2 = fopen(file2, 'r');
    tline = fgetl(fid2);
    function_invalid = false;
    while ischar(tline)
      if contains(tline, 'NOTDEFINED') % noch nicht implementiert
        function_invalid = true;
        break;
      end
      tline = fgetl(fid2);
    end
    fclose(fid2);
    if function_invalid
      fprintf('%d/%d: Datei %s konnte nicht erzeugt werden (Einige Variablen nicht definiert)\n', ...
        i, length(Names), tplf);
      delete(file2);
      continue;
    end
    num_files_written = num_files_written + 1;
  end
  fprintf('%d/%d: Vorlagen-Funktionen für %s erstellt (%d/%d).\n', i, ...
    length(Names), Name_i, num_files_written, length(function_list));
  
  % Testen: Kompilieren aller erzeugter Funktionen im Zielordner
  if mex_results
    cd(fcn_dir);
    matlabfcn2mex(function_list_mex);
  end
end
% Zurückwechseln in vorheriges Verzeichnis
cd(old_dir);
