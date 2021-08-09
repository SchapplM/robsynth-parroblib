% Setzen oder Lösen der Schreibsperre für die PKM-Datenbank.
% Hiermit wird verhindert, dass parallele Prozesse gleichzeitig Daten ver-
% ändern oder ein schreibender einen lesenden Prozess stört.
% Erzeugt eine Sperr-Datei (.lock) und entfernt sie wieder.
% 
% Eingabe:
% mode
%   Modus der Schreibsperre: 'free', 'check', 'lock'
% content
%   Ressourcen, die geschützt werden: 'csv', PKM-Kinematik (z.B. 'P3PRRRR7')
%   'csv' beinhaltet alle Tabellen für die betroffenen Freiheitsgrade
%   PKM betrifft das kompilieren der Funktionen für die jeweilige PKM.
%   'template' beinhaltet den Temp-Ordner für die Generierung aus Vorlagen
% EEFG
%   EE-FG des Roboters im Format [1 1 1 0 0 0] (Bool)
% patience
%   Zeitdauer, bis die Schreibsperre ignoriert und trotzdem geschrieben
%   wird. Kann sinnvoll sein, falls der sperrende Parallelprozess
%   abgestürzt ist. Bezieht sich auf den Zeitpunkt der letzten Sperrung.
% verbosity
%   Bei true werden Textausgaben gemacht.

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-08
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function parroblib_writelock(mode, content, EEFG, patience, verbosity)

assert(isa(EEFG, 'logical') && all(size(EEFG)==[1 6]), 'EEFG muss 1x6 logical sein');
if nargin < 4
  patience = inf; % Warte unendlich, bis Ressource frei ist. Bei vorherigem Programm-Absturz ewig.
else
  % Es wird immer bis zu 20% der Zeit zufällig draufgeschlagen. Wenn viele
  % Instanzen parallel starten, warten sonst alle auf die blockierte
  % Ressourcen und geben gleichzeitig das Warten auf
  patience = patience + 0.2*patience*rand(1,1);
end
if nargin < 5
  verbosity = false;
end
%% Initialisierung
folder = sprintf('sym_%dT%dR', sum(EEFG(1:3)), sum(EEFG(4:6)));
repopath = fileparts(which('parroblib_path_init.m'));
% Prüfe, ob die Eingabe ein Robotername ist.
expression_robname = 'P(\d)([RP]+)(\d+)[V]?';
[is_robname, ~] = regexp(content,expression_robname);
if strcmp(content, 'csv')
  filename = 'writelock_csv.lock';
  filedir = fullfile(repopath, folder);
  ressource_name = sprintf('%s/%dT%dR', content, sum(EEFG(1:3)), sum(EEFG(4:6)));
elseif strcmp(content, 'template')
  filename = 'writelock_template.lock';
  filedir = fullfile(repopath, 'template_functions');
  ressource_name = 'Template';
elseif ~isempty(is_robname)
  expression_robname2 = '(.*)[G]';
  [robname_kin, ~] = regexp(content,expression_robname2, 'tokens','match');
  if ~isempty(robname_kin)
    pkm_kin = robname_kin{1}{1};
  else
    pkm_kin = content;
  end
  filename = 'writelock_mex.lock';
  filedir = fullfile(repopath, folder, pkm_kin);
  ressource_name = pkm_kin;
else
  error('Fall content=%s nicht definiert', content);
end
lockfile = fullfile(filedir, filename);
%% Sperrung aufheben (falls gefordert)
if strcmp(mode, 'free')
  if verbosity
    fprintf('Ressource %s ist wieder frei.\n', ressource_name);
  end
  delete(lockfile);
  return
end
%% Sperrung prüfen
t_start = tic();
while true
  fid=fopen(lockfile, 'r');
  if fid == -1
    break; % Sperrdatei existiert nicht. Ressource ist frei.
  end
  % Lese aus, wann die Sperre zuletzt gesetzt wurde
  try
    l = textscan(fid, '%d-%d-%d %d:%d:%d');
    locktime = datenum(sprintf('%04d-%02d-%02d %02d:%02d:%02d',l{1},l{2},l{3},l{4},l{5},l{6}));
  catch
    % Fehler: Ursache entweder nicht erwarteter Inhalt oder Verschwinden
    % der gerade noch vorhandenen Datei.
    warning('Fehler beim Interpretieren der Sperrdatei %s', lockfile);
    locktime = now()+1; % Dadurch wird weiter gewartet.
  end
  fclose(fid);
  % Sperrdatei existiert. Warte ab, bis Sperre vorbei ist
  if (now()-locktime)*24*3600 > patience
    if verbosity
      fprintf(['Es wurde %1.1fs auf die Ressource %s gewartet. ', ...
        'Ignoriere Schreibsperre mit Zeitstempel %s.\n'], ...
        toc(t_start),ressource_name, datestr(locktime,'yyyy-mm-dd HH:MM:SS'));
    end
    break;
  elseif verbosity
    fprintf('Ressource %s ist gesperrt (Zeitstempel %s). Warte ab bis sie frei wird.\n', ...
      ressource_name, datestr(locktime,'yyyy-mm-dd HH:MM:SS'));
  end
  pause(2+10*rand(1,1)); % Wartezeit bis zur nächten Prüfung.
end
if strcmp(mode, 'check')
  % Nur prüfen, ob zum Schreiben gesperrt, wenn nicht. Erfolg.
  return
end
%% Neue Sperr-Datei anlegen
% Hoffe, dass nicht hierzwischen ein anderer Prozess das gleiche gemacht
% hat (keine saubere Synchronisation)
fid = fopen(lockfile, 'w');
if fid == -1
  warning('Sperrdatei %s konnte nicht geschrieben werden', lockfile);
  return
end
fprintf(fid, '%s', datestr(now,'yyyy-mm-dd HH:MM:SS'));
fclose(fid);
if verbosity
  fprintf('Ressource %s wurde gesperrt.\n', ressource_name);
end