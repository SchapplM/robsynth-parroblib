% Formatiere einen Roboternamen für Publikationen (in Latex)
% 
% Eingabe:
% Name:
%   Name des Roboters
% Modus:
%   Ausgabeformat
%   1=Latex, Akzente über Buchstaben
%   2=Technische Gelenke der Beinkette ("UPS" statt "RRPRRR")
%   3=Gelenkkette ohne Aktuierung ("RRPRRR")

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-11
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function FormatName = parroblib_format_robot_name(RobName, Modus)
FormatName = '';
%% Daten für Kinematik laden
repopath=fileparts(which('parroblib_path_init.m'));
[NLEG, LEG_Names, Actuation, ~, ~, ~, EE_dof0, ...
  PName_Kin, ~, ~, ~, JointParallelity] = parroblib_load_robot(RobName, 2);
Chain_Name = LEG_Names{1};
NLegJ = str2double(Chain_Name(2));
%% Zeichenkette für Parallelität generieren. Siehe Kong/Gosselin 2007, S.10
if Modus == 1
  if isnan(JointParallelity)
    warning('Keine Parallelität der Gelenke in Datenbank angegeben. Namensbestimmung unmöglich.')
    return
  end
  % Variablen mit Latex-Code für Roboter-Namen
  Chain_StructName = '';
  Chain_StructNameAct = '';
  pgroups = JointParallelity;
  % Setze die Nummer 0 für Gelenke, die zu keiner parallelen Gruppe gehören
  for j = 1:NLegJ
    if sum(pgroups == pgroups(j)) == 1
      pgroups(j) = 0; % Dadurch dann kein Akzent auf dem Buchstaben
    end
  end
  % Entferne nicht belegte Nummern
  j = 1;
  while j <= max(pgroups)
    if ~any(pgroups==j) % reduziere alle folgenden Nummern um 1
      pgroups(pgroups>j) = pgroups(pgroups>j) - 1;
      continue; % Diese Nummer nochmal prüfen
    end
    j = j + 1;
  end
  % Setze die Markierungen entsprechend der Gruppen
  for j = 1:NLegJ
    groupidx = pgroups(j); % hochzählen
    if groupidx == 0
      % Diese Gelenkausrichtung gibt es nur einmal. Es muss kein
      % Gruppensymbol darüber gelegt werden
      newsymbol = '{';
    elseif groupidx == 1
      newsymbol = '{\`';
    elseif groupidx == 2
      newsymbol = '{\''';
    elseif groupidx == 3
      newsymbol = '{\=';
    else
      error('mehr als drei Achsrichtungen nicht vorgesehen');
    end
    % Füge "P"/"R" hinzu
    newsymbol = [newsymbol, Chain_Name(2+j), '}']; %#ok<AGROW>
    Chain_StructName = [Chain_StructName, newsymbol]; %#ok<AGROW>
    if any(Actuation{1} == j)
      Chain_StructNameAct = [Chain_StructNameAct, '\underline']; %#ok<AGROW>
    end
    Chain_StructNameAct = [Chain_StructNameAct, newsymbol]; %#ok<AGROW>
  end
  FormatName = sprintf('%d-%s', NLEG, Chain_StructNameAct);
end
%% Zeichenkette für technische Gelenke einer Beinkette
if Modus == 2
  EEstr = sprintf('%dT%dR', sum(EE_dof0(1:3)), sum(EE_dof0(4:6)));
  kintabmatfile = fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_kin.mat']);
  tmp = load(kintabmatfile);
  KinTab = tmp.KinTab;
  SName_TechJoint = KinTab.Beinkette_Tech{strcmp(KinTab.Name, PName_Kin)};
  FormatName = SName_TechJoint;
end
%% Zeichenkette für Gelenkkette des Beins ohne Aktuierung
if Modus ==3
  FormatName = Chain_Name(3:3+NLegJ-1);
end
%% Ende
if isempty(FormatName)
  error('Modus nicht implementiert');
end
