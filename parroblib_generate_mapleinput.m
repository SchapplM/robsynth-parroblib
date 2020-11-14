% Generiere die Eingabedateien zur Code-Generierung für gegebene
% Roboterstrukturen
% 
% Eingabe:
% Names
%   Cell-Array mit Liste der Namen der Roboterstrukturen, für die der Code
%   erzeugt werden soll. Der Name entspricht dem Schema "PxRRRyGuPvAz" mit
%   x=Anzahl Beine, y laufende Nummer für "RRR" (Benennung der seriellen
%   Beinkette), u Nummer der Gestellkoppelgelenk-Orientierung, v Nummer der
%   Plattformkoppelgelenk-Orientierung und zzz Nummer der Aktuierung.
% 
% Siehe auch: serroblib_generate_mapleinput.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

function parroblib_generate_mapleinput(Names)
%% Init
repopath=fileparts(which('parroblib_path_init.m'));
serrobpath=fileparts(which('serroblib_path_init.m'));
%% Roboterstrukturen durchgehen
for i = 1:length(Names)
  n = Names{i};
  
  %% Daten für diesen Roboter laden
  [NLEG, LEG_Names, Actuation, Coupling, ActNr, ~, EE_FG0, ~, PName_Legs, AdditionalInfo_Akt] = parroblib_load_robot(n);
  EEstr = sprintf('%dT%dR', sum(EE_FG0(1:3)), sum(EE_FG0(4:6)));
  for j = 1:NLEG
    if isempty(Actuation{j})
      error('Beinkette %d ist nicht aktuiert. Das wird aktuell nicht unterstützt', j);
    end
  end
  % Robotereigenschaften der Beinketten prüfen. Siehe serroblib_gen_bitarrays
  legdata = load(fullfile(serrobpath, sprintf('mdl_%sdof', LEG_Names{1}(2)), ...
    sprintf('S%s_list', LEG_Names{1}(2))));
  AddInfo_Leg = legdata.AdditionalInfo(strcmp(legdata.Names_Ndof,LEG_Names{1}),:);
  if all(EE_FG0==[1 1 0 0 0 1]) && AddInfo_Leg(1) > 2 || AddInfo_Leg(1) > 3
    warning('Symbolischer Code kann nicht basierend auf %s-Beinkette gebildet werden. Zu viele positionsbeeinflussende Gelenke (%d)', LEG_Names{1}, AddInfo_Leg(1));
    continue
  end
  % Robotereigenschaft der PKM prüfen.
  if AdditionalInfo_Akt(1) > 0 % Rangverlust
    warning('PKM hat laut Datenbank/Struktursynthese Rangverlust. Symbolischer Code nicht sinnvoll generierbar.');
    % continue
  end
  %% Definition der Beinketten-Orientierung
  % Bei Definition neuer Beinketten-Orientierungen in align_base_coupling
  % (ParRob-Klasse) müssen die Werte hier angepasst werden.
  leg_frame_entries = {'xA', 'yA', '0', '0', '0', 'gammaLeg'};
  if ~all(EE_FG0 == [1 1 0 0 0 1])
    if Coupling(1)==1 || Coupling(1)==5
      % Für Koppelpunkt-Methode 1 wird das Basis-KS der Beinkette nur um
      % die z-Achse gedreht. Lasse obige Standardeinstellung
    elseif Coupling(1)==2  || Coupling(1)==6
      leg_frame_entries = {'xA', 'yA', '0', '-Pi/2', 'betaLeg', '-Pi/2'};
    elseif Coupling(1)==3  || Coupling(1)==7
      leg_frame_entries = {'xA', 'yA', '0', '-Pi/2', 'betaLeg', '0'};
    elseif Coupling(1)==4  || Coupling(1)==8
      leg_frame_entries = {'xA', 'yA', '0', 'alphaLeg', 'betaLeg', 'gammaLeg'};
    else
      warning('Methode %d noch nicht für Code-Generierung definiert.', Coupling(1));
      continue;
    end
  end

  %% Maple-Toolbox-Eingabe erzeugen
  % Zur Definition des Formats der Eingabedatei; siehe HybrDyn-Repo
  mapleinputfile=fullfile(repopath, ['sym_', EEstr], PName_Legs, ...
    sprintf('hd_G%dP%dA%d',Coupling(1),Coupling(2),ActNr), sprintf('robot_env_par_%s', n));
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
  % Aktuelle Annahme: Gestell ist eine Ebene, alle Beinketten liegen mit
  % Basis in dieser Ebene
  % TODO: xA und yA scheinen für die Dynamik gar nicht relevant zu sein.
  fprintf(fid, 'leg_frame := Matrix(7,1,[%s, %s, %s, %s, %s, %s, xyz]):\n', ...
    leg_frame_entries{1}, leg_frame_entries{2}, leg_frame_entries{3}, ...
    leg_frame_entries{4}, leg_frame_entries{5}, leg_frame_entries{6});
  
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
  % Wähle standardmäßig den maximalen Optimierungsgrad.
  fprintf(fid, 'codegen_opt := 2:\n');
    
  
  % Generierung der Dynamik nur für erstes Aktuierungs-Modell jeder PKM:
  % Die Dynamik in Plattform-Koordinaten ist immer gleich und die
  % Dynamik-Berechnung in Antriebs-Koordinaten wird numerisch mit Jacobi
  % berechnet. Nur die Jacobi wird dann noch neu berechnet
  % Hier wird die Generierung standardmäßig deaktiviert und in
  % parroblib_generate_code.m einmal aktiviert und mit Kennung "A0"
  % gespeichert.
  fprintf(fid, 'codeexport_invdyn := false:\n');
  fclose(fid);
end
