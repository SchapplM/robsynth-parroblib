% Ändere die Eigenschaften der Spalte für die Werte des zusätzlichen
% Winkels.
% Vorher: Nur für einen Winkel (ersten theta-Parameter einstellbar)
% Jetzt: Allgemeinere Form für beliebig viele freie theta/alpha-Parameter
% 
% Siehe auch: add_csv_column_theta1_values.m

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2020-08
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clear
clc

%% Liste von Seriellen Beinketten zusammenstellen
serroblibpath=fileparts(which('serroblib_path_init.m'));
Anz_Beinketten = zeros(1,7);
LegChainList = cell(0,2);
for i = 3:6
  mdllistfile_Ndof = fullfile(serroblibpath, sprintf('mdl_%ddof', i), sprintf('S%d_list.mat',i));
  l = load(mdllistfile_Ndof);
  for j = 1:length(l.Names_Ndof)
    csvline = serroblib_bits2csvline(l.BitArrays_Ndof(j,:));
    I_alpha = contains(csvline, 'alpha');
    I_theta = contains(csvline, 'theta');
    % Mögliche Kombinationen zur Synthese: Winkel 0, 90, egal.
    % Mögliche Speicherung in Datenbank mit Zuständen 0, 90, 0+90, egal
    nparam=sum(I_alpha|I_theta);
    ncomb=3^nparam;
    Anz_Beinketten(nparam+1) = Anz_Beinketten(nparam+1)+1;
    if sum(nparam) == 0
      continue;
    end
    fprintf('%d FG, Rob. %d/%d (%s): %d Freie Parameter, also %d Kombinationen\n', ...
      i, j, length(l.Names_Ndof), l.Names_Ndof{j}, sum(I_alpha|I_theta), ncomb);
    % disp([csvline(I_alpha),csvline(I_theta)]);
    paramstr = '';
    params = [csvline(I_alpha),csvline(I_theta)];
    for k = 1:length(params)
      paramstr = [paramstr, params{k}]; %#ok<AGROW>
      if k < length(params)
        paramstr=[paramstr,',']; %#ok<AGROW>
      end
    end
    % Merke, welche Beinkette welche freien Parameter hat
    LegChainList = [LegChainList; {l.Names_Ndof{j}, paramstr}]; %#ok<AGROW>
  end
end
for i =0:length(Anz_Beinketten)-1
  fprintf('%d Beinketten mit %d freien Parametern\n', Anz_Beinketten(i+1), i);
end

%% PKM-Datenbank bearbeiten
repopath=fileparts(which('parroblib_path_init.m'));
actfilelist = dir(fullfile(repopath, '**', 'actuation.csv'));

for i = 1:length(actfilelist)
  acttabfile = fullfile(actfilelist(i).folder, actfilelist(i).name);
  % Finde den Roboternamen heraus
  [tokens,match] = regexp(actfilelist(i).folder, '.*/(P[\d]+.*)', 'tokens', 'match');
  PName = tokens{1}{1};
  SName = ['S',PName(2),PName(3:end)];
  % Prüfe den Status der Variablen
  I = strcmp(LegChainList(:,1), SName);
  if ~any(I)
    fprintf('Überspringe %s\n', acttabfile);
    continue
  elseif sum(I) > 1
    error('Name der Beinkette nicht eindeutig');
  end
  paramstr = LegChainList{I,2};
  fprintf('Bearbeite %s\n', acttabfile);
  acttabfile_copy = [acttabfile, '.copy']; % Kopie der Tabelle zur Bearbeitung
  fid = fopen(acttabfile, 'r');
  fidc = fopen(acttabfile_copy, 'w');
  tline = fgetl(fid);
  iline=0;
  abort_migration = false;
  while ischar(tline)
    iline=iline+1;
    if iline == 1
      if ~contains(tline, 'Wertebereich Winkel1')
        warning('Überschrift in Datei %s nicht im alten Format. Überspringe.', acttabfile);
        abort_migration = true;
        break;
      end
      % Überschrift umbenennen
      tline_new = strrep(tline, 'Wertebereich Winkel1', ['Wertebereich ',paramstr]);
    else
      % Verändere die Daten in der Spalte
      csvline = strsplit(tline, ';', 'CollapseDelimiters', false);
      data_angle1 = csvline{end};
      if strcmp(data_angle1, '*')
        csvline{end} = 'a'; % "Arbitrary"
      elseif strcmp(data_angle1, '0')
        csvline{end} = 'p'; % "Parallel"
      elseif strcmp(data_angle1, '90')
        csvline{end} = '0'; % "Orthogonal"
      elseif strcmp(data_angle1, '090')
        csvline{end} = 'b'; % "Both" ("Parallel" und "Orthogonal")
      elseif isempty(data_angle1)
        csvline{end} = ''; % Nicht definiert
      else
        error('Unerwarteter Wert in Tabellenspalte: %s', data_angle1);
      end
      tline_new = csvline{1};
      for k = 2:length(csvline)
        tline_new = sprintf('%s;%s', tline_new, csvline{k});
      end
    end
    fwrite(fidc, [tline_new, newline]);
    tline = fgetl(fid); % nächste Zeile
  end
  fclose(fid);
  fclose(fidc);
  if ~abort_migration
    copyfile(acttabfile_copy, acttabfile);
  end
  delete(acttabfile_copy);
end