% Gebe an, wie viele Strukturen in der Datenbank enthalten sind

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2021-03
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

clc
clear

% serroblib_gen_bitarrays
% parroblib_gen_bitarrays

% Seriellroboter-Datenbank auslesen
serroblibpath=fileparts(which('serroblib_path_init.m'));
SerRob_DB_all = load(fullfile(serroblibpath, 'serrob_list.mat'));
% Indizes in Seriellroboter-Datenbank auf Varianten
I_SR_var = SerRob_DB_all.AdditionalInfo(:,2)==1;
% SerRob_DB_all.Names(I_SR_var)

% Gehe alle möglichen PKM-Freiheitsgrade durch
EEFG_Ges = [1 1 0 0 0 1; ...
            1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 0; ...
            1 1 1 1 1 1];
for i_FG = 1:size(EEFG_Ges,1)
  EE_FG = logical(EEFG_Ges(i_FG,:));
  fprintf('\nPrüfe Anzahl der PKM mit EE-FG %dT%dR\n', sum(EE_FG(1:3)), sum(EE_FG(4:6)));
  % Lese PKM-Datenbank aus
  [PNames_Kin, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots(EE_FG, 0);
  % Stelle verschiedene Eigenschaften der PKM zusammen. Zunächst bezogen
  % auf Kinematiken (Mechanismen) ohne Betrachtung der Aktuierung.
  expression = 'P(\d)([RP]+)(\d+)[V]?(\d*)G(\d*)P(\d*)';
  [tokens, ~] = regexp(PNames_Kin,expression,'tokens','match');
  LegTechJointString = cell(length(PNames_Kin),1);
  LegJointString = cell(length(PNames_Kin),1);
  LegJointNums = NaN(length(PNames_Kin),1);
  LegChains = cell(length(PNames_Kin),1);
  LegChainIndex = NaN(length(PNames_Kin),1);
  LegNumTechJoints = NaN(length(PNames_Kin),1);
  LegChainIsVariant = false(length(PNames_Kin),1);
  Couplings = NaN(length(PNames_Kin),2);
  PNames_Kin_TechJoints = cell(length(PNames_Kin),1);
  for k = 1:length(PNames_Kin)
    token_k = tokens{k};
    LegJointString{k} = token_k{1}{2};
    LegJointNums(k) = length(LegJointString{k});
    LegChains{k} = ['S',sprintf('%d',LegJointNums(k)),LegJointString{k},token_k{1}{3}];
    if ~isempty(token_k{1}{4})
      LegChains{k} = [LegChains{k}, 'V', token_k{1}{4}];
      LegChainIsVariant(k) = true;
    end
    I_SR_k = strcmp(SerRob_DB_all.Names, LegChains{k});
    assert(sum(I_SR_k)==1, sprintf('Kein eindeutiger Eintrag %s in SerRobLib', LegChains{k}));
    LegChainIndex(k) = find(I_SR_k);
    LegNumTechJoints(k) = SerRob_DB_all.AdditionalInfo(I_SR_k,6);
    Couplings(k,:) = [str2double(token_k{1}{5}),str2double(token_k{1}{6})];
    SName_TechJoint = fliplr(regexprep(num2str(SerRob_DB_all.AdditionalInfo(I_SR_k,7)), ...
      {'1','2','3','4','5'}, {'R','P','C','U','S'}));
    LegTechJointString{k} = SName_TechJoint;
    PNames_Kin_TechJoints{k} = [token_k{1}{1}, LegTechJointString{k}];
  end
  % Bestimme die Zuordnung zwischen Aktuierungs-Liste und Kinematik-Liste
  I_KinAct = NaN(length(PNames_Akt), 1); % Zähl-Indizes von Akt. in Kin.
  for k = 1:length(PNames_Akt)
    % Finde die Kinematik der aktuierten PKM
    PName_k = PNames_Akt{k};
    I_k = find(strcmp(PNames_Kin, PName_k(1:end-2)));
    assert(length(I_k)==1, sprintf('Kein eindeutiger Eintrag %s in Kinematik-DB', PNames_Akt{k}));
    I_KinAct(k) = I_k;
  end
  % Bestimme die Anzahl der möglichen Beinketten-Gelenke  
  LegNums_possible = [];
  for N = unique(LegJointNums)' % Gehe die Anzahl der Beingelenke durch
    % Bestimme alle prinzipiell möglichen Beinketten
    I_SR_N = SerRob_DB_all.N == N; % Indizes in SerRobLib mit dieser Gelenkzahl
    LegNums_possible = [LegNums_possible, sum(~I_SR_var&I_SR_N)]; %#ok<AGROW>
    fprintf('%d allgemeine Beinketten mit %d Gelenken\n', sum(~I_SR_var&I_SR_N), N);
    fprintf('%d Varianten-Beinketten mit %d Gelenken\n', sum(I_SR_var&I_SR_N), N);
    % Bestimme die tatsächlich genutzten Beinketten
    [~,I_unique] = unique(LegChains);
    I_SR_unique = false(length(I_SR_var),1);
    I_SR_unique(LegChainIndex(I_unique)) = true;
    % SerRob_DB_all.Names(I_SR_unique)
    fprintf('%d allgemeine Beinketten mit %d Gelenken tatsächlich benutzt\n', ...
      sum(I_SR_unique&~I_SR_var&I_SR_N), N);
    % SerRob_DB_all.Names(I_SR_unique&I_SR_var&I_SR_N)
    fprintf('%d Varianten-Beinketten mit %d Gelenken tatsächlich benutzt\n', ...
      sum(I_SR_unique&I_SR_var&I_SR_N), N);
    % SerRob_DB_all.Names(I_SR_unique&~I_SR_var&I_SR_N)
  end
  % Bestimme die möglichen und genutzten Koppelgelenk-Anordnungen
  Couplings_unique = unique(Couplings, 'rows');
  CouplingG_unique = unique(Couplings(:,1));
  CouplingP_unique = unique(Couplings(:,2));
  fprintf('%d Kombinationen von Koppelgelenkanordnungen: %d Gestell ([%s]), %d Plattform ([%s])\n', ...
    size(Couplings_unique,1), length(CouplingG_unique), disp_array(CouplingG_unique','%d'), ...
    length(CouplingP_unique), disp_array(CouplingP_unique','%d'));
  % Mögliche Nummer des aktuierten Gelenks
  n_akt = 2;
  if all(EE_FG == [1 1 1 1 1 1])
    n_akt = 3;
  end
  n_possib_max = length(CouplingG_unique)*length(CouplingP_unique)*sum(LegNums_possible)*n_akt;
  fprintf('Anzahl theoretische Kombinationen (allgemein): %d x %d x %d x %d = %d\n', ...
    length(CouplingG_unique), length(CouplingP_unique), sum(LegNums_possible), n_akt, n_possib_max);
  
  % Gebe die möglichen Beinketten an
  fprintf(['Insgesamt %d kinematisch unterschiedliche Mechanismen mit ', ...
    'allgemeinen Beinketten\n'], length(PNames_Kin(~LegChainIsVariant)));
  fprintf(['Insgesamt %d kinematisch unterschiedliche Mechanismen mit ', ...
    'Varianten-Beinketten\n'], length(PNames_Kin(LegChainIsVariant)));
  
  % Prüfe, welche Beinketten in den aktuierbaren PKM tatsächlich benutzt
  % werden
  LegNums_possible_Akt = [];
  for N = unique(LegJointNums(I_KinAct))'
    I_SR_N = SerRob_DB_all.N == N;
    LegNums_possible_Akt = [LegNums_possible_Akt, sum(~I_SR_var&I_SR_N)]; %#ok<AGROW>
    % Bestimme die verwendeten Beinketten
    I_unique = unique(I_KinAct); % Index bezogen auf PKM-Kinematiken (die aber auch zu einer Aktuierbarkeit führen)
    I_SR_unique = false(length(I_SR_var),1); % Index bezogen auf Beinketten
    I_SR_unique(LegChainIndex(I_unique)) = true;
    % SerRob_DB_all.Names(I_SR_unique)
    assert(length(unique(SerRob_DB_all.Names(I_SR_unique)))==sum(I_SR_unique), ...
      'Auswahl der seriellen Beinketten nicht eindeutig');
    
    fprintf('%d allgemeine Beinketten mit %d Gelenken tatsächlich benutzt (in steuerbaren PKM)\n', ...
      sum(I_SR_unique&~I_SR_var&I_SR_N), N);
    % SerRob_DB_all.Names(I_SR_unique&I_SR_var&I_SR_N)
    fprintf('%d Varianten-Beinketten mit %d Gelenken tatsächlich benutzt (in steuerbaren PKM)\n', ...
      sum(I_SR_unique&I_SR_var&I_SR_N), N);
    % SerRob_DB_all.Names(I_SR_unique&~I_SR_var&I_SR_N)
  end
  
  % Filtere die Kinematiken (ohne Aktuierung) nach Anzahl technischer
  % Gelenke
  for ntj = min(LegNumTechJoints):max(LegNumTechJoints)
    I_kin_ntj = LegNumTechJoints == ntj;
    fprintf('%d Mechanismen mit %d technischen Gelenken pro Beinkette\n', sum(I_kin_ntj), ntj);
    % PNames_Kin_TechJoints(I_kin_ntj)
    % PNames_Kin(I_kin_ntj)
  end
  % Gebe die Anzahl der allgemein aktuierbaren PKM an
  I_act_var = LegChainIsVariant(I_KinAct)==1;
  % Filtere dabei die Nummer des aktuierten Gelenks (gezählte technische
  % Gelenke). Gebe nur die PKM an, die eine halbwegs gestellnahe Aktuierung
  % haben.
  I_act_lowacttechjoint = AdditionalInfo_Akt(:,3) <= n_akt;
  fprintf('%d PKM mit allgemeinen Beinketten\n', sum(~I_act_var&I_act_lowacttechjoint));
  fprintf('%d PKM mit Varianten-Beinketten\n', sum(I_act_var&I_act_lowacttechjoint));
  % PNames_Akt(~I_act_var&I_act_lowactjoint)
  
  % Anzahl der technischen Gelenke bezogen auf die aktuierten PKM
  LegNumTechJoints_Akt = LegNumTechJoints(I_KinAct);
  % Filtere die PKM  nach Anzahl technischer Gelenke ("bis zu")
  for ntjmax = min(LegNumTechJoints):max(LegNumTechJoints)
    I_act_ntjmax = LegNumTechJoints_Akt <= ntjmax;
    fprintf('%d PKM mit max. %d technischen Gelenken\n', ...
      sum(I_act_ntjmax&I_act_lowacttechjoint), ntjmax);
  end
  % Filtere die PKM  nach Anzahl technischer Gelenke ("genaue Anzahl")
  for ntj = min(LegNumTechJoints):max(LegNumTechJoints)
    I_act_ntj = LegNumTechJoints_Akt == ntj;
    fprintf(['%d PKM mit genau %d technischen Gelenken (mit Koppelgelenk- ', ...
      'und Aktuierungs-Variationen)\n\t%s\n'], sum(I_act_ntj&I_act_lowacttechjoint), ntj, ...
      disp_array(unique(PNames_Kin_TechJoints(I_KinAct(I_act_ntj&I_act_lowacttechjoint)))'));
    % disp(PNames_Akt(I_act_ntj&I_act_lowacttechjoint)');
  end
end