% Teste die Funktionen des ParRob-MdlBib-Repos durch hinzufügen bekannter Roboter
% Zusätzlich wird mit der HybrDyn-Toolbox symbolisch Code generiert.
% Die einzelnen Zellen dieses Skriptes sollten manuell ausgeführt werden

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12
% (C) Institut für Mechatronische Systeme, Universität Hannover

clear
clc

%% Roboter erstellen (Bsp 1: 3RRR)

NLEG = 3;
LEG_Names = {'S3RRR1'};
Actuation = {1, 1, 1}; % erstes Gelenk jeder Kette
EEdof0 = [1 1 0 0 0 1];
Coupling = [1 1]; % Basis-KS der Beinketten mit z-Achse nach oben
[Name, new] = parroblib_add_robot(NLEG, LEG_Names, Actuation, Coupling, EEdof0);

parroblib_generate_mapleinput({Name})
parroblib_generate_code({Name})

%% Roboter erstellen (Bsp 2: 3RPR)
NLEG = 3;
LEG_Names = {'S3RPR1'};
Actuation = {2, 2, 2}; % zweites Gelenk jeder Kette
EEdof0 = [1 1 0 0 0 1];
Coupling = [1 1]; % Basis-KS der Beinketten mit z-Achse nach oben
[Name, new] = parroblib_add_robot(NLEG, LEG_Names, Actuation, Coupling, EEdof0);
parroblib_generate_mapleinput({Name})
parroblib_generate_code({Name})
% Zum Testen: Robotermodell vorher entfernen
% parroblib_remove_robot('P3RPR1G1P1A1')
% parroblib_remove_robot('P3RPR1G1P1A2')

%% Roboter erstellen (Bsp 3: 6UPS)
NLEG = 6;
LEG_Names = {'S6RRPRRR14V3'};
Actuation = {3, 3, 3, 3, 3, 3}; % drittes Gelenk jeder Kette
EEdof0 = [1 1 1 1 1 1];
Coupling = [1 1]; % Basis-KS zeigt nach oben. Funktioniert für U-Gelenk
[Name, new] = parroblib_add_robot(NLEG, LEG_Names, Actuation, Coupling, EEdof0);
parroblib_generate_mapleinput({Name})
parroblib_generate_code({Name})

%% Roboter Bsp 4: 4PUU
% P4PRRRR1G3P1A1
NLEG = 4;
LEG_Names = {'S5PRRRR1'};
Actuation = {1,1,1,1};
EEdof0 = [1 1 1 0 0 1];
Coupling = [3 1]; % Basis-KS zeigt nach oben. Funktioniert für U-Gelenk
[Name, new] = parroblib_add_robot(NLEG, LEG_Names, Actuation, Coupling, EEdof0);
parroblib_generate_mapleinput({Name})
parroblib_generate_code({Name})

%% Alle räumlichen PKM erstellen, deren Beinketten in Kugel- oder Kardangelenk enden
EEFG_Ges = [1 1 1 0 0 0; ...
            1 1 1 0 0 1; ...
            1 1 1 1 1 1];
for i_FG = 1:size(EEFG_Ges,1)
  EEdof0 = EEFG_Ges(i_FG,:);
  [PNames_Kin, PNames_Akt] = parroblib_filter_robots(sum(EEdof0), EEdof0, EEdof0);
  for j = 1:length(PNames_Akt)
    fprintf('%d/%d: %s\n', j, length(PNames_Akt), PNames_Akt{j});
    [NLEG, LEG_Names, Actuation, ActNr, symrob, EE_dof0, PName_Kin] = parroblib_load_robot(PNames_Akt{j});
    parroblib_generate_mapleinput({PNames_Akt{j}})
    parroblib_generate_code({PNames_Akt{j}})
  end
end

%% Alle 6FG-PKM erstellen, deren Beinketten mit Kugelgelenk enden
% Zulässige PKM sollten vorher durch die Struktur- und Maßsynthese erzeugt
% werden. Das Hinzufügen aller möglicher ungeprüfter PKM ist nicht sinnvoll
EEdof0 = [1 1 1 1 1 1];
[PNames_Kin, PNames_Akt] = parroblib_filter_robots(6, EEdof0, EEdof0);
for j = 1:length(PNames_Akt)
  fprintf('%d/%d: %s\n', j, length(PNames_Akt), PNames_Akt{j});
  [NLEG, LEG_Names, Actuation, ActNr, symrob, EE_dof0, PName_Kin] = parroblib_load_robot(PNames_Akt{j});
  parroblib_generate_mapleinput({PNames_Akt{j}})
  parroblib_generate_code({PNames_Akt{j}})
end