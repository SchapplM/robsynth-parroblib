% Finde PKM-Robotermodelle in der Datenbank mit gegebenen Filtern
% 
% Eingabe:
% NLEG
%   Anzahl der Beine der symmetrischen PKM
% EE_FG [1x6]
%   Vektor mit 1/0 für Belegung ob EE-FG aktiv ist
% EE_FG_Mask [1x6]
%   Maske, die festlegt ob die FG exakt wie in `EE_FG` sind, oder ob auch
%   gesperrte FG wirklich nicht beweglich sind
% max_rankdeficit [1x1]
%   Grad des möglichen Rangverlustes der PKM-Jacobi-Matrix (Standard: 0)
%   Wird ein Wert von 1 bis 6 (max.) eingesetzt, werden auch PKM
%   ausgegeben, deren gewählte Aktuierung keine Bewegung ergibt.
% 
% Rückgabe:
% PNames_Kin
%   Cell-Array mit Namen aller Roboterkinematiken ohne Hinweis auf
%   Aktuierung
% PNames_Akt
%   Cell-Array aller Roboter mit Aktuierung als Namensbestandteil.
%   Diese Namen können durch die Funktion parroblib_load_robot abgerufen
%   werden
% AdditionalInfo_Akt
%   Array mit zusätzlichen Infos für alle Strukturen aus PNames_Akt (in den Zeilen).
%   Spalten:
%   1: Rangverlust der Jacobi-Matrix (in den vorgesehenen FG der PKM)
% 
% TODO: Aktuell sind nur symmetrische PKM berücksichtigt.
% 
% Siehe auch: serroblib_filter_robots

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2019-01
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function [PNames_Kin, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots(NLEG, EE_FG0, EE_FG_Mask, max_rankdeficit)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
if nargin < 4
  max_rankdeficit = 0;
end
%% Kinematik-Tabelle durchsuchen
PNames_Kin = {};
PNames_Akt = {};
AdditionalInfo_Akt = [];
kintabfile=fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list_kin.mat', NLEG));
if ~exist(kintabfile, 'file')
  error('Datei %s existiert nicht. parroblib_gen_bitarrays ausführen!', kintabfile);
end
tmp = load(kintabfile);
KinTab = tmp.KinTab;
% Tabelle zeilenweise durchgehen
for i = 1:size(KinTab,1)
  % Gültige Tabellenzeile. Daten auslesen
  PName_Kin = KinTab.Name{i};
  LegName = KinTab.Beinkette{i}; %#ok<NASGU>
  EE_FG_i = KinTab.EEFG(i,:);
  % Daten filtern
  if ~all( (EE_FG_i & EE_FG_Mask) == (EE_FG0 & EE_FG_Mask))
    continue
  end
  % Roboter passt zu Filter. Auswahl in Liste
  PNames_Kin = {PNames_Kin{:}, PName_Kin}; %#ok<CCAT>
end
if nargout == 1
  return
end
%% Aktuierung heraussuchen
acttabfile=fullfile(repopath, sprintf('sym%dleg', NLEG), sprintf('sym%dleg_list_act.mat', NLEG));
if ~exist(acttabfile, 'file')
  error('Datei %s existiert nicht. parroblib_gen_bitarrays ausführen!', acttabfile);
end
tmp = load(acttabfile);
ActTab = tmp.ActTab;
% Wähle die Roboter in der Datenbank anhand der Filterkriterien aus
I_rankloss = ActTab.Rankloss_Platform <= max_rankdeficit;
% Wähle nur PKM, die auf den oben basierenden Grundstrukturen aufbauen (in
% der Datenbank mit 3-Beinketten-PKM stehen sonst 2T1R- und 3T0R-PKM)
I_kin = false(size(ActTab,1),1);
for i = 1:size(KinTab,1)
  PName_Kin = KinTab.Name{i};
  I_kin = I_kin | contains(ActTab.Name, PName_Kin);
end
% Wähle PKM, deren Rang noch nicht geprüft wurde (in der Tabelle steht ein
% "?"). Werden ausgegeben. Sonst funktioniert die Struktursynthese nicht.
I_unchecked = isnan(ActTab.Rankloss_Platform);
% Kombiniere Filterkriterien
I = (I_rankloss|I_unchecked)&I_kin;
% Ausgabevariablen zuweisen
PNames_Akt = ActTab.Name(I);
AdditionalInfo_Akt = ActTab.Rankloss_Platform(I);
