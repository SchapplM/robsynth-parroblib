% Finde PKM-Robotermodelle in der Datenbank mit gegebenen Filtern
% 
% Eingabe:
% EE_FG [1x6]
%   Vektor mit 1/0 für Belegung ob EE-FG aktiv ist
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

function [PNames_Kin, PNames_Akt, AdditionalInfo_Akt] = parroblib_filter_robots(EE_FG0, max_rankdeficit)

%% Initialisierung
repopath=fileparts(which('parroblib_path_init.m'));
assert(size(EE_FG0,2)==6, 'EE_FG0 muss 6 Spalten haben');
if nargin < 2
  max_rankdeficit = 0;
end
EEstr = sprintf('%dT%dR', sum(EE_FG0(1:3)), sum(EE_FG0(4:6)));
%% Kinematik-Tabelle durchsuchen
PNames_Kin = {};
PNames_Akt = {};
AdditionalInfo_Akt = [];
kintabfile=fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_kin.mat']);
if ~exist(kintabfile, 'file')
  error('Datei %s existiert nicht. parroblib_gen_bitarrays ausführen!', kintabfile);
end
tmp = load(kintabfile);
KinTab = tmp.KinTab;
PNames_Kin = KinTab.Name';

if nargout == 1
  return
end
%% Aktuierung heraussuchen
acttabfile=fullfile(repopath, ['sym_', EEstr], ['sym_',EEstr,'_list_act.mat']);
if ~exist(acttabfile, 'file')
  error('Datei %s existiert nicht. parroblib_gen_bitarrays ausführen!', acttabfile);
end
tmp = load(acttabfile);
ActTab = tmp.ActTab;
% Wähle die Roboter in der Datenbank anhand der Filterkriterien aus
I_rankloss = ActTab.Rankloss_Platform <= max_rankdeficit;
% Wähle PKM, deren Rang noch nicht geprüft wurde (in der Tabelle steht ein
% "?"). Werden ausgegeben. Sonst funktioniert die Struktursynthese nicht.
I_unchecked = isnan(ActTab.Rankloss_Platform);
% Kombiniere Filterkriterien
I = I_rankloss|I_unchecked;
% Ausgabevariablen zuweisen
PNames_Akt = ActTab.Name(I);
AdditionalInfo_Akt = ActTab.Rankloss_Platform(I);
