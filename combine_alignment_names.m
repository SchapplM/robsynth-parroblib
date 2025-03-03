% Kombiniere die Ausdrücke der Koppelgelenk-Varianten und gebe einen Text aus
% 
% Eingabe:
% GP_list
%   Array von GP-Kombinationen, die funktionieren. 
%   Bspw. [1 1]; [2 1]; [3 1]]
% 
% Ausgabe:
% GP_compr_str
%   Zeichenkette mit Kombination der Werte.
%   Bspw. '(v, t, r)-v'
% GP_str
%   Nicht komprimierte Zeichenkette.
%   Bspw. 'v-v, t-v, r-v

% Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2022-11
% (C) Institut für Mechatronische Systeme, Leibniz Universität Hannover

function [GP_compr_str, GP_str] = combine_alignment_names(GP_list)
GP_compr_str = '';
GP_str = '';
if isempty(GP_list), return; end
% Erzeuge einen einzelnen Ausdruck für die GP-Varianten
% durch multiplikative Verknüpfung symbolischer Variablen
tmp = 0;
for i = 1:size(GP_list,1)
  tmp = tmp + sym(sprintf('G%d', GP_list(i,1))) * sym(sprintf('P%d', GP_list(i,2)));
end
% Nicht-vereinfachten Ausdruck ausgeben
GP_str = sym2str_alignment(tmp);
% Vereinfache den Ausdruck: Ausklammern öfter vorkommender Begriffe
% Benutze die bekannte Tatsache der Struktur der Terme aus Produkten von G_i*P_j
% Simplify-Befehl funktioniert hier nicht wie gewünscht. Daher spezifische Lösung.
P=sym('P',[1 max(GP_list(:,2))]);
tmp2 = simplify_expr(tmp, P);
assert(simplify(tmp2 - tmp)==0, 'Koeffizienten nicht vereinfachbar');
GP_compr_str = sym2str_alignment(tmp2);
end


function str = sym2str_alignment(symexp)
% Eingabe:
% symexp
%   Symbolischer Ausdruck mit Gi*Pj-Termen
% Ausgabe:
% str
%   Zeichenkette für Auflistung in der Ergebnis-Tabelle

% Ändere Reihenfolge der Terme. Zuerst G dann P
if contains(cell2str(string(symexp)), '+') % Summanden dafür bilden
  summands = children(symexp+1);
  % Prüfe rekursiven Aufruf von children. Wenn ein Doppel-Klammer-Ausdruck
  % auftaucht, werden andere Summanden manchmal durch children-Befehl
  % zusammengefasst.
  for i = 1:length(summands)
    sc = children(summands{i});
    for j = 1:length(sc)
      if sum(cell2str(string(sc{j}))== '*') == 1
        % Kind-Elemente auflösen. Muss sich um ein Produkt handeln und die
        % Kind-Elemente selbst daher um Summanden, die durch die vorherige
        % Zusammenfassung nicht gefunden wurden
        summands = {summands{[1:i-1, i+1:end]}, sc{1:end}};
        break; % hiernach keine weitere Aufteilung möglich
      end
    end
  end
  % Prüfe, ob durch die manuell korrigierte Bildung der Summanden der
  % Ausdruck gleich geblieben ist
  symexp_test = -1; % da oben eine 1 zur Hilfe addiert wurde
  for i = 1:length(summands)
    symexp_test = symexp_test + summands{i};
  end
  test_symexp = symexp_test - symexp;
  if test_symexp ~= 0
    error('Zusammenfassung der Summanden fehlgeschlagen');
  end
else
  summands = {symexp};
end
str = ''; % Ausdruck neu zusammensetzen
for i = 1:length(summands)
  % Suche nach den einzelnen Faktoren
  % Alternative 1: Regulärer Ausdruck
  summands_str_i = cell2str(string(summands{i}));
  if strcmp(summands_str_i, '1'), continue; end % Hilfs-Term
  [~, match] = regexp(summands_str_i, '([\(\)A-Z0-9\s+])*', 'tokens', 'match');
  % Alternative 2: Faktorenzerlegung symbolisch (vergisst die Klammer)
  % factors = children(summands{i});
  % match = string(factors);
  if length(match) ~= 2
    error('Format passt nicht'); 
  end
  if contains(match{1}, 'P') % Reihenfolge tauschen
    summands_str_i = sprintf('%s*%s', match{2}, match{1});
  else
    summands_str_i = sprintf('%s*%s', match{1}, match{2});
  end
  if ~isempty(str) % Anhängen der Summanden mit Komma+Leerzeichen
    % Verhindere Zeilenumbruch in einem Ausdruck mit mbox
    summands_str_i = [', \mbox{', summands_str_i, '}'];
  else
    summands_str_i = ['\mbox{', summands_str_i, '}'];
  end
  str = [str, summands_str_i];
  % In der Klammer ohne Leerzeichen nach Komma aufzählen
  str = strrep(str, ' + ', ',');
end

% Sortiere die Reihenfolge der Terme so, dass erst die G-Faktoren und dann
% die P-Faktoren kommen
str = strrep(str, '*', '-');
% Übersetzung der Koppelgelenk-Varianten in Buchstaben
% Siehe align_base_coupling und align_platform_coupling.
base_str = ['v', 't', 'r', 'c', 'V', 'T', 'R', 'C', 'p', 'v'];
plf_str = ['v', 't', 'r', 'V', 'T', 'R', 'm', 'p', 'c', 'i'];
for i = length(base_str):-1:1 % umgekehrte Reihenfolge, sonst wird 10 mit 1 teilweise ersetzt
  str = strrep(str, sprintf('G%d', i), base_str(i));
end
for i = length(plf_str):-1:1
  str = strrep(str, sprintf('P%d', i), plf_str(i));
end
end

function Enew = simplify_expr(E,P)
% Vereinfache den Term
% https://de.mathworks.com/matlabcentral/answers/1858453-reduce-length-of-sum-of-products-of-symbolic-variables#comment_2478688
[co,terms] = coeffs(E,P);
uniqueco = unique(co);
Enew = 0;
index = 1:numel(co);
for ii = 1:numel(uniqueco)
    Enew = Enew + uniqueco(ii)*sum(terms(index(isAlways(uniqueco(ii) == co,'unknown','false'))));
end
end