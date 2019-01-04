# Bibliothek paralleler Roboter

Dieses Repo enthält eine Sammlung von Modellen paralleler Roboter. Die Modelle beinhalten symbolisch generierte Kinematik und Dynamik der Roboter.

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12  
(C) Institut für Mechatronische Systeme, Universität Hannover

## Benutzung

Für die Benutzung des Repos muss mit dem Skript `parroblib_path_init.m` aus dem Hauptordner die Pfadinitialisierung durchgeführt werden.

### Abhängigkeiten

Die Robotermodelle setzen die Toolbox imes-matlab, serrob_mdlbib und imes-matlab-robotics mit initialisiertem Pfad voraus.

## Inhalt

Die Robotermodelle sind nach der Anzahl der Beine Roboter aufgeteilt und in csv-Dateien tabellarisch aufgelistet.
Zum Eintragen, Auslesen und Filtern von Robotermodellen existieren Funktionen.

Im Folgenden wird unterschieden zwischen

* Kinematik-Name: Kennzeichnet die Anzahl und Art der Beinketten.
* Aktuierungs-Name: Kinematik-Name mit Zusatz "A" und laufender Nummer, mit der die Aktuierung in einer Tabelle gefunden werden kann.

### Funktionen

Die genauen Übergabeargumente und Funktionsweisen sind dem Quelltext der Funktionen und Beispielaufrufen zu entnehmen.
Die Funktionen orientieren sich inhaltlich und von der Benennung an denen für die Datenbank serieller Roboter

Datenbank aufbauen (z.B. mit Struktursynthese-Ergebnissen):
* `parroblib_add_robot`: Hinzufügen eines neuen Roboters zur Datenbank. Der Roboter wird durch die Anzahl und die Namen der Beinketten (aus der Datenbank serieller Roboter) und durch die Aktuierung gekennzeichnet

Datenbank auslesen und benutzen:
* `parroblib_filter_robots`: Liste von Robotern nach Filterkriterien (wie Anzahl und Art der FG) ausgeben
* `parroblib_load_robot.m`: Eigenschaften eines Roboters ausgeben
* `parroblib_create_robot_class`: Instanz einer Matlab-Klasse für einen Roboter erstellen. Damit sind dann weitere Berechnungen mit dem Modell des Roboters möglich
* `parroblib_generate_mapleinput`: Eingabedaten zur Code-Generierung für diesen Roboter erstellen
* `parroblib_generate_code`: Kinematik- und Dynamikfunktionen symbolisch mit der HybrDyn-Toolbox berechnen lassen.
* `parroblib_addtopath`: Matlab-Funktionen dieses Roboters zum Matlab-Pfad hinzufügen, damit die Funktionsaufrufe funktionieren.

### Übersichtsdateien

Alle Informationen der Roboterbibliothek liegen in csv-Dateien (Trennzeichen: Semikolon) vor. Dadurch ist eine automatische Verarbeitung und eine Versionsverwaltung möglich.

#### Liste der Robotermodelle für eine Gelenkfolge

Die Datei `symxleg/symxleg_list.csv` enthält eine Liste aller symmetrischer Robotermodelle mit x Beinen.
Der Kinematik-Name des Roboters ergibt sich aus

* P: Kennzeichnung für parallele Roboter,
* x: Anzahl der Beine
* R/P die Kennzeichnung für Dreh- oder Schubgelenk.
* y: Für laufende Nummer der Kinematiken mit dieser Reihenfolge von Dreh- und Schubgelenken.

Der letzte Teil RRPRRyyy des Namens entspricht den Systemen aus der Bibliothek serieller Roboter.

#### Liste der möglichen Aktuierungen für ein Robotermodell

Die Datei `symxleg/PxRRRRRRyyyyy/actuation.csv` enthält mögliche Aktuierungen für die Gelenke dieses Robotermodells. Jede kinematische Struktur kann unterschiedliche Aktuierungen haben, die dann einen anderen Roboter ergibt.
In der Tabelle sind die Namen der aktuierten Roboter abgelegt, die den Namen der Kinematik mit dem Zusatz A1, A2, ... tragen.

### Matlab-Code

Matlab-Funktionen für jeden Roboter liegen im Ordner `symxleg/PxRRRRRRyyyyy/hd_Azz`.
Wobei die Dynamik-Funktionen auf die Plattform bezogen im Ordner A0 liegen, in dem die Aktuierung der Gelenk noch keine Rolle spielt. In den Ordnern mit Kennung der Aktuierung (A1, A2, ...) liegen Matlab-Funktionen, die sich auf diese eine mögliche Aktuierung beziehen.
