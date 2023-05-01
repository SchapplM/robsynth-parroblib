# Bibliothek paralleler Roboter

Dieses Repo enthält eine Sammlung von Modellen paralleler Roboter. Die Modelle beinhalten die Beschreibung funktionierender Systeme und teilweise symbolisch generierte Kinematik und Dynamik der Roboter.

Moritz Schappler, moritz.schappler@imes.uni-hannover.de, 2018-12  
(C) Institut für Mechatronische Systeme, Universität Hannover

## Benutzung

Für die Benutzung des Repos muss mit dem Skript `parroblib_path_init.m` aus dem Hauptordner die Pfadinitialisierung durchgeführt werden.

### Abhängigkeiten

Die Robotermodelle setzen die Toolbox imes-matlab, serrob_mdlbib und imes-matlab-robotics mit initialisiertem Pfad voraus.
Für die Code-Generierung muss die Datei `maplerepo_path.m` aus der Vorlage `maplerepo_path.m.template` so erstellt werden, dass diese auf den Pfad zur 
Toolbox robot-dynamics-hybriddyn zeigt (diese muss ebenfalls heruntergeladen werden).

## Inhalt

Die Robotermodelle sind nach der Anzahl der Beine Roboter aufgeteilt und in csv-Dateien tabellarisch aufgelistet.
Zum Eintragen, Auslesen und Filtern von Robotermodellen existieren Funktionen.

Im Folgenden wird unterschieden zwischen

* Grundstruktur-Name: Kennzeichnet die Anzahl und Art der Beinketten (nicht aber ihre Orientierung)
* Kinematik-Name: Kennzeichnet die Anzahl und Art der Beinketten sowie die Art ihrer Befestigung an Gestell und Plattform.
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
* `parroblib_create_template_functions.m`: Matlab-Funktionen aus Vorlagen erzeugen. Die Funktionen sind kompilierbar und schneller als in der Matlab-Klasse implementierte Methoden
* `parroblib_update_template_functions.m`: Aktualisiere die Vorlagen-Funktionen auf die neuste Version

Datenbank ordnen:
* `parroblib_remove_robot`: Entfernen eines Roboters aus der Datenbank. Das ist nützlich, um im nachhinein unbrauchbare Roboter zu entfernen oder um Roboter zu Testzwecken temporär hinzuzufügen.

### Übersichtsdateien

Alle Informationen der Roboterbibliothek liegen in csv-Dateien (Trennzeichen: Semikolon) vor. Dadurch ist eine automatische Verarbeitung und eine Versionsverwaltung möglich.

#### Liste der Robotermodelle für eine Gelenkfolge

Die Datei `sym_aTbR/sym_aTbR_list.csv` enthält eine Liste aller symmetrischer Robotermodelle mit a translatorischen und b rotatorischen Plattform-FG.
Der Kinematik-Name des Roboters ergibt sich aus

* P: Kennzeichnung für parallele Roboter,
* x: Anzahl der Beine
* R/P die Kennzeichnung für Dreh- oder Schubgelenk.
* y: Für laufende Nummer der Kinematiken mit dieser Reihenfolge von Dreh- und Schubgelenken.
* V: Kennung für Variante der Kinematik der Beinkette (optional, siehe SerRobLib)
* zzz: laufende Nummer der Variante (optional)
* G: Kennung für Gestell-Koppelgelenkorientierung
* u: Nummer dazu
* P: Kennung für Plattform-Koppelgelenkorientierung
* v: Nummer dazu

Der Teil RRPRRyyy(Vzzz) des Namens entspricht den Systemen aus der Bibliothek serieller Roboter.

#### Liste der möglichen Aktuierungen für ein Robotermodell

Die Datei `sym_aTbR/PxRRRRRRyyyyy/actuation.csv` enthält mögliche Aktuierungen für die Gelenke der Grundstruktur des Robotermodells. Jede kinematische Struktur kann unterschiedliche Aktuierungen haben, die dann einen anderen Roboter ergibt.
In der Tabelle sind die Namen der aktuierten Roboter abgelegt, die den Namen der Kinematik mit dem Zusatz A1, A2, ... tragen.

#### Liste der Roboterparameter für ein Robotermodell (models.csv)

Die Datei `sym_aTbR/PxRRRRRRyyyyy/models.csv` enthält mögliche Roboterparameter für ein parametriertes Robotermodell. 
Die Spaltenüberschriften entsprechen den einzelnen Parametern aus der Struktur- und Maßsynthese.

### Matlab-Code

Matlab-Funktionen für jeden Roboter liegen im Ordner `sym_aTbR/PxRRRRRRyyyyy/hd_GuuPvvAzz`.
Wobei die Dynamik-Funktionen auf die Plattform bezogen im Ordner A0 liegen, in dem die Aktuierung der Gelenk noch keine Rolle spielt. In den Ordnern mit Kennung der Aktuierung (A1, A2, ...) liegen Matlab-Funktionen, die sich auf diese eine mögliche Aktuierung beziehen.
