# Simulatore Caduta Massi (Rockfall Simulator)

Questo repository contiene un'applicazione Python per la simulazione della caduta massi, sviluppata con Matplotlib per la visualizzazione e PyQt5 per l'interfaccia grafica (GUI). Il simulatore permette di modellare il movimento di uno o più massi su un pendio, calcolando le loro traiettorie, velocità ed energie, sia in caduta libera che a contatto con la superficie. Tutto il codice è stato generato da Gemini AI con minimo intervento umano.

## Caratteristiche Principali

* **Simulazione Fisica Dettagliata**: Modella la caduta dei massi considerando gravità, attrito di rotolamento e resistenza dell'aria.
* **Gestione Collisioni**: Implementa un modello di collisione con il pendio basato su coefficienti di restituzione normale e tangenziale.
* **Pendio Personalizzabile**: Possibilità di definire il profilo del pendio tramite un file CSV, specificando coordinate X, Y e i coefficienti di restituzione lungo il percorso.
* **Generazione Multi-Massa**: Simula la caduta di più massi contemporaneamente con parametri iniziali casuali (volume, densità, coefficienti di attrito e drag).
* **Visualizzazione Grafica Interattiva**: Mostra le traiettorie dei massi in tempo reale su un grafico matplotlib integrato nella GUI.
* **Tabella dei Risultati**: Presenta i dati finali di ciascun masso (posizione, altezza, energia cinetica finale) in una tabella riassuntiva.
* **Esportazione Dati**: Consente di salvare i dati di output dettagliati di tutte le simulazioni in un file CSV per analisi successive.
* **Interfaccia Utente Intuitiva (GUI)**: Facile da usare grazie a un'interfaccia grafica ben organizzata che permette di impostare i parametri e controllare la simulazione.
* **Criteri di Arresto Flessibili**: I massi si fermano se la loro velocità scende sotto una soglia, se escono dai limiti definiti del pendio o se il tempo massimo di simulazione viene raggiunto.
* **Generazione Report PDF**: Possibilità di generare un report PDF con grafici e dati riassuntivi della simulazione.
* **Grafico a Torta Stati Finali**: Visualizzazione della distribuzione degli stati finali dei massi (es. fermo, uscito dall'area, in movimento).

## Requisiti

Per eseguire il simulatore, è necessario avere installato Python e le seguenti librerie:

* `numpy`
* `matplotlib`
* `PyQt5`
* `scipy`
* `pandas`

## Puoi installare queste librerie tramite `pip`:

```bash
pip install numpy matplotlib pyqt5 scipy pandas\
```
## Simulatore di Caduta Massi
Questo progetto implementa un simulatore di caduta massi con interfaccia grafica (GUI) basato su PyQt5, che permette di modellare il movimento di blocchi rocciosi su un pendio. Il simulatore calcola le traiettorie, le energie e gli stati finali dei massi, considerando vari parametri fisici e la possibilità di caricare profili di pendio personalizzati da file CSV.

Caratteristiche Principali
Simulazione Fisica Dettagliata: Modella la caduta dei massi considerando gravità, attrito di rotolamento e resistenza dell'aria.

Gestione Collisioni: Implementa un modello di collisione con il pendio basato su coefficienti di restituzione normale e tangenziale.

Pendio Personalizzabile: Possibilità di definire il profilo del pendio tramite un file CSV, specificando coordinate X, Y e i coefficienti di restituzione lungo il percorso.

Generazione Multi-Massa: Simula la caduta di più massi contemporaneamente con parametri iniziali casuali (volume, densità, coefficienti di attrito e drag).

Visualizzazione Grafica Interattiva: Mostra le traiettorie dei massi in tempo reale su un grafico matplotlib integrato nella GUI.

Tabella dei Risultati: Presenta i dati finali di ciascun masso (posizione, altezza, energia cinetica finale) in una tabella riassuntiva.

Esportazione Dati: Consente di salvare i dati di output dettagliati di tutte le simulazioni in un file CSV per analisi successive.

Interfaccia Utente Intuitiva (GUI): Facile da usare grazie a un'interfaccia grafica ben organizzata che permette di impostare i parametri e controllare la simulazione.

Criteri di Arresto Flessibili: I massi si fermano se la loro velocità scende sotto una soglia, se escono dai limiti definiti del pendio o se il tempo massimo di simulazione viene raggiunto.

Requisiti
Python 3.x

numpy

matplotlib

PyQt5

scipy

pandas

Installazione
Clona il repository (o scarica il file simulatore_masso_06.py):

Bash
```
git clone https://github.com/tuo_utente/nome_repo.git
cd nome_repo
```
(Sostituisci tuo_utente/nome_repo con il percorso effettivo del tuo repository su GitHub)

##Crea un ambiente virtuale (raccomandato):

Bash

```
python -m venv venv
```
# Su Windowshe public, the best way to achieve this is to make it

```.\venv\Scripts\activate```
# Su macOS/Linux
```source venv/bin/activate```
## Installa le dipendenze:


Bash
```
pip install numpy matplotlib PyQt5 scipy pandas
```
Utilizzo
Per avviare l'applicazione, esegui il file Python:


Bash
```
python simulatore_masso_06.py
```
## Interfaccia Utente:

Pannello Parametri di Simulazione: Consente di modificare i parametri globali della simulazione come gravità, passo temporale, numero di massi, intervalli di volume e densità, ecc.

## Pulsanti di Controllo:

Avvia Simulazione: Inizia la simulazione con i parametri correnti.

Reset Simulazione: Resetta la simulazione e i suoi parametri ai valori di default.

Carica Pendio da CSV: Permette di importare un profilo di pendio personalizzato da un file CSV. Il CSV deve contenere le colonne x, y, r (coefficiente di restituzione normale) e et (coefficiente di restituzione tangenziale).

Salva Output CSV: Esporta tutti i dati di output della simulazione in un file CSV.

Genera Report PDF: Genera un report PDF con i grafici e i dati riassuntivi.

Genera Grafico Torta Stati Finali: Mostra un grafico a torta della distribuzione degli stati finali dei massi.

Visualizzazione Grafica: Mostra le traiettorie dei massi in tempo reale.

Tabella Risultati: Visualizza i risultati finali per ciascun masso simulato.

Tabella Pendio: Mostra i punti discreti del pendio caricato da CSV o generato di default.

## Struttura del Codice
BloccoMasso: Classe che rappresenta un singolo masso con le sue proprietà fisiche e stato.

SimulatoreCadutaMassi: Classe principale per la logica fisica della simulazione, la gestione dei parametri e l'aggiornamento dello stato dei massi.

MainWindow: Classe che gestisce l'interfaccia utente grafica (GUI) usando PyQt5, collegando gli input dell'utente alla logica di simulazione e visualizzando i risultati.

## Esempio di File CSV per il Pendio
Il file CSV per il pendio deve avere il seguente formato (intestazioni incluse):

Snippet di codice
```
x,y,r,et
0.0,50.0,0.7,0.8
10.0,45.0,0.68,0.78
20.0,40.0,0.65,0.75
30.0,35.0,0.62,0.72
40.0,30.0,0.6,0.7
50.0,25.0,0.58,0.68
60.0,20.0,0.55,0.65
70.0,18.0,0.55,0.65
80.0,15.0,0.55,0.65
90.0,10.0,0.55,0.65
100.0,5.0,0.55,0.65
```
Dove:

x: Coordinata X del punto sul pendio.

y: Coordinata Y del punto sul pendio.

r: Coefficiente di restituzione normale in quel punto (tra 0 e 1).

et: Coefficiente di restituzione tangenziale in quel punto (tra 0 e 1).

 
