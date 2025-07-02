# Simulatore Caduta Massi (Rockfall Simulator)

Questo repository contiene un'applicazione Python per la simulazione della caduta massi, sviluppata con Matplotlib per la visualizzazione e PyQt5 per l'interfaccia grafica (GUI). Il simulatore permette di modellare il movimento di uno o più massi su un pendio, calcolando le loro traiettorie, velocità ed energie, sia in caduta libera che a contatto con la superficie.

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

Puoi installare queste librerie tramite `pip`:

```bash
pip install numpy matplotlib pyqt5 scipy pandas

