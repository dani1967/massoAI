Come Eseguire il Simulatore
Clona il repository (o scarica il file simulatore_masso_06.py):

Bash

git clone [https://github.com/tuo_nome_utente/nome_del_tuo_repo.git](https://github.com/tuo_nome_utente/nome_del_tuo_repo.git)
cd nome_del_tuo_repo
(Sostituisci tuo_nome_utente e nome_del_tuo_repo con i tuoi dati reali su GitHub.)

Esegui il file simulatore_masso_06.py:

Bash

python simulatore_masso_06.py
Si aprirà l'interfaccia grafica del simulatore.

Utilizzo dell'Interfaccia Grafica (GUI)
Pannello Parametri di Simulazione: Modifica la gravità, il passo temporale della simulazione, il numero di massi da generare e gli intervalli per volume, densità, attrito di rotolamento e coefficiente di drag. Puoi anche impostare le coordinate di partenza dei massi e le forme (sferica, lenticolare, prismatica o casuale).

Carica Pendio da CSV: Clicca su questo pulsante per selezionare un file CSV contenente i dati del pendio (colonne x, y, r, et). r e et sono i coefficienti di restituzione normale e tangenziale rispettivamente. Se non carichi un file, verrà utilizzato un pendio inclinato di default.

Controllo Simulazione: Utilizza i pulsanti "Avvia Simulazione", "Reset Simulazione" per controllare il flusso della simulazione.

Salva Output CSV: Dopo una simulazione, puoi salvare tutti i dati dettagliati (posizione, velocità, energia cinetica nel tempo per ogni masso) in un file CSV per analisi esterne.

Genera Report PDF: Crea un report in formato PDF con i grafici delle traiettorie e una sintesi dei risultati.

Genera Grafico Torta Stati Finali: Visualizza un grafico a torta che mostra la percentuale di massi che si sono fermati, sono usciti dall'area di simulazione, o sono ancora in movimento alla fine della simulazione.

Struttura del File CSV per il Pendio
Il file CSV per il profilo del pendio deve essere formattato come segue, con le prime due righe che possono essere commentate con # (non sono obbligatorie, ma utili per la documentazione):

Snippet di codice

# Descrizione del profilo del pendio
# Esempio di profilo di un canale o di una scarpata
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
x: Coordinata orizzontale del punto del pendio.

y: Coordinata verticale del punto del pendio.

r: Coefficiente di restituzione normale per gli impatti in quel punto (valore tra 0 e 1).

et: Coefficiente di restituzione tangenziale per gli impatti in quel punto (valore tra 0 e 1).

Contributi
Contributi, segnalazioni di bug e suggerimenti per miglioramenti sono i benvenuti! Sentiti libero di aprire una "Issue" o inviare una "Pull Request".

Licenza
Questo progetto è distribuito sotto la licenza GNU General Public License v3.0. Vedi il file LICENSE per maggiori dettagli.

