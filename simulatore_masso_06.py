import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget,
    QPushButton, QHBoxLayout, QLabel, QLineEdit, QFormLayout, QFileDialog, QSizePolicy, QMessageBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea
)
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from scipy.interpolate import interp1d
import pandas as pd
import random
import os

# --- Parametri della simulazione (valori predefiniti) ---
G_DEFAULT = 9.81
DT_DEFAULT = 0.01 # Passo temporale per l'integrazione - MODIFICATO PER MIGLIORARE IL RIMBALZO
TEMPO_TOTALE_MAX_DEFAULT = 300 # Tempo massimo di simulazione (limite di sicurezza)

# Parametri per la generazione di massi multipli
NUM_MASSI_DEFAULT = 10
VOLUME_MIN_DEFAULT = 0.1  # m^3
VOLUME_MAX_DEFAULT = 1.0  # m^3
DENSITA_MIN_DEFAULT = 2700 # kg/m^3 (roccia comune)
DENSITA_MAX_DEFAULT = 2800 # kg/m^3

# Parametri per l'intervallo dell'attrito e drag
COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT = 0.03
COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT = 0.07
COEFF_DRAG_ARIA_MIN_DEFAULT = 0.4
COEFF_DRAG_ARIA_MAX_DEFAULT = 0.6 # Correzione qui: COEFF_DRAG_ARIA_MAX_DEFAULT anziché COEFF_DRAG_ARIA_MAX

DENSITA_ARIA_DEFAULT = 1.225
COEFF_RESTITUZIONE_NORMALE_DEFAULT = 0.7 # Coefficiente di restituzione normale di default - MODIFICATO PER MIGLIORARE IL RIMBALZO
COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT = 0.8 # NUOVO: Coefficiente di restituzione tangenziale di default

PENDENZA_GRADI_DEFAULT = 30 # Usato solo se non si carica un CSV
X0_DEFAULT = 0
Y_INIZIALE_OFFSET_DEFAULT = 0.5 # Altezza iniziale del masso sopra il pendio

# Nuovo parametro per il salvataggio dell'output
PASSO_SALVATAGGIO_OUTPUT_DEFAULT = 5.0 # Metri, ogni quanto salvare i dati di energia

# --- Soglie per la terminazione del singolo masso ---
SOGLIA_VELOCITA_ARRESTO = 0.01 # m/s
SOGLIA_DISTANZA_PENDIO_ARRESTO = 0.1 # m (tolleranza per considerare il masso "sul" pendio)
SOGLIA_MARGINE_FINE_PENDIO_X = 0.0 # m (quanto oltre la fine del CSV il masso può andare prima di fermarsi)
SOGLIA_PROFONDITA_FUORI_PENDIO_Y = 50.0 # m (quanto sotto il pendio il masso può andare prima di essere disattivato)

# --- Classe per rappresentare un singolo masso ---
class BloccoMasso:
    def __init__(self, id_blocco, x0, y0, raggio, massa, coeff_attrito_volvimento, coeff_drag_aria, colore=None):
        self.id = id_blocco
        self.raggio = raggio
        self.massa = massa
        self.coeff_attrito_volvimento = coeff_attrito_volvimento
        self.coeff_drag_aria = coeff_drag_aria
        self.colore = colore if colore is not None else (random.random(), random.random(), random.random()) # Colore casuale

        self.posizioni_x = [x0]
        self.posizioni_y = [y0]
        self.velocita_x = [0.0]
        self.velocita_y = [0.0]
        self.velocita_angolare = [0.0] # rad/s
        self.momento_inerzia = (2/5) * self.massa * self.raggio**2 # Per una sfera

        # Dati per il salvataggio in output
        self.output_data = []
        self.last_saved_x = -float('inf') # Per tenere traccia dell'ultima progressiva salvata
        self.last_saved_time = -float('inf') # Per tenere traccia dell'ultimo tempo salvato

        self.is_active = True # Flag per indicare se il masso è ancora in simulazione
        self.final_status = "In simulazione" # Stato finale del masso (es. "Fermato", "Uscito", "Precipitato")
        self.final_time = 0.0 # Tempo in cui il masso si è fermato/uscito
        self.initial_time_step = 0 # Tempo di simulazione al momento della creazione del blocco

# --- Classe per la simulazione fisica ---
class SimulatoreCadutaMassi:
    def __init__(self, parent_gui=None):
        self.parent_gui = parent_gui
        
        # Parametri di default (verranno sovrascritti da set_default_parameters)
        self.G = G_DEFAULT
        self.DT = DT_DEFAULT
        self.TEMPO_TOTALE = TEMPO_TOTALE_MAX_DEFAULT
        self.NUM_MASSI = NUM_MASSI_DEFAULT # Questo verrà gestito correttamente
        self.VOLUME_MIN = VOLUME_MIN_DEFAULT
        self.VOLUME_MAX = VOLUME_MAX_DEFAULT # Correzione qui
        self.DENSITA_MIN = DENSITA_MIN_DEFAULT
        self.DENSITA_MAX = DENSITA_MAX_DEFAULT # Correzione qui
        self.COEFF_ATTRITO_VOLVIMENTO_MIN = COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT
        self.COEFF_ATTRITO_VOLVIMENTO_MAX = COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT # Correzione qui
        self.COEFF_DRAG_ARIA_MIN = COEFF_DRAG_ARIA_MIN_DEFAULT
        self.COEFF_DRAG_ARIA_MAX = COEFF_DRAG_ARIA_MAX_DEFAULT # Correzione qui
        self.DENSITA_ARIA = DENSITA_ARIA_DEFAULT
        self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK = COEFF_RESTITUZIONE_NORMALE_DEFAULT
        self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK = COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT
        self.PENDENZA_GRADI = PENDENZA_GRADI_DEFAULT
        self.PENDENZA_RAD = np.deg2rad(self.PENDENZA_GRADI)
        self.x0_sim = X0_DEFAULT
        self.y_offset_sim = Y_INIZIALE_OFFSET_DEFAULT
        self.PASSO_SALVATAGGIO_OUTPUT = PASSO_SALVATAGGIO_OUTPUT_DEFAULT

        # Dati del pendio e interpolatori
        self.x_pendio_data = None
        self.y_pendio_data = None
        self.r_normale_pendio_data = None
        self.r_tangenziale_pendio_data = None
        self.profilo_pendio_interpolator = None
        self.coeff_restituzione_normale_interpolator = None
        self.coeff_restituzione_tangenziale_interpolator = None
        self.profilo_pendio_func = None 
        self.coeff_restituzione_normale_func = None
        self.coeff_restituzione_tangenziale_func = None 
        self.x_pendio_min_csv = None
        self.x_pendio_max_csv = None

        self.blocchi = [] 
        self.y_riferimento_energia = 0 
        self.tempo_corrente = 0
        self.frame_corrente = 0 

    def set_default_parameters(self):
        """Imposta tutti i parametri della simulazione ai loro valori di default."""
        self.G = G_DEFAULT
        self.DT = DT_DEFAULT
        self.TEMPO_TOTALE = TEMPO_TOTALE_MAX_DEFAULT 
        self.NUM_MASSI = NUM_MASSI_DEFAULT
        self.VOLUME_MIN = VOLUME_MIN_DEFAULT
        self.VOLUME_MAX = VOLUME_MAX_DEFAULT # Correzione qui
        self.DENSITA_MIN = DENSITA_MIN_DEFAULT
        self.DENSITA_MAX = DENSITA_MAX_DEFAULT # Correzione qui
        self.COEFF_ATTRITO_VOLVIMENTO_MIN = COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT
        self.COEFF_ATTRITO_VOLVIMENTO_MAX = COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT # Correzione qui
        self.COEFF_DRAG_ARIA_MIN = COEFF_DRAG_ARIA_MIN_DEFAULT
        self.COEFF_DRAG_ARIA_MAX = COEFF_DRAG_ARIA_MAX_DEFAULT # Correzione qui
        self.DENSITA_ARIA = DENSITA_ARIA_DEFAULT
        self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK = COEFF_RESTITUZIONE_NORMALE_DEFAULT
        self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK = COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT
        self.PENDENZA_GRADI = PENDENZA_GRADI_DEFAULT
        self.PENDENZA_RAD = np.deg2rad(self.PENDENZA_GRADI)
        self.x0_sim = X0_DEFAULT
        self.y_offset_sim = Y_INIZIALE_OFFSET_DEFAULT
        self.PASSO_SALVATAGGIO_OUTPUT = PASSO_SALVATAGGIO_OUTPUT_DEFAULT

        # Resetta anche i dati del pendio per usare quello di default
        self.x_pendio_data = None
        self.y_pendio_data = None
        self.r_normale_pendio_data = None
        self.r_tangenziale_pendio_data = None
        self.profilo_pendio_interpolator = None
        self.coeff_restituzione_normale_interpolator = None
        self.coeff_restituzione_tangenziale_interpolator = None
        self.x_pendio_min_csv = None
        self.x_pendio_max_csv = None
        self.y_riferimento_energia = 0 

    def set_parameters_from_dict(self, params):
        """Aggiorna i parametri della simulazione da un dizionario."""
        self.G = float(params.get('gravita', self.G))
        self.DT = float(params.get('passo_temporale', self.DT))
        self.TEMPO_TOTALE = float(params.get('tempo_totale', self.TEMPO_TOTALE))
        self.NUM_MASSI = int(params.get('num_massi', self.NUM_MASSI)) # Ora aggiorna il valore dell'istanza
        self.VOLUME_MIN = float(params.get('volume_min', self.VOLUME_MIN))
        self.VOLUME_MAX = float(params.get('volume_max', self.VOLUME_MAX))
        self.DENSITA_MIN = float(params.get('densita_min', self.DENSITA_MIN))
        self.DENSITA_MAX = float(params.get('densita_max', self.DENSITA_MAX))
        self.COEFF_ATTRITO_VOLVIMENTO_MIN = float(params.get('attrito_volvimento_min', self.COEFF_ATTRITO_VOLVIMENTO_MIN))
        self.COEFF_ATTRITO_VOLVIMENTO_MAX = float(params.get('attrito_volvimento_max', self.COEFF_ATTRITO_VOLVIMENTO_MAX))
        self.COEFF_DRAG_ARIA_MIN = float(params.get('drag_aria_min', self.COEFF_DRAG_ARIA_MIN))
        self.COEFF_DRAG_ARIA_MAX = float(params.get('drag_aria_max', self.COEFF_DRAG_ARIA_MAX))
        self.DENSITA_ARIA = float(params.get('densita_aria', self.DENSITA_ARIA))
        self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK = float(params.get('coeff_restituzione_normale', self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK))
        self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK = float(params.get('coeff_restituzione_tangenziale', self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK))
        self.PENDENZA_GRADI = float(params.get('pendenza_gradi', self.PENDENZA_GRADI))
        self.PENDENZA_RAD = np.deg2rad(self.PENDENZA_GRADI) # Aggiorna la pendenza in radianti
        self.x0_sim = float(params.get('pos_x_iniziale', self.x0_sim))
        self.y_offset_sim = float(params.get('offset_y_iniziale', self.y_offset_sim))
        self.PASSO_SALVATAGGIO_OUTPUT = float(params.get('passo_salvataggio_output', self.PASSO_SALVATAGGIO_OUTPUT))

    def initialize_simulation_state(self):
        """Prepara la simulazione per un nuovo avvio, utilizzando i parametri correnti."""
        self.blocchi = [] # Cancella i massi precedenti
        self.tempo_corrente = 0
        self.frame_corrente = 0 

        # Se non è stato caricato un CSV o è stato cancellato, usa il profilo rettilineo di default
        if self.profilo_pendio_interpolator is None:
            x_default = np.array([0.0, 100.0])
            y_default = self._crea_profilo_pendio_default_values(x_default) 
            r_normale_default = np.array([self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK, self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK])
            r_tangenziale_default = np.array([self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK, self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK])
            # Chiamo direttamente set_pendio_data per impostare l'interpolatore
            self.set_pendio_data(x_default, y_default, r_normale_default, r_tangenziale_default)
        # else: set_pendio_data è già stato chiamato e i dati sono validi

        # Assicurati che y_riferimento_energia sia impostato correttamente
        if self.y_pendio_data is not None and len(self.y_pendio_data) > 0:
            self.y_riferimento_energia = min(0, np.min(self.y_pendio_data))
        else:
            self.y_riferimento_energia = 0

        self.generare_blocchi() # Genera i massi con il NUM_MASSI corrente

        # Notifica la GUI di ripopolare la tabella del pendio se necessario
        if self.parent_gui and hasattr(self.parent_gui, 'populate_pendio_table'): 
            self.parent_gui.populate_pendio_table()

    def _crea_profilo_pendio_default_values(self, x_values):
        """Genera i valori Y per il profilo pendio di default basato sulla pendenza."""
        return -np.tan(self.PENDENZA_RAD) * x_values

    def set_pendio_data(self, x_data, y_data, r_normale_data, r_tangenziale_data): # Aggiornata la firma
        """
        Aggiorna i dati del pendio e rigenera gli interpolatori.
        Args:
            x_data (np.array): Array delle coordinate X del pendio.
            y_data (np.array): Array delle coordinate Y del pendio.
            r_normale_data (np.array): Array dei coefficienti di restituzione normale.
            r_tangenziale_data (np.array): Array dei coefficienti di restituzione tangenziale. # NUOVO
        """
        if not (len(x_data) == len(y_data) == len(r_normale_data) == len(r_tangenziale_data) and len(x_data) >= 2):
            raise ValueError("I dati del pendio devono avere la stessa lunghezza e almeno 2 punti.")
        
        sorted_indices = np.argsort(x_data)
        self.x_pendio_data = np.array(x_data)[sorted_indices]
        self.y_pendio_data = np.array(y_data)[sorted_indices]
        self.r_normale_pendio_data = np.array(r_normale_data)[sorted_indices]
        self.r_tangenziale_pendio_data = np.array(r_tangenziale_data)[sorted_indices] # NUOVO

        self.profilo_pendio_interpolator = interp1d(self.x_pendio_data, self.y_pendio_data, kind='linear', fill_value="extrapolate")
        self.profilo_pendio_func = self.profilo_pendio_interpolator

        raw_interp_r_normale_func = interp1d(self.x_pendio_data, self.r_normale_pendio_data, kind='linear', fill_value="extrapolate")
        self.coeff_restituzione_normale_interpolator = lambda x: np.clip(raw_interp_r_normale_func(x), 0, 1)
        self.coeff_restituzione_normale_func = self.coeff_restituzione_normale_interpolator 

        # NUOVO: Interpolatore per il coefficiente di restituzione tangenziale
        raw_interp_r_tangenziale_func = interp1d(self.x_pendio_data, self.r_tangenziale_pendio_data, kind='linear', fill_value="extrapolate")
        self.coeff_restituzione_tangenziale_interpolator = lambda x: np.clip(raw_interp_r_tangenziale_func(x), 0, 1)
        self.coeff_restituzione_tangenziale_func = self.coeff_restituzione_tangenziale_interpolator

        self.x_pendio_min_csv = np.min(self.x_pendio_data)
        self.x_pendio_max_csv = np.max(self.x_pendio_data)
        
        self.y_riferimento_energia = min(0, np.min(self.y_pendio_data))

    def carica_dati_csv(self, filepath):
        try:
            df = pd.read_csv(filepath)
            if 'x' not in df.columns or 'y' not in df.columns or 'r' not in df.columns or 'et' not in df.columns:
                raise ValueError("Il file CSV deve contenere le colonne 'x', 'y', 'r' (normale) e 'et' (tangenziale).")

            self.set_pendio_data(df['x'].values, df['y'].values, df['r'].values, df['et'].values)

            print(f"Dati CSV caricati da: {filepath}")
            # Non chiamare initialize_simulation_state o reset_simulazione qui, lo farà la MainWindow dopo il caricamento
            return True
        except Exception as e:
            print(f"Errore durante il caricamento del CSV: {e}")
            self.x_pendio_data = None
            self.y_pendio_data = None
            self.r_normale_pendio_data = None
            self.r_tangenziale_pendio_data = None
            
            self.profilo_pendio_interpolator = None
            self.coeff_restituzione_normale_interpolator = None
            self.coeff_restituzione_tangenziale_interpolator = None
            
            self.x_pendio_min_csv = None
            self.x_pendio_max_csv = None
            self.y_riferimento_energia = 0 

            self.profilo_pendio_func = None
            self.coeff_restituzione_normale_func = None
            self.coeff_restituzione_tangenziale_func = None
            return False

    def generare_blocchi(self):
        self.blocchi = []
        for i in range(self.NUM_MASSI): # Ora usa self.NUM_MASSI, che è stato aggiornato correttamente
            volume = random.uniform(self.VOLUME_MIN, self.VOLUME_MAX)
            densita = random.uniform(self.DENSITA_MIN, self.DENSITA_MAX)
            
            raggio = (volume * 3 / (4 * np.pi))**(1/3)
            massa = volume * densita

            coeff_attrito_volvimento = random.uniform(self.COEFF_ATTRITO_VOLVIMENTO_MIN, self.COEFF_ATTRITO_VOLVIMENTO_MAX)
            coeff_drag_aria = random.uniform(self.COEFF_DRAG_ARIA_MIN, self.COEFF_DRAG_ARIA_MAX)

            x_iniziale = self.x0_sim + random.uniform(-0.5, 0.5)
            y_iniziale = self.profilo_pendio_func(x_iniziale) + raggio + self.y_offset_sim

            blocco = BloccoMasso(i + 1, x_iniziale, y_iniziale, raggio, massa, coeff_attrito_volvimento, coeff_drag_aria)
            # Aggiungi una piccola velocità iniziale in X e Y per simulare una partenza più dinamica
            blocco.velocita_x = [random.uniform(0.1, 1.0)] 
            blocco.velocita_y = [random.uniform(-0.1, -0.5)] 
            blocco.velocita_angolare = [random.uniform(-0.5, 0.5)] # Piccola rotazione iniziale
            self.blocchi.append(blocco)

    def calcola_energie(self, blocco):
        Ek = 0.5 * blocco.massa * (blocco.velocita_x[-1]**2 + blocco.velocita_y[-1]**2)
        Ep = blocco.massa * self.G * (blocco.posizioni_y[-1] - self.y_riferimento_energia)
        Et = Ek + Ep
        return Ek, Ep, Et

    def calcola_un_passo(self):
        self.tempo_corrente += self.DT
        self.frame_corrente += 1
        
        any_mass_active_this_step = False

        for blocco in self.blocchi:
            if not blocco.is_active: 
                continue

            any_mass_active_this_step = True

            x_curr = blocco.posizioni_x[-1] 
            y_curr = blocco.posizioni_y[-1] 
            vx_curr = blocco.velocita_x[-1]
            vy_curr = blocco.velocita_y[-1]
            omega_curr = blocco.velocita_angolare[-1]
            
            # Determina la pendenza e l'angolo
            delta_x_pendenza = 0.1 

            if self.x_pendio_max_csv is not None:
                x_next_for_slope = min(x_curr + delta_x_pendenza, self.x_pendio_max_csv)
            else:
                # Per il pendio di default, estrapola
                x_next_for_slope = x_curr + delta_x_pendenza
            
            delta_x_pendenza_actual = x_next_for_slope - x_curr

            y_curr_pendio_val = self.profilo_pendio_func(x_curr)
            y_next_pendio_val = self.profilo_pendio_func(x_next_for_slope)
            
            if abs(delta_x_pendenza_actual) < 1e-6:
                pendenza_locale = 0.0
            else:
                pendenza_locale = (y_next_pendio_val - y_curr_pendio_val) / delta_x_pendenza_actual
            
            angolo_pendio_rad = np.arctan(pendenza_locale)

            tangente_pendio = np.array([np.cos(angolo_pendio_rad), np.sin(angolo_pendio_rad)])
            normale_pendio = np.array([-np.sin(angolo_pendio_rad), np.cos(angolo_pendio_rad)]) 

            coeff_restituzione_normale_locale = self.coeff_restituzione_normale_func(x_curr)
            coeff_restituzione_tangenziale_locale = self.coeff_restituzione_tangenziale_func(x_curr) 

            # --- Calcolo delle Forze ---
            forza_gravita_x = 0
            forza_gravita_y = -blocco.massa * self.G

            forza_normale_gravita = blocco.massa * self.G * np.cos(angolo_pendio_rad)
            forza_normale_gravita = max(0, forza_normale_gravita)

            forza_attrito_volvimento_magnitudine = blocco.coeff_attrito_volvimento * forza_normale_gravita
            
            velocita_orizzontale_reale = np.sqrt(vx_curr**2 + vy_curr**2)
            if velocita_orizzontale_reale > SOGLIA_VELOCITA_ARRESTO:
                direzione_movimento_x = vx_curr / velocita_orizzontale_reale
                direzione_movimento_y = vy_curr / velocita_orizzontale_reale
            else:
                direzione_movimento_x = 0
                direzione_movimento_y = 0

            forza_attrito_volvimento_x = -forza_attrito_volvimento_magnitudine * direzione_movimento_x
            forza_attrito_volvimento_y = -forza_attrito_volvimento_magnitudine * direzione_movimento_y

            velocita_totale_blocco = np.sqrt(vx_curr**2 + vy_curr**2)
            area_frontale = np.pi * blocco.raggio**2
            forza_drag_magnitudine = 0.5 * self.DENSITA_ARIA * blocco.coeff_drag_aria * area_frontale * velocita_totale_blocco**2
            
            forza_drag_x = -forza_drag_magnitudine * direzione_movimento_x
            forza_drag_y = -forza_drag_magnitudine * direzione_movimento_y

            forza_totale_x = forza_gravita_x + forza_attrito_volvimento_x + forza_drag_x
            forza_totale_y = forza_gravita_y + forza_attrito_volvimento_y + forza_drag_y

            ax = forza_totale_x / blocco.massa
            ay = forza_totale_y / blocco.massa

            vx_new = vx_curr + ax * self.DT
            vy_new = vy_curr + ay * self.DT

            x_new = x_curr + vx_new * self.DT
            y_new = y_curr + vy_new * self.DT

            momento_attrito_volvimento = forza_attrito_volvimento_magnitudine * blocco.raggio
            
            if omega_curr > 0:
                momento_attrito_volvimento = -momento_attrito_volvimento
            elif omega_curr < 0:
                momento_attrito_volvimento = abs(momento_attrito_volvimento)
            else:
                momento_attrito_volvimento = 0.0

            omega_new = omega_curr + (momento_attrito_volvimento / blocco.momento_inerzia) * self.DT


            # --- Gestione collisione con il pendio ---
            y_pendio_sotto_masso = self.profilo_pendio_func(x_new)
            
            if y_new < y_pendio_sotto_masso + blocco.raggio:
                y_new = y_pendio_sotto_masso + blocco.raggio

                vx_vec_pre_impatto = np.array([vx_new, vy_new])

                v_normale_pre_impatto = np.dot(vx_vec_pre_impatto, normale_pendio)
                v_tangenziale_pre_impatto = np.dot(vx_vec_pre_impatto, tangente_pendio)

                coeff_restituzione_normale_locale = self.coeff_restituzione_normale_func(x_new)
                v_normale_post_impatto = -v_normale_pre_impatto * coeff_restituzione_normale_locale

                coeff_restituzione_tangenziale_locale = self.coeff_restituzione_tangenziale_func(x_new) 
                v_tangenziale_post_impatto = v_tangenziale_pre_impatto * coeff_restituzione_tangenziale_locale

                vx_new = v_tangenziale_post_impatto * tangente_pendio[0] + v_normale_post_impatto * normale_pendio[0]
                vy_new = v_tangenziale_post_impatto * tangente_pendio[1] + v_normale_post_impatto * normale_pendio[1]

            # --- Criteri di arresto per il singolo masso ---
            current_speed = np.sqrt(vx_new**2 + vy_new**2)
            if current_speed < SOGLIA_VELOCITA_ARRESTO and \
               abs(y_new - (self.profilo_pendio_func(x_new) + blocco.raggio)) < SOGLIA_DISTANZA_PENDIO_ARRESTO:
                blocco.is_active = False
                blocco.final_status = "Fermato sul Pendio"
                blocco.final_time = self.tempo_corrente
                vx_new = 0.0
                vy_new = 0.0
                omega_new = 0.0
            elif self.x_pendio_max_csv is not None and x_new > self.x_pendio_max_csv + SOGLIA_MARGINE_FINE_PENDIO_X:
                blocco.is_active = False
                blocco.final_status = "Uscito dal Pendio"
                blocco.final_time = self.tempo_corrente
            elif y_new < self.profilo_pendio_func(x_new) - SOGLIA_PROFONDITA_FUORI_PENDIO_Y:
                blocco.is_active = False
                blocco.final_status = "Precipitato"
                blocco.final_time = self.tempo_corrente

            blocco.posizioni_x.append(x_new)
            blocco.posizioni_y.append(y_new)
            blocco.velocita_x.append(vx_new)
            blocco.velocita_y.append(vy_new)
            blocco.velocita_angolare.append(omega_new)

            if (not blocco.output_data) or \
               (x_new >= blocco.last_saved_x + self.PASSO_SALVATAGGIO_OUTPUT) or \
               (not blocco.is_active and blocco.last_saved_x != x_new):
                
                if not blocco.output_data or \
                   (blocco.output_data[-1]['Progressiva_X'] != x_new and blocco.output_data[-1]['Tempo'] != self.tempo_corrente):
                    Ek, Ep, Et = self.calcola_energie(blocco)
                    blocco.output_data.append({
                        'ID_Masso': blocco.id,
                        'Tempo': self.tempo_corrente,
                        'Progressiva_X': x_new,
                        'Altezza_Y': y_new,
                        'Energia_Cinetica': Ek,
                        'Energia_Potenziale': Ep,
                        'Energia_Totale': Et
                    })
                    blocco.last_saved_x = x_new 
                    blocco.last_saved_time = self.tempo_corrente

        if not any_mass_active_this_step and self.tempo_corrente > 0:
            return False
        if self.tempo_corrente >= self.TEMPO_TOTALE:
            for blocco in self.blocchi:
                if blocco.is_active:
                    blocco.is_active = False
                    blocco.final_status = "Tempo Massimo Raggiunto"
                    blocco.final_time = self.tempo_corrente
            return False
        return True

    def salva_output_csv(self, filename):
        all_output_data = []
        for blocco in self.blocchi:
            all_output_data.extend(blocco.output_data)
        
        if not all_output_data:
            print("Nessun dato di simulazione da salvare.")
            return False

        df_output = pd.DataFrame(all_output_data)
        df_output = df_output.sort_values(by=['ID_Masso', 'Tempo']).reset_index(drop=True)
        try:
            df_output.to_csv(filename, index=False)
            print(f"Dati di output salvati in: {filename}")
            return True
        except Exception as e:
            print(f"Errore durante il salvataggio del file CSV: {e}")
            return False


class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.simulatore = SimulatoreCadutaMassi(self)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_animation)
        self.is_simulation_running = False

        self.massi_grafici = []
        self.tracce_grafiche = []

        self.initUI()
        self.simulatore.set_default_parameters() # Inizializza i parametri di default
        self.simulatore.initialize_simulation_state() # Poi inizializza lo stato della simulazione
        self.draw_initial_state()

    def initUI(self):
        self.setWindowTitle('Simulazione Caduta Massi')
        self.setGeometry(100, 100, 1400, 900)

        main_widget = QWidget()
        self.setCentralWidget(main_widget)
        main_layout = QHBoxLayout(main_widget)

        # --- Pannello di Controllo ---
        control_panel_layout = QVBoxLayout()
        control_widget = QWidget()
        control_widget.setLayout(control_panel_layout)
        control_widget.setFixedWidth(350)
        main_layout.addWidget(control_widget)

        title_label = QLabel("Parametri Simulazione")
        title_label.setStyleSheet("font-weight: bold; font-size: 16px; margin-bottom: 5px;")
        control_panel_layout.addWidget(title_label)

        scroll_area = QScrollArea()
        scroll_area.setWidgetResizable(True)
        scroll_content_widget = QWidget()
        self.form_layout = QFormLayout(scroll_content_widget)
        scroll_area.setWidget(scroll_content_widget)
        
        self.param_inputs = {}

        params_data = {
            "Gravità (m/s^2)": (G_DEFAULT, 'gravita'),
            "Passo Temporale (s)": (DT_DEFAULT, 'passo_temporale'),
            "Tempo Max Simulazione (s)": (TEMPO_TOTALE_MAX_DEFAULT, 'tempo_totale'),
            "Densità Aria (kg/m^3)": (DENSITA_ARIA_DEFAULT, 'densita_aria'),
            "Coeff. Restituzione Normale (Default)": (COEFF_RESTITUZIONE_NORMALE_DEFAULT, 'coeff_restituzione_normale'), 
            "Coeff. Restituzione Tangenziale (Default)": (COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT, 'coeff_restituzione_tangenziale'), 
            "Pendenza (gradi) (Default)": (PENDENZA_GRADI_DEFAULT, 'pendenza_gradi'),
            "Pos X Iniziale (m)": (X0_DEFAULT, 'pos_x_iniziale'),
            "Offset Y Iniziale (m)": (Y_INIZIALE_OFFSET_DEFAULT, 'offset_y_iniziale'),
            "-- Parametri Generazione Massi --": None,
            "Numero Massi": (NUM_MASSI_DEFAULT, 'num_massi'),
            "Volume Min (m^3)": (VOLUME_MIN_DEFAULT, 'volume_min'),
            "Volume Max (m^3)": (VOLUME_MAX_DEFAULT, 'volume_max'), # Correzione qui
            "Densità Min (kg/m^3)": (DENSITA_MIN_DEFAULT, 'densita_min'),
            "Densità Max (kg/m^3)": (DENSITA_MAX_DEFAULT, 'densita_max'), # Correzione qui
            "Attrito Rotolamento Min": (COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT, 'attrito_volvimento_min'),
            "Attrito Rotolamento Max": (COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT, 'attrito_volvimento_max'), # Correzione qui
            "Drag Aria Min": (COEFF_DRAG_ARIA_MIN_DEFAULT, 'drag_aria_min'),
            "Drag Aria Max": (COEFF_DRAG_ARIA_MAX_DEFAULT, 'drag_aria_max'), # Correzione qui
            "-- Parametri Output --": None,
            "Passo Salvataggio Output (m)": (PASSO_SALVATAGGIO_OUTPUT_DEFAULT, 'passo_salvataggio_output'),
        }

        for label_text, value_tuple_or_none in params_data.items():
            if value_tuple_or_none is None:
                label = QLabel(label_text)
                label.setStyleSheet("font-weight: bold; margin-top: 10px; margin-bottom: 5px;")
                self.form_layout.addRow(label)
            else:
                default_value, key = value_tuple_or_none
                line_edit = QLineEdit(str(default_value))
                validator = QDoubleValidator()
                validator.setBottom(0.0)
                if key in ['num_massi']:
                    validator = QIntValidator()
                    validator.setBottom(1)
                line_edit.setValidator(validator)
                line_edit.setAccessibleName(label_text)
                self.param_inputs[key] = line_edit
                self.form_layout.addRow(label_text, line_edit) # Aggiunge qui per l'ordine corretto

        control_panel_layout.addWidget(scroll_area)

        button_layout = QHBoxLayout()
        self.start_button = QPushButton("Avvia Simulazione")
        self.start_button.clicked.connect(self.start_simulazione)
        button_layout.addWidget(self.start_button)

        self.pause_button = QPushButton("Ferma Simulazione") 
        self.pause_button.clicked.connect(self.stop_simulation) 
        self.pause_button.setEnabled(False)
        button_layout.addWidget(self.pause_button)

        self.reset_button = QPushButton("Reset Simulazione")
        self.reset_button.clicked.connect(self.reset_simulazione)
        button_layout.addWidget(self.reset_button)
        
        control_panel_layout.addLayout(button_layout)

        self.apply_params_button = QPushButton("Applica Parametri e Riavvia")
        self.apply_params_button.clicked.connect(self.apply_parameters)
        control_panel_layout.addWidget(self.apply_params_button)

        csv_layout = QVBoxLayout()
        csv_label = QLabel("Carica o Modifica Profilo Pendio (CSV)")
        csv_label.setStyleSheet("font-weight: bold; margin-top: 10px; margin-bottom: 5px;")
        csv_layout.addWidget(csv_label)

        self.csv_status_label = QLabel("Nessun CSV caricato, usando pendio di default (30°).")
        csv_layout.addWidget(self.csv_status_label)

        load_csv_button_layout = QHBoxLayout()
        self.load_csv_button = QPushButton("Carica CSV Pendio")
        self.load_csv_button.clicked.connect(self.load_pendio_csv)
        load_csv_button_layout.addWidget(self.load_csv_button)
        
        self.clear_csv_button = QPushButton("Cancella CSV")
        self.clear_csv_button.clicked.connect(self.clear_pendio_csv)
        self.clear_csv_button.setEnabled(False)
        load_csv_button_layout.addWidget(self.clear_csv_button)

        csv_layout.addLayout(load_csv_button_layout)
        
        self.pendio_table_widget = QTableWidget()
        self.pendio_table_widget.setColumnCount(4)
        self.pendio_table_widget.setHorizontalHeaderLabels(['X (m)', 'Y (m)', 'Rn', 'Et'])
        self.pendio_table_widget.horizontalHeader().setSectionResizeMode(QHeaderView.Stretch)
        self.pendio_table_widget.setMinimumHeight(150)
        csv_layout.addWidget(self.pendio_table_widget)

        table_buttons_layout = QHBoxLayout()
        self.add_row_button = QPushButton("Aggiungi Riga")
        self.add_row_button.clicked.connect(self.add_pendio_row)
        table_buttons_layout.addWidget(self.add_row_button)

        self.remove_row_button = QPushButton("Rimuovi Riga Selezionata")
        self.remove_row_button.clicked.connect(self.remove_pendio_row)
        table_buttons_layout.addWidget(self.remove_row_button)

        self.apply_pendio_button = QPushButton("Applica Modifiche Pendio")
        self.apply_pendio_button.clicked.connect(self.apply_pendio_changes)
        table_buttons_layout.addWidget(self.apply_pendio_button)

        csv_layout.addLayout(table_buttons_layout)
        control_panel_layout.addLayout(csv_layout)

        output_button_layout = QHBoxLayout()
        self.show_summary_plot_button = QPushButton("Mostra Riepilogo Finale")
        self.show_summary_plot_button.clicked.connect(self.show_summary_plot)
        self.show_summary_plot_button.setEnabled(False)
        output_button_layout.addWidget(self.show_summary_plot_button)

        self.save_output_csv_button = QPushButton("Salva Dati Output CSV")
        self.save_output_csv_button.clicked.connect(self.save_output_csv_file)
        self.save_output_csv_button.setEnabled(False)
        output_button_layout.addWidget(self.save_output_csv_button)

        control_panel_layout.addLayout(output_button_layout)
        control_panel_layout.addStretch(1)

        self.figure, self.ax = plt.subplots(figsize=(10, 8))
        self.canvas = FigureCanvas(self.figure)
        main_layout.addWidget(self.canvas)

        self.ax.set_aspect('equal', adjustable='box')
        self.ax.set_title('Simulazione Caduta Massi')
        self.ax.set_xlabel('Posizione X (m)')
        self.ax.set_ylabel('Posizione Y (m)')
        self.ax.grid(True)
        self.linea_pendio, = self.ax.plot([], [], 'k-', label='Profilo Pendio', zorder=1)
        self.ax.legend()
        
        self.populate_pendio_table()
        self.draw_initial_state()


    def create_float_validator(self):
        validator = QDoubleValidator()
        validator.setBottom(0.0)
        return validator

    def draw_initial_state(self):
        for masso_g in self.massi_grafici:
            masso_g.remove()
        for traccia_g in self.tracce_grafiche:
            traccia_g.remove()
        self.massi_grafici = []
        self.tracce_grafiche = []

        x_plot_range = []
        if self.simulatore.x_pendio_data is not None and len(self.simulatore.x_pendio_data) >= 2:
            x_plot_range = np.linspace(min(self.simulatore.x_pendio_data), max(self.simulatore.x_pendio_data), 200)
        else: # Usa il pendio di default basato sulla pendenza corrente
            # Assicurati che il simulatore abbia il profilo pendio di default impostato prima di disegnare
            if self.simulatore.profilo_pendio_func is None:
                x_default = np.array([0.0, 100.0])
                self.simulatore.set_pendio_data(x_default, self.simulatore._crea_profilo_pendio_default_values(x_default),
                                                np.array([self.simulatore.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK, self.simulatore.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK]),
                                                np.array([self.simulatore.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK, self.simulatore.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK]))
            x_plot_range = np.linspace(self.simulatore.x0_sim - 10, self.simulatore.x0_sim + 100, 200)
            
        y_plot_range = self.simulatore.profilo_pendio_func(x_plot_range)
        self.linea_pendio.set_data(x_plot_range, y_plot_range)

        for i, blocco in enumerate(self.simulatore.blocchi):
            masso_g, = self.ax.plot([blocco.posizioni_x[0]], [blocco.posizioni_y[0]], 'o', markersize=blocco.raggio * 20, color=blocco.colore, zorder=5)
            self.massi_grafici.append(masso_g)
            traccia_g, = self.ax.plot([], [], '--', linewidth=1, alpha=0.7, color=blocco.colore, zorder=2)
            self.tracce_grafiche.append(traccia_g)

        all_x_data = list(x_plot_range)
        all_y_data = list(y_plot_range)
        for blocco in self.simulatore.blocchi:
            all_x_data.append(blocco.posizioni_x[0])
            all_y_data.append(blocco.posizioni_y[0])

        if all_x_data and all_y_data:
            min_x, max_x = np.min(all_x_data), np.max(all_x_data)
            min_y, max_y = np.min(all_y_data), np.max(all_y_data)
            
            margin_x = (max_x - min_x) * 0.1
            margin_y = (max_y - min_y) * 0.1
            
            self.ax.set_xlim(min_x - margin_x, max_x + margin_x)
            self.ax.set_ylim(min_y - margin_y - 2, max_y + margin_y + 2)

        self.canvas.draw()


    def start_simulazione(self):
        if not self.is_simulation_running:
            self.timer.start(int(self.simulatore.DT * 1000))
            self.is_simulation_running = True
            self.start_button.setEnabled(False)
            self.pause_button.setEnabled(True)
            self.reset_button.setEnabled(False)
            self.load_csv_button.setEnabled(False)
            self.clear_csv_button.setEnabled(False)
            self.pendio_table_widget.setEnabled(False)
            self.add_row_button.setEnabled(False)
            self.remove_row_button.setEnabled(False)
            self.apply_pendio_button.setEnabled(False)
            self.apply_params_button.setEnabled(False)
            for key in self.param_inputs:
                self.param_inputs[key].setEnabled(False)


    def stop_simulation(self):
        if self.is_simulation_running:
            self.timer.stop()
            self.is_simulation_running = False
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.load_csv_button.setEnabled(True)
            if self.simulatore.x_pendio_data is not None: # Se c'è un CSV caricato, riabilita clear
                 self.clear_csv_button.setEnabled(True)
            self.pendio_table_widget.setEnabled(True)
            self.add_row_button.setEnabled(True)
            self.remove_row_button.setEnabled(True)
            self.apply_pendio_button.setEnabled(True)
            self.apply_params_button.setEnabled(True)
            for key in self.param_inputs:
                self.param_inputs[key].setEnabled(True)
            self.on_simulation_finished()

    def reset_simulazione(self):
        self.timer.stop()
        self.is_simulation_running = False
        self.simulatore.set_default_parameters() # Resetta tutti i parametri ai default
        self.simulatore.initialize_simulation_state() # Poi inizializza lo stato con questi default
        self.draw_initial_state()
        self.populate_pendio_table() # Aggiorna la tabella con il pendio di default
        # Imposta lo stato iniziale dei pulsanti e input
        self.start_button.setEnabled(True)
        self.pause_button.setEnabled(False)
        self.reset_button.setEnabled(True)
        self.load_csv_button.setEnabled(True)
        self.clear_csv_button.setEnabled(False) # Dopo un reset, nessun CSV caricato
        self.show_summary_plot_button.setEnabled(False)
        self.save_output_csv_button.setEnabled(False)
        self.pendio_table_widget.setEnabled(True)
        self.add_row_button.setEnabled(True)
        self.remove_row_button.setEnabled(True)
        self.apply_pendio_button.setEnabled(True)
        self.apply_params_button.setEnabled(True)
        for key in self.param_inputs:
            self.param_inputs[key].setEnabled(True)
        self.param_inputs['pendenza_gradi'].setEnabled(True) # Abilita pendenza gradi dopo reset

    def on_simulation_finished(self):
        self.stop_simulation()
        self.show_summary_plot_button.setEnabled(True)
        self.save_output_csv_button.setEnabled(True)

    def apply_parameters(self):
        params = {}
        try:
            for key, line_edit in self.param_inputs.items():
                if line_edit.text():
                    if key == 'num_massi':
                        params[key] = int(line_edit.text())
                        if params[key] <= 0:
                            raise ValueError(f"'{line_edit.accessibleName()}' deve essere maggiore di zero.")
                    else:
                        params[key] = float(line_edit.text())
                        if params[key] < 0:
                            raise ValueError(f"'{line_edit.accessibleName()}' non può essere negativo.")
            
            if params.get('passo_salvataggio_output') is not None and params.get('passo_salvataggio_output') <= 0:
                raise ValueError("Il 'Passo Salvataggio Output' deve essere maggiore di zero.")

            self.simulatore.set_parameters_from_dict(params) # Imposta i nuovi parametri
            self.simulatore.initialize_simulation_state() # Inizializza lo stato con i nuovi parametri
            self.draw_initial_state() # Disegna lo stato iniziale
            QMessageBox.information(self, "Parametri Applicati", "I nuovi parametri sono stati applicati. La simulazione è stata resettata con i nuovi valori.")
            self.start_button.setEnabled(True)
            self.pause_button.setEnabled(False)
            self.reset_button.setEnabled(True)
            self.show_summary_plot_button.setEnabled(False)
            self.save_output_csv_button.setEnabled(False)

        except ValueError as e:
            QMessageBox.warning(self, "Errore Parametri", f"Errore nei parametri inseriti: {e}. Assicurati che siano valori numerici validi e positivi dove richiesto.")


    def update_animation(self):
        sim_continuing = self.simulatore.calcola_un_passo()
        
        if sim_continuing:
            min_x_visible, max_x_visible = float('inf'), float('-inf')
            min_y_visible, max_y_visible = float('inf'), float('-inf')

            for i, blocco in enumerate(self.simulatore.blocchi):
                if blocco.is_active:
                    self.massi_grafici[i].set_data([blocco.posizioni_x[-1]], [blocco.posizioni_y[-1]])
                    self.tracce_grafiche[i].set_data(blocco.posizioni_x, blocco.posizioni_y)
                    
                    min_x_visible = min(min_x_visible, blocco.posizioni_x[-1])
                    max_x_visible = max(max_x_visible, blocco.posizioni_x[-1])
                    min_y_visible = min(min_y_visible, blocco.posizioni_y[-1])
                    max_y_visible = max(max_y_visible, blocco.posizioni_y[-1])
                else:
                    if len(blocco.posizioni_x) > 0:
                        min_x_visible = min(min_x_visible, blocco.posizioni_x[-1])
                        max_x_visible = max(max_x_visible, blocco.posizioni_x[-1])
                        min_y_visible = min(min_y_visible, blocco.posizioni_y[-1])
                        max_y_visible = max(max_y_visible, blocco.posizioni_y[-1])

            if max_x_visible > -float('inf'):
                current_xlim = self.ax.get_xlim()
                current_ylim = self.ax.get_ylim()

                if max_x_visible + 10 > current_xlim[1]:
                    self.ax.set_xlim(current_xlim[0], max_x_visible + 50)
                if min_x_visible - 10 < current_xlim[0]:
                    self.ax.set_xlim(min_x_visible - 50, current_xlim[1])
                
                x_min_bound = self.simulatore.x_pendio_data.min() if self.simulatore.x_pendio_data is not None else 0
                x_max_bound = self.simulatore.x_pendio_data.max() if self.simulatore.x_pendio_data is not None else 100

                pendio_y_at_visible_min_x = self.simulatore.profilo_pendio_func(max(min_x_visible - 20, x_min_bound))
                pendio_y_at_visible_max_x = self.simulatore.profilo_pendio_func(min(max_x_visible + 20, x_max_bound))

                overall_min_y = min(min_y_visible, pendio_y_at_visible_min_x, pendio_y_at_visible_max_x)
                overall_max_y = max(max_y_visible, pendio_y_at_visible_min_x, pendio_y_at_visible_max_x)

                if overall_min_y - 10 < current_ylim[0] or overall_max_y + 10 > current_ylim[1]:
                    self.ax.set_ylim(overall_min_y - 20, overall_max_y + 20)

            self.canvas.draw()
            self.canvas.flush_events()
        else:
            self.on_simulation_finished()

    def load_pendio_csv(self):
        options = QFileDialog.Options()
        fileName, _ = QFileDialog.getOpenFileName(self, "Carica File CSV del Pendio", "", "CSV Files (*.csv);;All Files (*)", options=options)
        if fileName:
            if self.simulatore.carica_dati_csv(fileName):
                self.csv_status_label.setText(f"CSV caricato: {os.path.basename(fileName)}")
                self.clear_csv_button.setEnabled(True)
                self.param_inputs['pendenza_gradi'].setEnabled(False)
                self.simulatore.initialize_simulation_state() # Inizializza lo stato dopo il caricamento del CSV
            else:
                self.csv_status_label.setText("Errore caricamento CSV.")
                self.clear_csv_button.setEnabled(False)
                self.param_inputs['pendenza_gradi'].setEnabled(True)
            self.populate_pendio_table()
            self.draw_initial_state()

    def clear_pendio_csv(self):
        # Questo farà sì che initialize_simulation_state ricrei il pendio di default
        self.simulatore.x_pendio_data = None
        self.simulatore.y_pendio_data = None
        self.simulatore.r_normale_pendio_data = None
        self.simulatore.r_tangenziale_pendio_data = None
        
        self.simulatore.profilo_pendio_interpolator = None
        self.simulatore.coeff_restituzione_normale_interpolator = None
        self.simulatore.coeff_restituzione_tangenziale_interpolator = None
        
        self.simulatore.x_pendio_min_csv = None
        self.simulatore.x_pendio_max_csv = None
        self.simulatore.y_riferimento_energia = 0 

        self.simulatore.profilo_pendio_func = None
        self.simulatore.coeff_restituzione_normale_func = None
        self.simulatore.coeff_restituzione_tangenziale_func = None
        
        # Ora inizializza lo stato. initialize_simulation_state rileverà che non c'è un CSV e userà i valori di default
        self.simulatore.initialize_simulation_state() 
        self.csv_status_label.setText("Nessun CSV caricato, usando pendio di default (30°).")
        self.clear_csv_button.setEnabled(False)
        self.param_inputs['pendenza_gradi'].setEnabled(True)
        self.populate_pendio_table()
        self.draw_initial_state()
        QMessageBox.information(self, "CSV Cancellato", "I dati del pendio dal CSV sono stati rimossi. Verrà utilizzato il profilo di default.")


    def populate_pendio_table(self):
        self.pendio_table_widget.setRowCount(0)
        # Assicurati che self.simulatore.x_pendio_data esista e sia popolato, se no non fa nulla
        if self.simulatore.x_pendio_data is not None and len(self.simulatore.x_pendio_data) > 0:
            for i in range(len(self.simulatore.x_pendio_data)):
                row_position = self.pendio_table_widget.rowCount()
                self.pendio_table_widget.insertRow(row_position)
                self.pendio_table_widget.setItem(row_position, 0, QTableWidgetItem(str(self.simulatore.x_pendio_data[i])))
                self.pendio_table_widget.setItem(row_position, 1, QTableWidgetItem(str(self.simulatore.y_pendio_data[i])))
                self.pendio_table_widget.setItem(row_position, 2, QTableWidgetItem(str(self.simulatore.r_normale_pendio_data[i])))
                self.pendio_table_widget.setItem(row_position, 3, QTableWidgetItem(str(self.simulatore.r_tangenziale_pendio_data[i])))

    def add_pendio_row(self):
        row_position = self.pendio_table_widget.rowCount()
        self.pendio_table_widget.insertRow(row_position)
        self.pendio_table_widget.setItem(row_position, 0, QTableWidgetItem(""))
        self.pendio_table_widget.setItem(row_position, 1, QTableWidgetItem(""))
        self.pendio_table_widget.setItem(row_position, 2, QTableWidgetItem(str(COEFF_RESTITUZIONE_NORMALE_DEFAULT)))
        self.pendio_table_widget.setItem(row_position, 3, QTableWidgetItem(str(COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT)))

    def remove_pendio_row(self):
        selected_rows = self.pendio_table_widget.selectionModel().selectedRows()
        if selected_rows:
            for row in sorted(selected_rows, reverse=True):
                self.pendio_table_widget.removeRow(row.row())
        else:
            QMessageBox.information(self, "Nessuna Riga Selezionata", "Seleziona una riga da rimuovere.")

    def apply_pendio_changes(self):
        x_data = []
        y_data = []
        r_normale_data = []
        r_tangenziale_data = []
        try:
            if self.pendio_table_widget.rowCount() < 2:
                raise ValueError("Il pendio deve avere almeno 2 punti definiti nella tabella.")

            for row in range(self.pendio_table_widget.rowCount()):
                x_item = self.pendio_table_widget.item(row, 0)
                y_item = self.pendio_table_widget.item(row, 1)
                r_n_item = self.pendio_table_widget.item(row, 2)
                r_t_item = self.pendio_table_widget.item(row, 3)

                if x_item is None or y_item is None or r_n_item is None or r_t_item is None or \
                   not x_item.text() or not y_item.text() or not r_n_item.text() or not r_t_item.text():
                    raise ValueError(f"Riga {row+1} incompleta. Assicurati che tutte le celle siano popolate.")

                x = float(x_item.text())
                y = float(y_item.text())
                r_n = float(r_n_item.text())
                r_t = float(r_t_item.text())

                if not (0.0 <= r_n <= 1.0):
                    raise ValueError(f"Il coefficiente Rn nella riga {row+1} ({r_n}) deve essere tra 0 e 1.")
                if not (0.0 <= r_t <= 1.0):
                    raise ValueError(f"Il coefficiente Et nella riga {row+1} ({r_t}) deve essere tra 0 e 1.")

                x_data.append(x)
                y_data.append(y)
                r_normale_data.append(r_n)
                r_tangenziale_data.append(r_t)

            for i in range(1, len(x_data)):
                if x_data[i] <= x_data[i-1]:
                    raise ValueError(f"I valori X devono essere strettamente crescenti. Problema alla riga {i+1} (X={x_data[i]} <= X={x_data[i-1]}).")

            self.simulatore.set_pendio_data(np.array(x_data), np.array(y_data), np.array(r_normale_data), np.array(r_tangenziale_data))
            self.simulatore.initialize_simulation_state() # Inizializza lo stato con i nuovi dati del pendio
            self.draw_initial_state()
            self.csv_status_label.setText("Pendio modificato manualmente.")
            self.clear_csv_button.setEnabled(True)
            self.param_inputs['pendenza_gradi'].setEnabled(False)
            QMessageBox.information(self, "Modifiche Pendio Applicate", "Le modifiche al profilo del pendio sono state applicate con successo. La simulazione verrà resettata.")
        except ValueError as e:
            QMessageBox.warning(self, "Errore Modifica Pendio", f"Errore durante l'applicazione delle modifiche: {e}")
        except Exception as e:
            QMessageBox.critical(self, "Errore Inatteso", f"Si è verificato un errore inatteso: {e}")

    def show_summary_plot(self):
        plt.close('all') # Chiude tutte le finestre di plot precedenti

        summary_table_data = []
        stopped_on_slope_count = 0
        reached_end_count = 0
        precipitated_count = 0
        total_masses = len(self.simulatore.blocchi)

        for blocco in self.simulatore.blocchi:
            final_x = blocco.posizioni_x[-1] if blocco.posizioni_x else np.nan
            final_y = blocco.posizioni_y[-1] if blocco.posizioni_y else np.nan
            final_ek = blocco.output_data[-1]['Energia_Cinetica'] / 1000 if blocco.output_data else 0.0
            final_time = blocco.final_time if blocco.final_time > 0 else self.simulatore.tempo_corrente

            # Calcola l'altezza dal suolo
            # Se il masso è fermo, l'altezza dal suolo dovrebbe essere raggio
            altezza_dal_suolo = final_y - self.simulatore.profilo_pendio_func(final_x)

            status = blocco.final_status
            if status == "Fermato sul Pendio":
                stopped_on_slope_count += 1
                # Se è fermo sul pendio, la sua "altezza dal suolo" è il suo raggio
                # Usiamo il raggio del masso per indicare che è a contatto con il pendio
                altezza_dal_suolo = blocco.raggio 
            elif status == "Uscito dal Pendio":
                reached_end_count += 1
                # Per i massi usciti dal pendio, l'altezza dal suolo è la loro altezza attuale - altezza del pendio
                # Non aggiungiamo il raggio, perché non sono "rotolanti" sul pendio ma "caduti" dalla fine
                pass # altezza_dal_suolo è già calcolata correttamente
            elif status == "Tempo Massimo Raggiunto":
                # Se il masso è attivo fino alla fine della simulazione, consideriamolo come "arrivato in fondo"
                # a meno che non sia rimasto molto indietro rispetto alla fine del pendio CSV
                if self.simulatore.x_pendio_max_csv is not None and final_x < self.simulatore.x_pendio_max_csv - SOGLIA_MARGINE_FINE_PENDIO_X:
                     status = "Fermato per Tempo Massimo (non a fine pendio)"
                     stopped_on_slope_count += 1 # Conta come fermato se non è arrivato alla fine
                     altezza_dal_suolo = blocco.raggio # Assume che sia atterrato sul pendio
                else:
                    reached_end_count += 1
                    # altezza_dal_suolo già calcolata
            elif status == "Precipitato":
                precipitated_count += 1
                # altezza_dal_suolo già calcolata, sarà un valore negativo significativo
            elif status == "In simulazione":
                 status = "N/D (Simulazione incompleta)"
                 altezza_dal_suolo = "N/A" # Non applicabile se la simulazione è incompleta per questo masso


            summary_table_data.append([
                blocco.id,
                f"{final_x:.2f}",
                f"{altezza_dal_suolo:.2f}" if isinstance(altezza_dal_suolo, (float, int)) else altezza_dal_suolo,
                f"{final_ek:.2f}",
                status,
                f"{final_time:.2f}"
            ])
        
        if not summary_table_data:
            QMessageBox.information(self, "Nessun Dato", "Nessun dato di simulazione disponibile per il riepilogo.")
            return

        # --- Tabella Dettagliata ---
        fig_table, ax_table = plt.subplots(figsize=(12, min(8, 0.5 * total_masses + 2)))
        ax_table.axis('off')
        ax_table.set_title('Riepilogo Dettagliato Massi', y=1.05) 

        col_labels = ["ID Masso", "X Finale (m)", "Altezza dal Suolo (m)", "Ek Finale (kJ)", "Stato Finale", "Tempo Finale (s)"] # Aggiornata l'etichetta

        table = ax_table.table(
            cellText=summary_table_data,
            colLabels=col_labels,
            loc='center', 
            cellLoc='center' 
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10) 
        table.scale(1.0, 1.2)

        # --- Riepilogo Percentuale (testuale e grafico a torta) ---
        fig_summary, (ax_text, ax_pie) = plt.subplots(1, 2, figsize=(14, 7)) # Due subplot per testo e torta

        # Testo di riepilogo
        ax_text.axis('off')
        ax_text.set_title('Riepilogo Percentuale', y=0.8)

        # Data for the pie chart
        labels = []
        sizes = []
        
        if reached_end_count > 0:
            labels.append(f'Arrivati a fine pendio ({reached_end_count})')
            sizes.append(reached_end_count)
        if stopped_on_slope_count > 0:
            labels.append(f'Fermati sul pendio ({stopped_on_slope_count})')
            sizes.append(stopped_on_slope_count)
        if precipitated_count > 0:
            labels.append(f'Precipitati ({precipitated_count})')
            sizes.append(precipitated_count)

        if not sizes: # Se non ci sono massi o nessuno ha raggiunto uno stato finale chiaro
            labels = ['Nessun dato']
            sizes = [1] # Per disegnare una torta vuota
            colors = ['lightgray']
        else:
            # Colori predefiniti o generati casualmente
            colors = plt.cm.Paired(np.arange(len(labels))) # Usa una colormap per colori distinti

        ax_pie.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors)
        ax_pie.axis('equal')  # Assicura che il grafico a torta sia circolare.
        ax_pie.set_title('Distribuzione Stati Finali Massi')

        plt.tight_layout()
        plt.show()

    def save_output_csv_file(self):
        if not any(blocco.output_data for blocco in self.simulatore.blocchi):
            QMessageBox.information(self, "Nessun Dato", "Nessun dato di simulazione da salvare.")
            return

        options = QFileDialog.Options()
        default_filename = f"simulazione_massi_output_{pd.Timestamp.now().strftime('%Y%m%d_%H%M%S')}.csv"
        fileName, _ = QFileDialog.getSaveFileName(self, "Salva Dati Output CSV", default_filename, "CSV Files (*.csv);;All Files (*)", options=options)
        if fileName:
            if self.simulatore.salva_output_csv(fileName):
                QMessageBox.information(self, "Salvataggio Completato", f"Dati di output salvati con successo in:\n{fileName}")
            else:
                QMessageBox.warning(self, "Errore Salvataggio", f"Si è verificato un errore durante il salvataggio dei dati.")

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
