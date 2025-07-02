import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QVBoxLayout, QWidget,
    QPushButton, QHBoxLayout, QLabel, QLineEdit, QFormLayout, QFileDialog, QMessageBox,
    QTableWidget, QTableWidgetItem, QHeaderView, QScrollArea, QComboBox
)
from PyQt5.QtCore import QTimer, Qt
from PyQt5.QtGui import QDoubleValidator, QIntValidator
from scipy.interpolate import interp1d
import pandas as pd
import random
import os

# --- Parametri di Simulazione (Valori Predefiniti) ---
G_DEFAULT = 9.81
DT_DEFAULT = 0.01
TEMPO_TOTALE_MAX_DEFAULT = 300

# Parametri per la generazione dei massi
NUM_MASSI_DEFAULT = 10
DENSITA_DEFAULT = 2700  # kg/m^3

# --- NUOVI PARAMETRI PER LA FORMA ---
FORMA_BLOCCO_DEFAULT = 'Prisma'  # 'Sfera' o 'Prisma'
# Dimensioni per Sfera (Volume in m^3)
VOLUME_MIN_DEFAULT = 0.1
VOLUME_MAX_DEFAULT = 1.0
# Dimensioni per Prisma (in metri)
LARGHEZZA_MIN_DEFAULT = 0.5
LARGHEZZA_MAX_DEFAULT = 1.5
ALTEZZA_MIN_DEFAULT = 0.3
ALTEZZA_MAX_DEFAULT = 1.0
PROFONDITA_MIN_DEFAULT = 0.4
PROFONDITA_MAX_DEFAULT = 1.2

# Parametri fisici
COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT = 0.03
COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT = 0.07
COEFF_DRAG_ARIA_MIN_DEFAULT = 0.8  # Valore tipico per un prisma
COEFF_DRAG_ARIA_MAX_DEFAULT = 1.2
DENSITA_ARIA_DEFAULT = 1.225
COEFF_RESTITUZIONE_NORMALE_DEFAULT = 0.7
COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT = 0.8

# Parametri del pendio di default
PENDENZA_GRADI_DEFAULT = 30
X0_DEFAULT = 0
Y_INIZIALE_OFFSET_DEFAULT = 0.5

PASSO_SALVATAGGIO_OUTPUT_DEFAULT = 5.0

# --- Soglie di Arresto ---
SOGLIA_VELOCITA_ARRESTO = 0.01
SOGLIA_DISTANZA_PENDIO_ARRESTO = 0.1
SOGLIA_MARGINE_FINE_PENDIO_X = 0.0
SOGLIA_PROFONDITA_FUORI_PENDIO_Y = 50.0

# --- Classe per rappresentare un singolo blocco (MODIFICATA) ---
class BloccoMasso:
    def __init__(self, id_blocco, x0, y0, densita, forma, dimensioni, coeff_attrito_volvimento, coeff_drag_aria, colore=None):
        self.id = id_blocco
        self.densita = densita
        self.forma = forma
        self.dimensioni = dimensioni
        self.coeff_attrito_volvimento = coeff_attrito_volvimento
        self.coeff_drag_aria = coeff_drag_aria
        self.colore = colore if colore is not None else (random.random(), random.random(), random.random())

        self.calcola_proprieta_geometriche()

        self.posizioni_x = [x0]
        self.posizioni_y = [y0]
        self.velocita_x = [0.0]
        self.velocita_y = [0.0]
        self.angolo = [0.0]
        self.velocita_angolare = [0.0]

        self.output_data = []
        self.last_saved_x = -float('inf')
        self.last_saved_time = -float('inf')

        self.is_active = True
        self.final_status = "In simulazione"
        self.final_time = 0.0

    def calcola_proprieta_geometriche(self):
        if self.forma == 'Sfera':
            self.raggio = self.dimensioni[0]
            self.volume = (4 / 3) * np.pi * self.raggio**3
            self.massa = self.volume * self.densita
            self.momento_inerzia = (2 / 5) * self.massa * self.raggio**2
            self.area_frontale = np.pi * self.raggio**2
            self.dist_centro_vertice = self.raggio
        elif self.forma == 'Prisma':
            l, h, p = self.dimensioni
            self.volume = l * h * p
            self.massa = self.volume * self.densita
            # Momento d'inerzia per un prisma rispetto a un asse passante per il centro di massa e parallelo alla profondità
            self.momento_inerzia = (1 / 12) * self.massa * (l**2 + h**2)
            self.area_frontale = l * h # Assumiamo area proiettata come larghezza * altezza
            self.dist_centro_vertice = np.sqrt((l / 2)**2 + (h / 2)**2)
        else:
            raise ValueError(f"Forma '{self.forma}' non supportata.")

# --- Classe per la simulazione fisica ---
class SimulatoreCadutaMassi:
    def __init__(self, parent_gui=None):
        self.parent_gui = parent_gui
        self.x_pendio_data = None
        self.y_pendio_data = None
        self.profilo_pendio_interpolator = None
        self.blocchi = []
        
        self.set_default_parameters() # Imposta i parametri e crea il pendio di default inizialmente
        self.reset_state() # Resetta solo lo stato della simulazione (tempo, frame, ecc.)


    def reset_state(self):
        self.tempo_corrente = 0
        self.frame_corrente = 0
        
        # Aggiorna y_riferimento_energia solo se ci sono dati validi del pendio
        if self.x_pendio_data is not None and len(self.y_pendio_data) > 0:
            self.y_riferimento_energia = min(0, np.min(self.y_pendio_data))
        else:
            self.y_riferimento_energia = 0 # Valore di fallback
        
        # NON chiamare _crea_profilo_pendio_default() qui.
        # Il pendio dovrebbe essere già impostato da __init__, set_default_parameters o carica_dati_csv.


    def set_default_parameters(self):
        self.G = G_DEFAULT
        self.DT = DT_DEFAULT
        self.TEMPO_TOTALE = TEMPO_TOTALE_MAX_DEFAULT
        self.NUM_MASSI = NUM_MASSI_DEFAULT
        self.DENSITA = DENSITA_DEFAULT
        self.FORMA_BLOCCO = FORMA_BLOCCO_DEFAULT
        self.VOLUME_MIN = VOLUME_MIN_DEFAULT
        self.VOLUME_MAX = VOLUME_MAX_DEFAULT
        self.LARGHEZZA_MIN = LARGHEZZA_MIN_DEFAULT
        self.LARGHEZZA_MAX = LARGHEZZA_MAX_DEFAULT
        self.ALTEZZA_MIN = ALTEZZA_MIN_DEFAULT
        self.ALTEZZA_MAX = ALTEZZA_MAX_DEFAULT
        self.PROFONDITA_MIN = PROFONDITA_MIN_DEFAULT
        self.PROFONDITA_MAX = PROFONDITA_MAX_DEFAULT
        self.COEFF_ATTRITO_VOLVIMENTO_MIN = COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT
        self.COEFF_ATTRITO_VOLVIMENTO_MAX = COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT
        self.COEFF_DRAG_ARIA_MIN = COEFF_DRAG_ARIA_MIN_DEFAULT
        self.COEFF_DRAG_ARIA_MAX = COEFF_DRAG_ARIA_MAX_DEFAULT
        self.DENSITA_ARIA = DENSITA_ARIA_DEFAULT
        self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK = COEFF_RESTITUZIONE_NORMALE_DEFAULT
        self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK = COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT
        self.PENDENZA_GRADI = PENDENZA_GRADI_DEFAULT
        self.PENDENZA_RAD = np.deg2rad(self.PENDENZA_GRADI)
        self.x0_sim = X0_DEFAULT
        self.y_offset_sim = Y_INIZIALE_OFFSET_DEFAULT
        self.PASSO_SALVATAGGIO_OUTPUT = PASSO_SALVATAGGIO_OUTPUT_DEFAULT
        self.x_pendio_min_csv = 0 # Valori di fallback
        self.x_pendio_max_csv = 200 # Valori di fallback
        
        self._crea_profilo_pendio_default() # Assicurati che il pendio di default sia sempre impostato


    def set_parameters_from_dict(self, params):
        self.G = float(params.get('gravita', self.G))
        self.DT = float(params.get('passo_temporale', self.DT))
        self.TEMPO_TOTALE = float(params.get('tempo_totale', self.TEMPO_TOTALE))
        self.NUM_MASSI = int(params.get('num_massi', self.NUM_MASSI))
        self.DENSITA = float(params.get('densita', self.DENSITA))
        self.FORMA_BLOCCO = params.get('forma_blocco', self.FORMA_BLOCCO)
        self.VOLUME_MIN = float(params.get('volume_min', self.VOLUME_MIN))
        self.VOLUME_MAX = float(params.get('volume_max', self.VOLUME_MAX))
        self.LARGHEZZA_MIN = float(params.get('larghezza_min', self.LARGHEZZA_MIN))
        self.LARGHEZZA_MAX = float(params.get('larghezza_max', self.LARGHEZZA_MAX))
        self.ALTEZZA_MIN = float(params.get('altezza_min', self.ALTEZZA_MIN))
        self.ALTEZZA_MAX = float(params.get('altezza_max', self.ALTEZZA_MAX))
        self.PROFONDITA_MIN = float(params.get('profondita_min', self.PROFONDITA_MIN))
        self.PROFONDITA_MAX = float(params.get('profondita_max', self.PROFONDITA_MAX))
        self.COEFF_ATTRITO_VOLVIMENTO_MIN = float(params.get('attrito_volvimento_min', self.COEFF_ATTRITO_VOLVIMENTO_MIN))
        self.COEFF_ATTRITO_VOLVIMENTO_MAX = float(params.get('attrito_volvimento_max', self.COEFF_ATTRITO_VOLVIMENTO_MAX))
        self.COEFF_DRAG_ARIA_MIN = float(params.get('drag_aria_min', self.COEFF_DRAG_ARIA_MIN))
        self.COEFF_DRAG_ARIA_MAX = float(params.get('drag_aria_max', self.COEFF_DRAG_ARIA_MAX))
        self.DENSITA_ARIA = float(params.get('densita_aria', self.DENSITA_ARIA))
        self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK = float(params.get('coeff_restituzione_normale', self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK))
        self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK = float(params.get('coeff_restituzione_tangenziale', self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK))
        self.PENDENZA_GRADI = float(params.get('pendenza_gradi', self.PENDENZA_GRADI))
        self.PENDENZA_RAD = np.deg2rad(self.PENDENZA_GRADI)
        self.x0_sim = float(params.get('pos_x_iniziale', self.x0_sim))
        self.y_offset_sim = float(params.get('offset_y_iniziale', self.y_offset_sim))
        self.PASSO_SALVATAGGIO_OUTPUT = float(params.get('passo_salvataggio_output', self.PASSO_SALVATAGGIO_OUTPUT))
        
        # Se i parametri di pendenza sono cambiati E il pendio corrente è quello di default,
        # allora ricrea il pendio di default con la nuova pendenza.
        # Altrimenti, se un CSV è stato caricato, i suoi dati persistono.
        if self.x_pendio_data is not None and np.allclose(self.x_pendio_data, [0.0, 200.0]):
            self._crea_profilo_pendio_default()


    def initialize_simulation_state(self):
        self.reset_state()
        # Assicurati che il pendio sia sempre creato o caricato prima di generare i blocchi.
        # Se profilo_pendio_interpolator è None qui, significa che non è stato impostato
        # né da __init__, né da set_default_parameters, né da carica_dati_csv.
        if self.profilo_pendio_interpolator is None:
            self._crea_profilo_pendio_default()
        self.generare_blocchi()

    def _crea_profilo_pendio_default(self):
        x_default = np.array([0.0, 200.0])
        y_default = -np.tan(self.PENDENZA_RAD) * x_default
        r_normale_default = np.array([self.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK] * 2)
        r_tangenziale_default = np.array([self.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK] * 2)
        self.set_pendio_data(x_default, y_default, r_normale_default, r_tangenziale_default)
        self.x_pendio_min_csv = np.min(self.x_pendio_data) # Aggiorna i limiti anche per il default
        self.x_pendio_max_csv = np.max(self.x_pendio_data)


    def set_pendio_data(self, x_data, y_data, r_normale_data, r_tangenziale_data):
        if not (len(x_data) == len(y_data) == len(r_normale_data) == len(r_tangenziale_data) and len(x_data) >= 2):
            raise ValueError("I dati del pendio devono avere la stessa lunghezza e almeno 2 punti.")
        
        sorted_indices = np.argsort(x_data)
        self.x_pendio_data = np.array(x_data)[sorted_indices]
        self.y_pendio_data = np.array(y_data)[sorted_indices]
        self.r_normale_pendio_data = np.array(r_normale_data)[sorted_indices]
        self.r_tangenziale_pendio_data = np.array(r_tangenziale_data)[sorted_indices]

        self.profilo_pendio_interpolator = interp1d(self.x_pendio_data, self.y_pendio_data, kind='linear', fill_value="extrapolate")
        
        raw_interp_r_normale = interp1d(self.x_pendio_data, self.r_normale_pendio_data, kind='linear', fill_value="extrapolate")
        self.coeff_restituzione_normale_interpolator = lambda x: np.clip(raw_interp_r_normale(x), 0, 1)
        
        raw_interp_r_tangenziale = interp1d(self.x_pendio_data, self.r_tangenziale_pendio_data, kind='linear', fill_value="extrapolate")
        self.coeff_restituzione_tangenziale_interpolator = lambda x: np.clip(raw_interp_r_tangenziale(x), 0, 1)

        self.x_pendio_min_csv = np.min(self.x_pendio_data)
        self.x_pendio_max_csv = np.max(self.x_pendio_data)
        self.y_riferimento_energia = min(0, np.min(self.y_pendio_data))

    def carica_dati_csv(self, filepath):
        try:
            df = pd.read_csv(filepath)
            required_cols = ['x', 'y', 'r', 'et']
            if not all(col in df.columns for col in required_cols):
                raise ValueError(f"Il file CSV deve contenere le colonne: {', '.join(required_cols)}.")
            self.set_pendio_data(df['x'].values, df['y'].values, df['r'].values, df['et'].values)
            return True
        except Exception as e:
            print(f"Errore caricamento CSV: {e}")
            self.profilo_pendio_interpolator = None
            self.x_pendio_data = None # Resetta anche i dati del pendio
            self.y_pendio_data = None
            return False

    def generare_blocchi(self):
        self.blocchi = []
        for i in range(self.NUM_MASSI):
            if self.FORMA_BLOCCO == 'Sfera':
                volume = random.uniform(self.VOLUME_MIN, self.VOLUME_MAX)
                raggio = (volume * 3 / (4 * np.pi))**(1 / 3)
                dimensioni = [raggio]
            elif self.FORMA_BLOCCO == 'Prisma':
                l = random.uniform(self.LARGHEZZA_MIN, self.LARGHEZZA_MAX)
                h = random.uniform(self.ALTEZZA_MIN, self.ALTEZZA_MAX)
                p = random.uniform(self.PROFONDITA_MIN, self.PROFONDITA_MAX)
                dimensioni = [l, h, p]
            else:
                raise ValueError(f"Forma blocco '{self.FORMA_BLOCCO}' non riconosciuta.")

            coeff_attrito_volvimento = random.uniform(self.COEFF_ATTRITO_VOLVIMENTO_MIN, self.COEFF_ATTRITO_VOLVIMENTO_MAX)
            coeff_drag_aria = random.uniform(self.COEFF_DRAG_ARIA_MIN, self.COEFF_DRAG_ARIA_MAX)
            
            # Crea un blocco temporaneo per ottenere la dimensione di collisione
            blocco_temp = BloccoMasso(0, 0, 0, self.DENSITA, self.FORMA_BLOCCO, dimensioni, 0, 0)
            dist_collisione = blocco_temp.dist_centro_vertice

            # Modifica per posizionare i massi in cima o vicino alla posizione iniziale definita
            # Aggiungi una piccola variazione per evitare sovrapposizioni iniziali esatte
            x_iniziale = self.x0_sim + random.uniform(-0.5, 0.5) 
            
            # Calcola y_iniziale basandosi sull'interpolazione del pendio e l'offset
            y_pendio_at_x = self.profilo_pendio_interpolator(x_iniziale)
            y_iniziale = y_pendio_at_x + dist_collisione + self.y_offset_sim

            blocco = BloccoMasso(i + 1, x_iniziale, y_iniziale, self.DENSITA, self.FORMA_BLOCCO, dimensioni, coeff_attrito_volvimento, coeff_drag_aria)
            
            # Velocità iniziali casuali o predefinite (puoi affinare queste se necessario)
            blocco.velocita_x = [random.uniform(0.1, 1.0)]
            blocco.velocita_y = [random.uniform(-0.1, -0.5)]
            blocco.velocita_angolare = [random.uniform(-0.5, 0.5)]
            
            self.blocchi.append(blocco)

    def calcola_energie(self, blocco):
        Ek_trasl = 0.5 * blocco.massa * (blocco.velocita_x[-1]**2 + blocco.velocita_y[-1]**2)
        Ek_rot = 0.5 * blocco.momento_inerzia * blocco.velocita_angolare[-1]**2
        Ek_tot = Ek_trasl + Ek_rot
        Ep = blocco.massa * self.G * (blocco.posizioni_y[-1] - self.y_riferimento_energia)
        return Ek_tot, Ep

    def step_simulazione(self):
        self.tempo_corrente += self.DT
        self.frame_corrente += 1

        active_blocks_count = 0
        for blocco in self.blocchi:
            if not blocco.is_active:
                continue

            active_blocks_count += 1
            x, y = blocco.posizioni_x[-1], blocco.posizioni_y[-1]
            vx, vy = blocco.velocita_x[-1], blocco.velocita_y[-1]
            omega = blocco.velocita_angolare[-1]

            # Calcolo delle forze
            # Forza di gravità
            Fx_g = 0
            Fy_g = -blocco.massa * self.G

            # Forza di resistenza dell'aria (Drag)
            # Dipende dall'area frontale e dalla velocità relativa all'aria
            v_rel = np.sqrt(vx**2 + vy**2)
            if v_rel > 0.01: # Evita divisione per zero e drag insignificante a basse velocità
                # Drag proporzionale al quadrato della velocità
                # Direzione opposta al vettore velocità
                Fd_x = -0.5 * self.DENSITA_ARIA * v_rel * vx * blocco.area_frontale * blocco.coeff_drag_aria
                Fd_y = -0.5 * self.DENSITA_ARIA * v_rel * vy * blocco.area_frontale * blocco.coeff_drag_aria
            else:
                Fd_x = 0
                Fd_y = 0

            Fx_net = Fx_g + Fd_x
            Fy_net = Fy_g + Fd_y

            # Integrazione delle equazioni del moto (Euler semi-implicito o Velocity Verlet)
            ax = Fx_net / blocco.massa
            ay = Fy_net / blocco.massa

            vx_new = vx + ax * self.DT
            vy_new = vy + ay * self.DT

            x_new = x + vx_new * self.DT
            y_new = y + vy_new * self.DT

            # Contatto con il pendio
            y_pendio = self.profilo_pendio_interpolator(x_new)
            
            # Calcola la pendenza locale del pendio nel punto x_new
            # Questo è un valore approssimato per pendii discreti/interpolati linearmente
            # Per una vera derivata, servirebbe una funzione di interpolazione che permetta la derivazione
            # o calcolare la pendenza del segmento su cui si trova x_new.
            # Per semplicità, usiamo la pendenza della linea che collega il punto attuale e un punto leggermente avanti.
            # Metodo più robusto: trovare i due punti del pendio tra cui si trova x_new
            idx = np.searchsorted(self.x_pendio_data, x_new)
            if idx == 0:
                pendenza_locale_rad = np.arctan2(self.y_pendio_data[1] - self.y_pendio_data[0], self.x_pendio_data[1] - self.x_pendio_data[0])
            elif idx == len(self.x_pendio_data):
                pendenza_locale_rad = np.arctan2(self.y_pendio_data[-1] - self.y_pendio_data[-2], self.x_pendio_data[-1] - self.x_pendio_data[-2])
            else:
                # Interpolazione lineare, la pendenza è costante tra due punti dati
                pendenza_locale_rad = np.arctan2(self.y_pendio_data[idx] - self.y_pendio_data[idx-1], self.x_pendio_data[idx] - self.x_pendio_data[idx-1])
            
            # Vettore normale al pendio
            norm_x = np.sin(pendenza_locale_rad)
            norm_y = -np.cos(pendenza_locale_rad)

            # Vettore tangenziale al pendio
            tan_x = np.cos(pendenza_locale_rad)
            tan_y = np.sin(pendenza_locale_rad)

            # Distanza dal pendio
            dist_dal_pendio = y_new - y_pendio

            # Raggio del blocco o dimensione rilevante per collisione
            raggio_collisione = blocco.dist_centro_vertice

            # Condizione di collisione o penetrazione
            if dist_dal_pendio < raggio_collisione: # Se il blocco penetra il pendio
                # Rimbalzo
                # Componenti della velocità normale e tangenziale al pendio
                v_dot_n = vx_new * norm_x + vy_new * norm_y
                v_dot_t = vx_new * tan_x + vy_new * tan_y

                # Coefficienti di restituzione (possono variare lungo il pendio se caricati da CSV)
                coeff_restituzione_normale = self.coeff_restituzione_normale_interpolator(x_new)
                coeff_restituzione_tangenziale = self.coeff_restituzione_tangenziale_interpolator(x_new)

                # Nuove velocità dopo il rimbalzo
                v_n_new = -coeff_restituzione_normale * v_dot_n
                v_t_new = coeff_restituzione_tangenziale * v_dot_t

                # Ricombina le componenti per ottenere le nuove velocità cartesiane
                vx_new = v_n_new * norm_x + v_t_new * tan_x
                vy_new = v_n_new * norm_y + v_t_new * tan_y
                
                # Applica attrito al rotolamento (semplificato)
                # Forza di attrito opposta alla direzione del moto
                # Non è un attrito di rotolamento completo, ma una decelerazione angolare
                blocco.velocita_angolare.append(omega * (1 - blocco.coeff_attrito_volvimento))

                # Evita che il blocco sprofondi nel pendio (lo riposiziona sulla superficie)
                y_new = y_pendio + raggio_collisione + 0.001 # Piccolo offset per evitare nuove collisioni immediate

            blocco.posizioni_x.append(x_new)
            blocco.posizioni_y.append(y_new)
            blocco.velocita_x.append(vx_new)
            blocco.velocita_y.append(vy_new)

            # Aggiungi un punto per l'angolo e la velocità angolare anche se non cambiano ad ogni passo
            if len(blocco.angolo) < len(blocco.posizioni_x):
                blocco.angolo.append(blocco.angolo[-1] + blocco.velocita_angolare[-1] * self.DT)
            if len(blocco.velocita_angolare) < len(blocco.posizioni_x):
                blocco.velocita_angolare.append(omega) # Mantieni costante se non c'è interazione

            # Controllo delle condizioni di arresto
            velocita_mag = np.sqrt(vx_new**2 + vy_new**2)
            if (velocita_mag < SOGLIA_VELOCITA_ARRESTO and abs(omega) < 0.1 and dist_dal_pendio < SOGLIA_DISTANZA_PENDIO_ARRESTO):
                blocco.is_active = False
                blocco.final_status = "Arrestato sul pendio"
                blocco.final_time = self.tempo_corrente
            elif y_new < y_pendio - SOGLIA_PROFONDITA_FUORI_PENDIO_Y: # Blocco molto sotto il pendio
                blocco.is_active = False
                blocco.final_status = "Fuori dal pendio (profondità)"
                blocco.final_time = self.tempo_corrente
            elif x_new > self.x_pendio_max_csv + SOGLIA_MARGINE_FINE_PENDIO_X: # Blocco oltre la fine del pendio
                 blocco.is_active = False
                 blocco.final_status = "Fuori dal pendio (fine orizzontale)"
                 blocco.final_time = self.tempo_corrente

            # Salva dati di output a intervalli regolari
            if self.tempo_corrente - blocco.last_saved_time >= self.PASSO_SALVATAGGIO_OUTPUT or not blocco.is_active:
                Ek_tot, Ep = self.calcola_energie(blocco)
                blocco.output_data.append({
                    'Tempo (s)': self.tempo_corrente,
                    'Pos_X (m)': x_new,
                    'Pos_Y (m)': y_new,
                    'Vel_X (m/s)': vx_new,
                    'Vel_Y (m/s)': vy_new,
                    'Vel_Ang (rad/s)': omega,
                    'Ek_Tot (J)': Ek_tot,
                    'Ep (J)': Ep,
                    'E_Tot (J)': Ek_tot + Ep,
                    'Stato Finale': blocco.final_status
                })
                blocco.last_saved_time = self.tempo_corrente
                blocco.last_saved_x = x_new # Aggiorna anche l'ultima posizione X salvata

        if active_blocks_count == 0 or self.tempo_corrente >= self.TEMPO_TOTALE:
            for blocco in self.blocchi:
                if blocco.is_active: # Cattura i massi ancora attivi alla fine del tempo
                    blocco.is_active = False
                    blocco.final_status = "Tempo massimo raggiunto"
                    blocco.final_time = self.tempo_corrente
                    # Assicurati che l'ultimo punto sia salvato
                    Ek_tot, Ep = self.calcola_energie(blocco)
                    blocco.output_data.append({
                        'Tempo (s)': self.tempo_corrente,
                        'Pos_X (m)': blocco.posizioni_x[-1],
                        'Pos_Y (m)': blocco.posizioni_y[-1],
                        'Vel_X (m/s)': blocco.velocita_x[-1],
                        'Vel_Y (m/s)': blocco.velocita_y[-1],
                        'Vel_Ang (rad/s)': blocco.velocita_angolare[-1],
                        'Ek_Tot (J)': Ek_tot,
                        'Ep (J)': Ep,
                        'E_Tot (J)': Ek_tot + Ep,
                        'Stato Finale': blocco.final_status
                    })
            return False # La simulazione è terminata
        return True # La simulazione continua

    def salva_output_csv(self, filename):
        all_blocks_data = []
        for blocco in self.blocchi:
            df_block = pd.DataFrame(blocco.output_data)
            df_block['ID_Masso'] = blocco.id
            df_block['Forma'] = blocco.forma
            df_block['Massa (kg)'] = blocco.massa
            all_blocks_data.append(df_block)
        
        if not all_blocks_data:
            return False

        df_final = pd.concat(all_blocks_data, ignore_index=True)
        try:
            df_final.to_csv(filename, index=False)
            return True
        except Exception as e:
            print(f"Errore durante il salvataggio del CSV: {e}")
            return False

# --- Interfaccia Utente (GUI) ---
class MainWindow(QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Simulatore Caduta Massi")
        self.setGeometry(100, 100, 1200, 800)

        self.simulatore = SimulatoreCadutaMassi(parent_gui=self)
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_simulation)
        self.is_running = False

        self.central_widget = QWidget()
        self.setCentralWidget(self.central_widget)
        self.main_layout = QHBoxLayout(self.central_widget)

        self.create_control_panel()
        self.create_simulation_plot()
        
        # Inizializza la simulazione e disegna lo stato iniziale
        self.simulatore.initialize_simulation_state()
        self.draw_initial_state()
        self.update_input_fields() # Aggiorna i campi input con i valori di default

    def create_control_panel(self):
        self.control_panel = QScrollArea()
        self.control_panel.setWidgetResizable(True)
        self.control_panel.setFixedWidth(400)
        
        self.control_content_widget = QWidget()
        self.control_panel.setWidget(self.control_content_widget)
        self.control_layout = QVBoxLayout(self.control_content_widget)

        # Sezione Parametri di Simulazione
        self.control_layout.addWidget(QLabel("<h2>Parametri di Simulazione</h2>"))
        self.form_layout = QFormLayout()
        self.input_fields = {}

        # Parametri Generali
        self.add_input_field("gravita", "Gravità (m/s²):", str(G_DEFAULT), QDoubleValidator())
        self.add_input_field("passo_temporale", "Passo Temporale (s):", str(DT_DEFAULT), QDoubleValidator())
        self.add_input_field("tempo_totale", "Tempo Totale Max (s):", str(TEMPO_TOTALE_MAX_DEFAULT), QDoubleValidator())
        self.add_input_field("passo_salvataggio_output", "Passo Salv. Output (s):", str(PASSO_SALVATAGGIO_OUTPUT_DEFAULT), QDoubleValidator())

        # Parametri Massi
        self.add_input_field("num_massi", "Numero Massi:", str(NUM_MASSI_DEFAULT), QIntValidator())
        self.add_input_field("densita", "Densità Masso (kg/m³):", str(DENSITA_DEFAULT), QDoubleValidator())
        
        self.form_layout.addRow(QLabel("Forma Blocco:"), self.create_forma_blocco_selector())
        self.add_input_field("volume_min", "Volume Min (m³):", str(VOLUME_MIN_DEFAULT), QDoubleValidator())
        self.add_input_field("volume_max", "Volume Max (m³):", str(VOLUME_MAX_DEFAULT), QDoubleValidator())
        self.add_input_field("larghezza_min", "Larghezza Min (m):", str(LARGHEZZA_MIN_DEFAULT), QDoubleValidator())
        self.add_input_field("larghezza_max", "Larghezza Max (m):", str(LARGHEZZA_MAX_DEFAULT), QDoubleValidator())
        self.add_input_field("altezza_min", "Altezza Min (m):", str(ALTEZZA_MIN_DEFAULT), QDoubleValidator())
        self.add_input_field("altezza_max", "Altezza Max (m):", str(ALTEZZA_MAX_DEFAULT), QDoubleValidator())
        self.add_input_field("profondita_min", "Profondità Min (m):", str(PROFONDITA_MIN_DEFAULT), QDoubleValidator())
        self.add_input_field("profondita_max", "Profondità Max (m):", str(PROFONDITA_MAX_DEFAULT), QDoubleValidator())

        self.add_input_field("attrito_volvimento_min", "Attrito Volvimento Min:", str(COEFF_ATTRITO_VOLVIMENTO_MIN_DEFAULT), QDoubleValidator())
        self.add_input_field("attrito_volvimento_max", "Attrito Volvimento Max:", str(COEFF_ATTRITO_VOLVIMENTO_MAX_DEFAULT), QDoubleValidator())
        self.add_input_field("drag_aria_min", "Drag Aria Min:", str(COEFF_DRAG_ARIA_MIN_DEFAULT), QDoubleValidator())
        self.add_input_field("drag_aria_max", "Drag Aria Max:", str(COEFF_DRAG_ARIA_MAX_DEFAULT), QDoubleValidator())
        self.add_input_field("densita_aria", "Densità Aria (kg/m³):", str(DENSITA_ARIA_DEFAULT), QDoubleValidator())
        
        # Parametri Pendio (per default o fallback)
        self.add_input_field("pendenza_gradi", "Pendenza Pendio (°):", str(PENDENZA_GRADI_DEFAULT), QDoubleValidator())
        self.add_input_field("pos_x_iniziale", "Pos. X Iniziale (m):", str(X0_DEFAULT), QDoubleValidator())
        self.add_input_field("offset_y_iniziale", "Offset Y Iniziale (m):", str(Y_INIZIALE_OFFSET_DEFAULT), QDoubleValidator())
        self.add_input_field("coeff_restituzione_normale", "Coeff. Rest. Normale:", str(COEFF_RESTITUZIONE_NORMALE_DEFAULT), QDoubleValidator())
        self.add_input_field("coeff_restituzione_tangenziale", "Coeff. Rest. Tangenziale:", str(COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT), QDoubleValidator())


        self.control_layout.addLayout(self.form_layout)

        # Pulsanti di controllo
        self.button_layout = QHBoxLayout()
        
        # Nuovo pulsante Applica Parametri
        self.apply_params_button = QPushButton("Applica Parametri")
        self.apply_params_button.clicked.connect(self.apply_parameters_and_reset_view)
        self.button_layout.addWidget(self.apply_params_button)

        self.start_button = QPushButton("Avvia Simulazione")
        self.start_button.clicked.connect(self.start_simulation)
        self.button_layout.addWidget(self.start_button)

        self.stop_button = QPushButton("Ferma Simulazione")
        self.stop_button.clicked.connect(self.stop_simulation)
        self.stop_button.setEnabled(False)
        self.button_layout.addWidget(self.stop_button)

        self.reset_button = QPushButton("Reset Simulazione")
        self.reset_button.clicked.connect(self.reset_simulation)
        self.button_layout.addWidget(self.reset_button)
        
        self.control_layout.addLayout(self.button_layout)

        # Pulsanti per Caricamento/Salvataggio CSV
        self.file_button_layout = QHBoxLayout()
        self.load_csv_button = QPushButton("Carica Pendio da CSV")
        self.load_csv_button.clicked.connect(self.load_pendio_csv)
        self.file_button_layout.addWidget(self.load_csv_button)
        
        self.csv_status_label = QLabel("Nessun CSV caricato")
        self.file_button_layout.addWidget(self.csv_status_label)

        self.save_output_button = QPushButton("Salva Dati Output CSV")
        self.save_output_button.clicked.connect(self.save_output_csv_file)
        self.file_button_layout.addWidget(self.save_output_button)

        self.control_layout.addLayout(self.file_button_layout)

        # Nuovo pulsante per il report finale
        self.report_button = QPushButton("Mostra Report Finale")
        self.report_button.clicked.connect(self.show_final_results_summary)
        self.control_layout.addWidget(self.report_button)

        self.control_layout.addStretch() # Spinge tutto in alto

        self.main_layout.addWidget(self.control_panel)

    def add_input_field(self, name, label_text, default_value, validator=None):
        line_edit = QLineEdit(default_value)
        if validator:
            line_edit.setValidator(validator)
        self.form_layout.addRow(QLabel(label_text), line_edit)
        self.input_fields[name] = line_edit

    def create_forma_blocco_selector(self):
        self.forma_blocco_combo = QComboBox()
        self.forma_blocco_combo.addItems(['Sfera', 'Prisma'])
        self.forma_blocco_combo.setCurrentText(FORMA_BLOCCO_DEFAULT)
        self.forma_blocco_combo.currentTextChanged.connect(self.update_dimensions_visibility)
        return self.forma_blocco_combo

    def update_dimensions_visibility(self):
        selected_shape = self.forma_blocco_combo.currentText()
        
        widgets_to_toggle = {
            'volume_min': [], 'volume_max': [],
            'larghezza_min': [], 'larghezza_max': [],
            'altezza_min': [], 'altezza_max': [],
            'profondita_min': [], 'profondita_max': []
        }

        # Popola il dizionario con tutti i widget della riga
        for i in range(self.form_layout.rowCount()):
            label_item = self.form_layout.itemAt(i, QFormLayout.LabelRole)
            field_item = self.form_layout.itemAt(i, QFormLayout.FieldRole)
            
            if label_item and field_item:
                label_widget = label_item.widget()
                field_widget = field_item.widget()
                
                # Cerca quale campo input corrisponde a questo field_widget
                for name, line_edit_widget in self.input_fields.items():
                    if line_edit_widget == field_widget:
                        if name in widgets_to_toggle:
                            widgets_to_toggle[name].extend([label_widget, field_widget])
                        break

        # Ora applica la visibilità in base alla forma selezionata
        for name, widgets in widgets_to_toggle.items():
            if 'volume' in name:
                visible = (selected_shape == 'Sfera')
            elif 'larghezz' in name or 'altezza' in name or 'profondita' in name:
                visible = (selected_shape == 'Prisma')
            else:
                visible = True # Di default visibile se non specificato

            for widget in widgets:
                if widget: # Assicurati che il widget esista
                    widget.setVisible(visible)
        
        self.control_content_widget.adjustSize()
        self.control_panel.widget().layout().invalidate() # Forza il ricalcolo del layout
        
    def update_input_fields(self):
        # Popola i campi di input con i valori correnti del simulatore
        self.input_fields['gravita'].setText(str(self.simulatore.G))
        self.input_fields['passo_temporale'].setText(str(self.simulatore.DT))
        self.input_fields['tempo_totale'].setText(str(self.simulatore.TEMPO_TOTALE))
        self.input_fields['passo_salvataggio_output'].setText(str(self.simulatore.PASSO_SALVATAGGIO_OUTPUT))
        self.input_fields['num_massi'].setText(str(self.simulatore.NUM_MASSI))
        self.input_fields['densita'].setText(str(self.simulatore.DENSITA))
        self.forma_blocco_combo.setCurrentText(self.simulatore.FORMA_BLOCCO)
        self.input_fields['volume_min'].setText(str(self.simulatore.VOLUME_MIN))
        self.input_fields['volume_max'].setText(str(self.simulatore.VOLUME_MAX))
        self.input_fields['larghezza_min'].setText(str(self.simulatore.LARGHEZZA_MIN))
        self.input_fields['larghezza_max'].setText(str(self.simulatore.LARGHEZZA_MAX))
        self.input_fields['altezza_min'].setText(str(self.simulatore.ALTEZZA_MIN))
        self.input_fields['altezza_max'].setText(str(self.simulatore.ALTEZZA_MAX))
        self.input_fields['profondita_min'].setText(str(self.simulatore.PROFONDITA_MIN))
        self.input_fields['profondita_max'].setText(str(self.simulatore.PROFONDITA_MAX))
        self.input_fields['attrito_volvimento_min'].setText(str(self.simulatore.COEFF_ATTRITO_VOLVIMENTO_MIN))
        self.input_fields['attrito_volvimento_max'].setText(str(self.simulatore.COEFF_ATTRITO_VOLVIMENTO_MAX))
        self.input_fields['drag_aria_min'].setText(str(self.simulatore.COEFF_DRAG_ARIA_MIN))
        self.input_fields['drag_aria_max'].setText(str(self.simulatore.COEFF_DRAG_ARIA_MAX))
        self.input_fields['densita_aria'].setText(str(self.simulatore.DENSITA_ARIA))
        self.input_fields['pendenza_gradi'].setText(str(self.simulatore.PENDENZA_GRADI))
        self.input_fields['pos_x_iniziale'].setText(str(self.simulatore.x0_sim))
        self.input_fields['offset_y_iniziale'].setText(str(self.simulatore.y_offset_sim))
        self.input_fields['coeff_restituzione_normale'].setText(str(self.simulatore.COEFF_RESTITUZIONE_NORMALE_DEFAULT_FALLBACK))
        self.input_fields['coeff_restituzione_tangenziale'].setText(str(self.simulatore.COEFF_RESTITUZIONE_TANGENZIALE_DEFAULT_FALLBACK))
        
        self.update_dimensions_visibility() # Assicurati che i campi di dimensione siano visibili correttamente

    def apply_parameters_and_reset_view(self):
        try:
            params = {
                'gravita': self.input_fields['gravita'].text(),
                'passo_temporale': self.input_fields['passo_temporale'].text(),
                'tempo_totale': self.input_fields['tempo_totale'].text(),
                'passo_salvataggio_output': self.input_fields['passo_salvataggio_output'].text(),
                'num_massi': self.input_fields['num_massi'].text(),
                'densita': self.input_fields['densita'].text(),
                'forma_blocco': self.forma_blocco_combo.currentText(),
                'volume_min': self.input_fields['volume_min'].text(),
                'volume_max': self.input_fields['volume_max'].text(),
                'larghezza_min': self.input_fields['larghezza_min'].text(),
                'larghezza_max': self.input_fields['larghezza_max'].text(),
                'altezza_min': self.input_fields['altezza_min'].text(),
                'altezza_max': self.input_fields['altezza_max'].text(),
                'profondita_min': self.input_fields['profondita_min'].text(),
                'profondita_max': self.input_fields['profondita_max'].text(),
                'attrito_volvimento_min': self.input_fields['attrito_volvimento_min'].text(),
                'attrito_volvimento_max': self.input_fields['attrito_volvimento_max'].text(),
                'drag_aria_min': self.input_fields['drag_aria_min'].text(),
                'drag_aria_max': self.input_fields['drag_aria_max'].text(),
                'densita_aria': self.input_fields['densita_aria'].text(),
                'pendenza_gradi': self.input_fields['pendenza_gradi'].text(),
                'pos_x_iniziale': self.input_fields['pos_x_iniziale'].text(),
                'offset_y_iniziale': self.input_fields['offset_y_iniziale'].text(),
                'coeff_restituzione_normale': self.input_fields['coeff_restituzione_normale'].text(),
                'coeff_restituzione_tangenziale': self.input_fields['coeff_restituzione_tangenziale'].text()
            }
            self.simulatore.set_parameters_from_dict(params)
            self.simulatore.initialize_simulation_state() # Rigenera i massi con i nuovi parametri
            self.draw_initial_state() # Aggiorna la vista
            QMessageBox.information(self, "Parametri Aggiornati", "I parametri della simulazione sono stati applicati e i massi rigenerati.")
        except ValueError as e:
            QMessageBox.critical(self, "Errore Parametri", f"Controlla i valori inseriti: {e}")
        except Exception as e:
            QMessageBox.critical(self, "Errore Applicazione", f"Si è verificato un errore inatteso durante l'applicazione dei parametri: {e}")


    def create_simulation_plot(self):
        self.figure, self.ax = plt.subplots(figsize=(10, 8))
        self.canvas = FigureCanvas(self.figure)
        self.main_layout.addWidget(self.canvas)
        self.ax.set_xlabel("Posizione X (m)")
        self.ax.set_ylabel("Posizione Y (m)")
        self.ax.set_title("Simulazione Caduta Massi")
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.grid(True)

        self.pendio_line, = self.ax.plot([], [], 'k-', linewidth=2, label="Pendio")
        self.massi_grafici = []
        self.tracce_grafiche = []

    def draw_initial_state(self):
        self.ax.clear()
        self.ax.set_xlabel("Posizione X (m)")
        self.ax.set_ylabel("Posizione Y (m)")
        self.ax.set_title("Simulazione Caduta Massi")
        self.ax.set_aspect('equal', adjustable='box')
        self.ax.grid(True)

        # Disegna il pendio
        if self.simulatore.x_pendio_data is not None:
            self.pendio_line, = self.ax.plot(self.simulatore.x_pendio_data, self.simulatore.y_pendio_data, 'k-', linewidth=2, label="Pendio")
        else: # Disegna il pendio di default se non caricato da CSV
            x_default = np.array([0.0, 200.0])
            y_default = -np.tan(self.simulatore.PENDENZA_RAD) * x_default
            self.pendio_line, = self.ax.plot(x_default, y_default, 'k-', linewidth=2, label="Pendio Default")

        # Disegna i massi iniziali
        self.massi_grafici = []
        self.tracce_grafiche = []
        for blocco in self.simulatore.blocchi:
            marker_size = blocco.dist_centro_vertice * 10  # Scala per visibilità
            masso_g, = self.ax.plot([blocco.posizioni_x[0]], [blocco.posizioni_y[0]], 'o', markersize=marker_size, color=blocco.colore)
            traccia_g, = self.ax.plot([], [], '--', color=blocco.colore, linewidth=1)
            self.massi_grafici.append(masso_g)
            self.tracce_grafiche.append(traccia_g)
        
        self.ax.relim()
        self.ax.autoscale_view()
        self.canvas.draw()
        
    def start_simulation(self):
        try:
            # Assicurati che i parametri siano applicati e i massi generati prima di iniziare
            self.apply_parameters_and_reset_view() # Questo chiamerà initialize_simulation_state e draw_initial_state

            if not self.simulatore.blocchi:
                QMessageBox.warning(self, "Nessun Masso", "Genera i massi o carica un pendio prima di avviare la simulazione.")
                return

            self.start_button.setEnabled(False)
            self.stop_button.setEnabled(True)
            self.reset_button.setEnabled(False)
            self.load_csv_button.setEnabled(False)
            self.save_output_button.setEnabled(False)
            self.report_button.setEnabled(False) # Disabilita il report durante la simulazione
            self.apply_params_button.setEnabled(False) # Disabilita applica parametri durante la simulazione


            self.simulatore.reset_state() # Resetta solo lo stato di runtime (tempo, frame)
            # Nota: i massi sono già stati rigenerati da apply_parameters_and_reset_view
            
            self.timer.start(int(self.simulatore.DT * 1000)) # Converti secondi in millisecondi
            self.is_running = True
        except ValueError as e:
            QMessageBox.critical(self, "Errore Parametri", f"Controlla i valori inseriti: {e}")
        except Exception as e:
            QMessageBox.critical(self, "Errore Avvio", f"Si è verificato un errore inatteso all'avvio: {e}")

    def stop_simulation(self):
        self.timer.stop()
        self.is_running = False
        self.start_button.setEnabled(True)
        self.stop_button.setEnabled(False)
        self.reset_button.setEnabled(True)
        self.load_csv_button.setEnabled(True)
        self.save_output_button.setEnabled(True)
        self.report_button.setEnabled(True) # Riabilita il report a simulazione ferma
        self.apply_params_button.setEnabled(True) # Riabilita applica parametri


    def reset_simulation(self):
        self.stop_simulation()
        self.simulatore.set_default_parameters() # Resetta i parametri interni del simulatore
        self.update_input_fields() # Aggiorna i campi input della GUI ai default
        self.csv_status_label.setText("Nessun CSV caricato")
        
        # Resetta lo stato interno del simulatore e rigenera i massi
        # set_default_parameters già chiama _crea_profilo_pendio_default
        self.simulatore.initialize_simulation_state() # Questo genererà nuovi massi con i default
        self.draw_initial_state() # E disegnerà il nuovo stato iniziale
        QMessageBox.information(self, "Reset Completato", "La simulazione è stata resettata ai valori predefiniti.")

    def update_simulation(self):
        if self.simulatore.step_simulazione():
            for i, blocco in enumerate(self.simulatore.blocchi):
                if blocco.is_active:
                    # FIX: Passare liste a set_data anche per singoli valori
                    self.massi_grafici[i].set_data([blocco.posizioni_x[-1]], [blocco.posizioni_y[-1]])
                    self.tracce_grafiche[i].set_data(blocco.posizioni_x, blocco.posizioni_y)
            self.ax.relim()
            self.ax.autoscale_view()
            self.canvas.draw()
        else:
            self.stop_simulation()
            QMessageBox.information(self, "Simulazione Terminata", "La simulazione è giunta al termine.")


    def load_pendio_csv(self):
        fname, _ = QFileDialog.getOpenFileName(self, "Carica CSV Pendio", "", "CSV Files (*.csv)")
        if fname:
            if self.simulatore.carica_dati_csv(fname):
                self.csv_status_label.setText(f"Caricato: {os.path.basename(fname)}")
                # Dopo aver caricato il CSV, reinizializza lo stato per usare il nuovo pendio
                self.simulatore.initialize_simulation_state() 
                self.draw_initial_state()
                QMessageBox.information(self, "Successo", "File CSV del pendio caricato con successo.")
            else:
                QMessageBox.warning(self, "Errore", "Impossibile caricare o leggere il file CSV del pendio.")

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
                QMessageBox.warning(self, "Errore Salvataggio", f"Si è verificato un errore durante il salvataggio dei dati in:\n{fileName}")

    def show_final_results_summary(self):
        if not self.simulatore.blocchi or all(blocco.final_status == "In simulazione" for blocco in self.simulatore.blocchi):
            QMessageBox.information(self, "Nessun Risultato", "La simulazione non è stata ancora eseguita o completata.")
            return

        final_statuses = [blocco.final_status for blocco in self.simulatore.blocchi]
        status_counts = pd.Series(final_statuses).value_counts()

        labels = status_counts.index.tolist()
        sizes = status_counts.values.tolist()

        # Dati per la tabella riassuntiva
        table_data = []
        for blocco in self.simulatore.blocchi:
            last_x = blocco.posizioni_x[-1] if blocco.posizioni_x else 'N/A'
            last_y = blocco.posizioni_y[-1] if blocco.posizioni_y else 'N/A'
            final_time = f"{blocco.final_time:.2f} s" if blocco.final_time is not None else 'N/A'
            
            # Recupera l'ultima energia totale se disponibile, altrimenti N/A
            e_tot_finale = 'N/A'
            if blocco.output_data:
                # Trova l'ultimo record che contiene 'E_Tot (J)'
                for record in reversed(blocco.output_data):
                    if 'E_Tot (J)' in record and record['E_Tot (J)'] is not None:
                        e_tot_finale = f"{record['E_Tot (J)']:.2f} J"
                        break

            table_data.append([
                str(blocco.id),
                blocco.forma,
                f"{blocco.massa:.2f} kg",
                f"{last_x:.2f}",
                f"{last_y:.2f}",
                final_time,
                blocco.final_status,
                e_tot_finale
            ])
        
        # Ordina la tabella per ID masso
        table_data.sort(key=lambda x: int(x[0]))

        # Creazione dei subplot: uno per il grafico a torta, uno per la tabella
        fig_summary, (ax_pie, ax_table) = plt.subplots(1, 2, figsize=(15, 7))

        # Grafico a torta
        if not labels:
            ax_pie.text(0.5, 0.5, "Nessun dato per il grafico", horizontalalignment='center', verticalalignment='center', transform=ax_pie.transAxes)
            ax_pie.axis('off')
        else:
            colors = plt.cm.Paired(np.arange(len(labels)))
            ax_pie.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90, colors=colors)
            ax_pie.axis('equal')  # Assicura che il grafico a torta sia circolare.
        ax_pie.set_title('Distribuzione Stati Finali Massi')

        # Tabella riassuntiva
        ax_table.axis('off') # Nasconde gli assi per la tabella
        ax_table.set_title('Riepilogo Risultati Finali Massi', y=1.05) 

        if not table_data:
            ax_table.text(0.5, 0.5, "Nessun dato di simulazione da mostrare", horizontalalignment='center', verticalalignment='center', transform=ax_table.transAxes)
        else:
            col_labels = ["ID Masso", "Forma", "Massa", "Pos X Finale (m)", "Pos Y Finale (m)", "Tempo Finale", "Stato Finale", "E Totale Finale"]
            table = ax_table.table(
                cellText=table_data,
                colLabels=col_labels,
                loc='center', 
                cellLoc='center' 
            )
            table.auto_set_font_size(False)
            table.set_fontsize(10) 
            table.scale(1, 1.2) 

        plt.tight_layout() 
        plt.show()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main_window = MainWindow()
    main_window.show()
    sys.exit(app.exec_())
