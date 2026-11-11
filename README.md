# Biomass_production — crescita di *Nostoc* in coltura continua

Analisi bayesiana di un esperimento sulla produzione di biomassa di un cianobatterio
del genere *Nostoc*: come i fattori di coltivazione influenzano la crescita nel tempo,
misurata attraverso la densità ottica della coltura.

Progetto di gruppo di **Giorgia Caicchiolo, Marco Nascimben, Riccardo Pati e Luca
Varotto**, M.Sc. in Scienze Statistiche, Università di Padova, A.A. 2025/2026.

## Struttura

```
Biomass_production/
├── Codice/
│   ├── Analisi_esplorativa.R      preparazione dei dati, esplorativa, GAM -> Figures/
│   └── Stan/
│       ├── modello.R              stima del modello finale con cmdstanr e diagnostica
│       ├── modello_finale.stan    il modello consegnato
│       ├── ottimizzazione.R       condizioni ottimali dai campioni a posteriori
│       ├── README.md              il percorso dai cinque modelli di prova a quello finale
│       ├── Prove/                 i cinque modelli provati, numerati in ordine
│       │   ├── 1_base_ar1.stan ... 5_eteroschedastico_quadratico.stan
│       │   └── fit_model.R, diagnostica.R
│       └── fits/                  oggetti stimati, non versionati
├── Figures/                       grafici: esplorativa, diagnostica MCMC e modello
├── Note/                          appunti e verbali di gruppo
├── Report/
│   ├── Relazione.pdf              la relazione finale
│   └── Slides.pdf                 le slide di presentazione
├── Consegna.pdf                   la traccia del progetto
├── .gitignore
└── README.md
```

## Come si esegue

Gli script ricavano i percorsi dalla posizione del file aperto in RStudio, quindi
funzionano senza impostare la working directory a mano. L'ordine è:

1. `Codice/Analisi_esplorativa.R` — legge `nostoc.xlsx`, scrive `dati_modificati.Rdata` e i
   grafici esplorativi in `Figures/`;
2. `Codice/Stan/modello.R` — compila `modello_finale.stan` con **cmdstanr**, stima il
   modello, salva l'oggetto in `Codice/Stan/fits/` e la diagnostica in
   `Figures/`;
3. `Codice/Stan/ottimizzazione.R` — rilegge il modello stimato e cerca le condizioni
   che massimizzano la resa.

Gli script in `Codice/Stan/Prove/` ripercorrono i modelli intermedi e servono solo per
riprodurre quel confronto: `Codice/Stan/README.md` racconta cosa aggiunge ciascuno e
perché si è arrivati al modello finale.

## Risultati

**L'esperimento.** 53 colture di *Nostoc*, seguite settimana per settimana, con tre
fattori di coltivazione: l'intensità luminosa `I`, la durata dell'esposizione `D` in
giorni e la concentrazione di fosforo `P`. La risposta è la densità ottica (OD), una
misura adimensionale della biomassa prodotta. Le 17 combinazioni sperimentali seguono
la geometria di un *Face-Centered Central Composite Design*, con repliche da 2 a 4 per
punto e 8 al centro del dominio — l'iper-replicazione del punto centrale serve a
stimare bene la varianza dell'errore.

**Il modello.** Risposta normale, media descritta da una curva di **Gompertz**, termini
quadratici nei predittori lineari, un effetto casuale sul parametro di velocità di
crescita ed eteroschedasticità *power-law*. Stimato in Stan con cmdstanr. La scelta
bayesiana non è decorativa: serve a propagare l'incertezza nelle previsioni e a
confrontare direttamente le combinazioni sperimentali senza incorrere nel problema dei
test multipli, cosa che con la *Response Surface Methodology* classica non sarebbe
stata possibile, dato che qui la risposta è una serie storica e non uno scalare.

**Cosa guida la crescita.** La biomassa massima raggiungibile dipende soprattutto da
`I` e `P` e dalla loro interazione, mentre la velocità con cui quella biomassa viene
raggiunta è governata in larga parte dall'intensità luminosa. Le diagnostiche MCMC e i
*posterior predictive check* mostrano che il modello riproduce bene le caratteristiche
dei dati, struttura eteroschedastica compresa.

**L'ottimo.** Alla settima settimana la resa è massimizzata da `I = 470`, `D = 21.1`,
`P = 31.4`, con un intervallo di previsione bayesiano al 90% per la OD pari a
**[11.12, 14.9]**. L'ottimizzazione restituisce una distribuzione a posteriori
sull'ottimo, quindi la raccomandazione viene con la sua incertezza invece che come un
singolo punto.

**Limiti.** L'ottimo cade **sulla frontiera** del dominio sperimentale (`I = 470`):
un esperimento che spinga più in alto l'intensità luminosa potrebbe migliorare ancora
la resa, e finché non si fa, quel valore va letto come "almeno 470". L'effetto casuale
è solo sul parametro di velocità, perché metterlo su tutti faceva sovrastimare la
varianza e dava problemi di convergenza. La griglia di ottimizzazione ha 30 punti per
lato, quindi c'è un errore di discretizzazione. Infine le previsioni oltre la sesta o
settima settimana escono dall'intervallo osservato e poggiano interamente sull'assunto
parametrico: la Gompertz non modella la decrescita e per costruzione converge
all'asintoto.

## Cosa non è pubblicato

| Fuori dal repository | Perché |
| --- | --- |
| `Codice/nostoc.xlsx` e `dati_modificati.Rdata` | i dati sperimentali restano in locale |
| `Codice/Stan/fits/` e le catene MCMC | pesanti e rigenerabili dagli script |
| `Paper/` | articoli di altri autori, citati nella relazione |
| i sorgenti LaTeX della relazione | della relazione si pubblica il PDF |

Chi clona il repository trova il codice, i grafici e la relazione, ma per rieseguire
l'analisi servono i dati dell'esperimento.
