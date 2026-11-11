# I modelli Stan

Tutti i modelli condividono la stessa media, una curva di crescita di **Gompertz**

```
mu_ij = Asym_i * exp( -b2_i * b3_i^t_ij )
```

dove i tre parametri individuali dipendono dai fattori sperimentali `I`, `D`, `P`
attraverso predittori lineari. Quello che cambia da un modello all'altro è **come si
tratta la varianza** e **su quali parametri si mette un effetto casuale**: la cartella
`Prove/` conserva il percorso che ha portato a `modello_finale.stan`.

| # | File | Cosa aggiunge rispetto al precedente | Esito |
| --- | --- | --- | --- |
| 1 | `Prove/1_base_ar1.stan` | matrice di varianze e covarianze **AR(1)**, effetto casuale sul solo asintoto `Asym` | punto di partenza: cattura la correlazione fra misure ripetute, ma non l'aumento della varianza nel tempo |
| 2 | `Prove/2_re_su_tutti_i_parametri.stan` | effetti casuali anche su `b2` e `b3`, sempre con AR(1) | scartato: con gli effetti casuali su tutti i parametri il modello **sovrastima la varianza** dei dati e dà problemi di convergenza |
| 3 | `Prove/3_eteroschedastico.stan` | abbandona l'AR(1) e modella la **varianza crescente nel tempo** | la strada giusta: l'eteroschedasticità è il fenomeno dominante nei dati, più della correlazione seriale |
| 4 | `Prove/4_horseshoe.stan` | prior **horseshoe** sui coefficienti, con AR(1) | tentativo di selezione automatica dei predittori; non è entrato nel modello finale |
| 5 | `Prove/5_eteroschedastico_quadratico.stan` | varianza eteroschedastica a due parametri e termini quadratici nei predittori | i termini quadratici restano, la parametrizzazione della varianza viene semplificata |
| ★ | `../Stan/modello_finale.stan` | effetto casuale sul solo `b3`, errore **i.i.d. eteroschedastico** `sigma^2 * mu^delta` (power-law), prior parsimoniose | **il modello consegnato**: riproduce la struttura eteroschedastica dei dati e converge senza forzature |

La scelta finale, in una riga: l'AR(1) serviva a spiegare una dipendenza che in realtà
era **varianza non costante**, e una volta modellata quella con la power-law della media
la struttura autoregressiva diventa superflua. L'effetto casuale resta solo dove serve,
sulla velocità di crescita.

## Come si eseguono

`modello.R`, nella cartella superiore, compila `modello_finale.stan` con **cmdstanr**,
salva l'oggetto in `fits/` e la diagnostica in `Figures/`. `ottimizzazione.R` rilegge
quell'oggetto e cerca le condizioni che massimizzano la resa.

I cinque modelli di `Prove/` si rieseguono con `Prove/fit_model.R`, che li compila uno
alla volta, e `Prove/diagnostica.R` ne confronta le catene. Servono solo per ripercorrere
il confronto: per riprodurre i risultati della relazione basta `modello.R`.

## Se un domani il percorso riprende

Due cose renderebbero il confronto più solido di una cronologia. La prima è aggiungere
`log_lik` alle `generated quantities` di ogni modello e confrontarli con **`loo`**: una
tabella di ELPD trasforma "abbiamo provato cinque modelli" in "abbiamo scelto questo
perché". La seconda è spostare in un file incluso con `#include` le parti identiche —
la dichiarazione dei dati, la Gompertz, le quantità generate per i *posterior predictive
check* — dato che fra `1_base_ar1.stan` e `2_re_su_tutti_i_parametri.stan` ci sono già
110 righe uguali su 193.
