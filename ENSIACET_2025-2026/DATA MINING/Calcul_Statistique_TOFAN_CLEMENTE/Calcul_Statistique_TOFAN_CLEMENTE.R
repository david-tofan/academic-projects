# =================================================================================
# David TOFAN, Esteban CLEMENTE | GROUPE 4
# Date MAJ 01/06/2026
# INP - ENSIACET | Calculs Statistiques et incertitudes
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# L'objectif est d'analyser la tendance du marché, d'y appliquer Holt-Winters 
# pour prendre la décision d'acheter, attendre ou vendre en fonction du prix estimé. 
# Les paramètres sont estimés par des incertitudes sur les troix varibles aléatoires
#  : tendance, résidu et saisonnalité
#*Note* : Certaines portions du code R ont été développées 
# avec l'assistance de l'IA, puis vérifiées 
# et adaptées par les auteurs.
# ==================================================================================


#-----------------------------------------------------------------------------------


# ==================================================================
# ÉTAPE 1 — LIBRAIRIES ET LECTURE DU FICHIER TXT (COURS DE L'ACTION)
# ==================================================================

library(forecast)
library(mvtnorm)

#conversion pour chemin relatif du fichier txt
script_path <- normalizePath(sys.frames()[[1]]$ofile, mustWork = FALSE)

#Declaration des variables : à modifier par l'utilisateur
FICHIER_TXT   <- "donnees_air_liq.txt"
RATIO_TRAIN   <- 0.75
PERIODE       <- 5
N_SIMULATIONS <- 2000
NIVEAU_INVEST <- 0.6



#lecture du fichier
donnees_brutes <- read.table(FICHIER_TXT, sep = ";", dec = ",", header = TRUE,
                             stringsAsFactors = FALSE, fill = TRUE)

#inversion des données extraites du fichier texte
donnees_brutes <- donnees_brutes[nrow(donnees_brutes):1, ] 

#on recupere la donnée Dernier qui correspond au prix a une date donnée du cours
cours <- ts(donnees_brutes$Dernier, frequency = PERIODE)

#On entraine d'abord 75% de notre modèle, on ne recupere que cette partie pour l'instant : validation hors-échantillon
n_train <- floor(RATIO_TRAIN * length(cours))
train <- ts(head(as.numeric(cours), n_train), frequency = PERIODE)
test <- ts(tail(as.numeric(cours), length(cours) - n_train),
              frequency = PERIODE, start = tsp(train)[2] + 1/PERIODE)

cat("Série :", length(cours), "observables | Train :", length(train), "| Test :", length(test), "\n")

# ------------------------------------------------------------------
# FIGURE 1 — Figure train / test du cours de l'action
# ------------------------------------------------------------------

plot(cours,
     main = "Cours Air Liquide (AI.PA) — Historique complet",
     xlab = "Temps (semaines)",
     ylab = "Cours (€)",
     col  = "steelblue",
     lwd  = 1.5)
abline(v   = tsp(test)[1],
       col = "red", lty = 2, lwd = 2)
legend("topleft",
       legend = c("Cours observé", "Limite train / test"),
       col    = c("steelblue", "red"),
       lty    = c(1, 2), lwd = 2)

# ==================================================================
# ÉTAPE 2 — TEST DE NORMALITÉ (Shapiro-Wilk) + Boxplot
# ==================================================================

#provient de la librairie forecast
modele_hw <- hw(train, seasonal = "additive", h = length(test))

decomp   <- decompose(train, type = "additive")

#T-t
tendance <- na.omit(as.numeric(decomp$trend))

#S-t
saison   <- as.numeric(decomp$seasonal)

#R_t
residu   <- na.omit(as.numeric(decomp$random))

tester_normalite <- function(x, nom_variable) {
  cat("\n=====", nom_variable, "=====\n")
  tes <- shapiro.test(x)
  cat("p-value Shapiro-Wilk :", round(tes$p.value, 4), "\n")
  if (tes$p.value < 0.05) cat("→ Normalité REJETÉE à 95%\n") else
    cat("→ Normalité NON rejetée à 95%\n")
  
# ------------------------------------------------------------------
# FIGURE 2, 3, 4 — Figures Shapiro + Error bar
# ------------------------------------------------------------------
  # 2 graphiques sur la même feuille : PDF | Boxplot
  par(mfrow = c(1, 2))
  
  # PDF observée vs normale théorique
  plot(density(x), main = paste("PDF -", nom_variable), xlab = "x")
  lines(density(rnorm(100000, mean(x), sd(x))), lty = 2, col = "red")
  legend("topright", legend = c("Observé", "N(mu,sigma²)"),
         lty = c(1,2), col = c("black","red"), cex = 0.7)
  
  # Boîte à moustaches
  boxplot(x,
          main   = paste("Boxplot -", nom_variable),
          ylab   = "Valeur",
          col    = "lightblue",
          border = "steelblue",
          notch  = FALSE)
  # Ajout de la moyenne (triangle rouge) en plus de la médiane
  points(1, mean(x), pch = 17, col = "red", cex = 1.5)
  legend("topright", legend = c("Médiane", "Moyenne"),
         pch = c(NA, 17), lty = c(1, NA),
         col = c("black", "red"), cex = 0.7)
  
  par(mfrow = c(1, 1))
}

tester_normalite(tendance, "Tendance T_t")
tester_normalite(saison,   "Saisonnalité S_t")
tester_normalite(residu,   "Résidu R_t")

# ==================================================================
# ÉTAPE 3 — LES 6 FONCTIONS [d,p,q,r,e,v]Y()
# ==================================================================
n_min       <- min(length(tendance), length(saison), length(residu))
composantes <- cbind(tendance[1:n_min], saison[1:n_min], residu[1:n_min])
colnames(composantes) <- c("T", "S", "R")

mu_vecteur <- colMeans(composantes)
Sigma      <- cov(composantes)

cat("\nVecteur mu (T, S, R) :", round(mu_vecteur, 4), "\n")
cat("Matrice Sigma :\n"); print(round(Sigma, 4))

# les 6 fonctions
rY <- function(n = N_SIMULATIONS) {
  sims <- rmvnorm(n, mean = mu_vecteur, sigma = Sigma)
  rowSums(sims)
}
dY <- function(n = N_SIMULATIONS)    density(rY(n), n = 512)
pY <- function(y, n = N_SIMULATIONS) ecdf(rY(n))(y)
qY <- function(p, n = N_SIMULATIONS) quantile(rY(n), probs = p)
eY <- function(n = N_SIMULATIONS)    mean(rY(n))
vY <- function(n = N_SIMULATIONS)    var(rY(n))

# ==================================================================
# ÉTAPE 4 — PROPAGATION MONTE CARLO + IC 95%
# ==================================================================

# un tirage, toutes les stats calculées dessus
cours_simules <- rY(N_SIMULATIONS)
esperance_Y   <- mean(cours_simules)
variance_Y    <- var(cours_simules)
IC_95         <- quantile(cours_simules, probs = c(0.025, 0.975))

cat("\n--- Cours moyen prédit E[Ŷ] :", round(esperance_Y, 2), "€\n")
cat("--- Variance Var[Ŷ]         :", round(variance_Y, 4), "\n")
cat("--- IC 95%                  :", round(IC_95, 2), "\n")

# Barre d'erreur de l'estimateur de la moyenne (loi de Student)
mu_Y  <- esperance_Y
sd_Y  <- sqrt(variance_Y)
IC_mu <- mu_Y + qt(c(0.025, 0.975), df = N_SIMULATIONS - 1) * sd_Y / sqrt(N_SIMULATIONS)
cat("--- IC 95% de E[Ŷ]         :", round(IC_mu, 4), "€\n")

# ------------------------------------------------------------------
# FIGURE 6 — Distribution Monte Carlo du cours prédit
# ------------------------------------------------------------------
dens <- density(cours_simules, n = 512)
plot(dens,
     main = "Distribution du cours prédit Ŷ (Monte Carlo, N=2000)",
     xlab = "Cours (€)", ylab = "Densité")
abline(v   = IC_95,       lty = 2, col = "red",       lwd = 2)
abline(v   = esperance_Y, lty = 1, col = "darkgreen", lwd = 2)
legend("topright",
       legend = c("IC 95%", "Espérance E[Ŷ]"),
       col    = c("red", "darkgreen"),
       lty    = c(2, 1), lwd = 2)

# ==================================================================
# ÉTAPE 5 — PRISE DE DÉCISION
# ==================================================================
cours_actuel <- as.numeric(tail(cours, 1))
proba_gain   <- 1 - pY(cours_actuel)
proba_perte  <- pY(cours_actuel)
q_bas        <- qY(1 - NIVEAU_INVEST)
q_haut       <- qY(NIVEAU_INVEST)

cat("\n========================================\n")
cat("  ANALYSE DE DÉCISION D'INVESTISSEMENT\n")
cat("========================================\n")
cat("Cours actuel               :", round(cours_actuel, 2), "€\n")
cat("Cours moyen prédit E[Ŷ]   :", round(esperance_Y, 2),  "€\n")
cat("IC 95%                     :", round(IC_95, 2),        "\n")
cat("P(gain | investissement)   :", round(proba_gain  * 100, 1), "%\n")
cat("P(perte | investissement)  :", round(proba_perte * 100, 1), "%\n")
cat("Seuil de décision          :", NIVEAU_INVEST * 100,         "%\n")

if (cours_actuel < q_bas) {
  decision    <- "ACHETER"
  explication <- paste0("Cours (", round(cours_actuel, 2), " €) < Q",
                        round((1-NIVEAU_INVEST)*100), "% (",
                        round(q_bas, 2), " €) → P(gain) = ",
                        round(proba_gain*100, 1), "%")
} else if (cours_actuel > q_haut) {
  decision    <- "VENDRE"
  explication <- paste0("Cours (", round(cours_actuel, 2), " €) > Q",
                        round(NIVEAU_INVEST*100), "% (",
                        round(q_haut, 2), " €) → P(perte) = ",
                        round(proba_perte*100, 1), "%")
} else {
  decision    <- "ATTENDRE"
  explication <- paste0("Cours dans zone neutre [", round(q_bas, 2),
                        " ; ", round(q_haut, 2), " €] → H0 non rejetée")
}

cat("\n--- DÉCISION :", decision, "\n")
cat(explication, "\n")
cat("========================================\n")


# ==================================================================
# ANNEXES - Calculs et affichage du modèle et de ses parametres
# ==================================================================

# Paramètres optimisés
cat("=== Paramètres Holt-Winters ===\n")
cat("alpha :", modele_hw$model$par["alpha"], "\n")
cat("beta  :", modele_hw$model$par["beta"], "\n")
cat("gamma :", modele_hw$model$par["gamma"], "\n")

# ------------------------------------------------------------------
# FIGURE 7 — visualisation du train
# ------------------------------------------------------------------

plot(modele_hw,
     main = "Holt-Winters — Train / Test / Prédiction",
     xlab = "Temps", ylab = "Cours (€)")

# Ajouter le test réel en rouge
lines(test, col = "red", lwd = 2)

# Ajouter la prédiction en bleu
lines(modele_hw$mean, col = "blue", lwd = 2, lty = 2)

legend("topleft",
       legend = c("Train observé", "Test réel", "Prédiction HW"),
       col    = c("black", "red", "blue"),
       lty    = c(1, 1, 2),
       lwd    = 2)
