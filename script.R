library(strucchange) # Chow and Cusum
library(car) # VIF
library(lmtest) # Durbin-Watson, Ljung-Box
library(tseries) # Test d'autocorrélation
library(nortest) # Normalité des résidus


rm(list=ls())
graphics.off()

data <- read.csv('data.csv', sep = ',')


# Premières corrélations

y_col <- "Transport.routier"
correlations <- sapply(names(data), function(col_name) {
  if (col_name != y_col) {
    cor(data[[y_col]], data[[col_name]], use = "complete.obs")
  } else {
    NA
  }
})
correlations <- correlations[!is.na(correlations)]
cor_table <- data.frame(
  Variable_Explicative = names(correlations),
  Coef_Correlation = correlations
)
cor_table_sorted <- cor_table[order(-abs(cor_table$Coef_Correlation)), ]
print(cor_table_sorted)

# Définition des variables explicatives

y <- data$Transport.routier 
x1 <- data$GDP/data$CPI
x2 <- 100*data$Ind_Prix_Transport/data$CPI
x3 <- 100*data$Pdiesel/data$CPI
x4 <- data$Qdieselvu
x5 <- data$Transport.ferroviaire.5.+data$Transport.fluvial.6.+data$Oléoducs
x51 <- data$Transport.ferroviaire.5.
x52 <- data$Transport.fluvial.6.
x53 <- data$Oléoducs

nstd=length(y)

nobs=cbind(1:nstd)
namind=cbind(as.character(nobs))

# Variables rupture
dummy1 <- as.numeric(nobs < 24)  # Variable muette pour avant la première rupture (2007)
dummy2 <- as.numeric(nobs >= 24 & nobs < 32)  # Variable muette pour entre les deux ruptures (2013)
dummy3 <- as.numeric(nobs >= 32) # Variable muette pour après la deuxième rupture



################################
# Transformation des variables #
################################

# Normalisation
ys <- data$Transport.routier
x1s <- scale(x1)
x2s <- scale(x2)
x3s <- scale(x3)
x5s <- scale(x5)

# Orthogonalisation
x1_1 <- x1 * dummy1
x1_2 <- x1 * dummy2
x1_3 <- x1 * dummy3

x2_1 <- x2 * dummy1
x2_2 <- x2 * dummy2
x2_3 <- x2 * dummy3

x3_1 <- x3 * dummy1
x3_2 <- x3 * dummy2
x3_3 <- x3 * dummy3

x5_1 <- x5 * dummy1
x5_2 <- x5 * dummy2
x5_3 <- x5 * dummy3

# Segmentation
segment1 <- which(nobs < 25)   # Avant 2008
segment2 <- which(nobs >= 25 & nobs < 32)  # Entre 2008 et 2015
segment3 <- which(nobs >= 32)  # Après 2015

yseg1 <- y[segment1]
yseg2 <- y[segment2]
yseg3 <- y[segment3]

x1seg1 <- x1[segment1]
x1seg2 <- x1[segment2]
x1seg3 <- x1[segment3]

x2seg1 <- x2[segment1]
x2seg2 <- x2[segment2]
x2seg3 <- x2[segment3]

x3seg1 <- x3[segment1]
x3seg2 <- x3[segment2]
x3seg3 <- x3[segment3]

x4seg1 <- x4[segment1]
x4seg2 <- x4[segment2]
x4seg3 <- x4[segment3]

x5seg1 <- x5[segment1]
x5seg2 <- x5[segment2]
x5seg3 <- x5[segment3]

nseg1 = length(yseg1)
nseg2 = length(yseg2)
nseg3 = length(yseg3)



###################
# Choix du modèle #
###################

# Modèle de régression simple
vec <- c(x1, x2, x5)
X <- matrix(vec, ncol=3)
Y = matrix(y,nstd,1)
model_simple = lm(Y ~ x1 + x2 + x5)

# Modèle de régression simple avec log
Y_log=log(Y)
X_log=log(X)
model_log = lm(formula = Y_log ~ X_log)


# Modèle de régression standardisé
model_std = lm(Y ~ x1s + x2s + x5s)

# Modèle de régression avec variables orthogonales et 1 cassure (2007)
Y_ortho1 = Y
X_ortho1 = matrix( c(x1_1, x1_3,
                     x2_1, x2_3,
                     x5_1, x5_3), ncol=6)

model_ortho1 <- lm(Y ~ x1_1 + x1_3 +
                       x2_1 + x2_3 +
                       x5_1 + x5_3)

# Modèle de régression avec variables orthogonales et 2 cassures
Y_ortho2 = Y
X_ortho2 = matrix( c(x1_1, x1_2, x1_3,
                    x2_1, x2_2, x2_3,
                    x5_1, x5_2, x5_3), ncol=9)

model_ortho2 <- lm(Y ~ x1_1 + x1_2 + x1_3 +
                      x2_1 + x2_2 + x2_3 +
                      x5_1 + x5_2 + x5_3)

# Modèle de régression par segment
model_seg1 <- lm(yseg1 ~ x1seg1 + x2seg1 + x5seg1)
X_seg1 = matrix( c(x1seg1, x2seg1, x5seg1), ncol=3)
Y_seg1 = matrix(yseg1,nseg1,1)
subset_data_seg1 <- data[1:24, ]

model_seg2 <- lm(yseg2 ~ x1seg2 + x2seg2 + x5seg2)
X_seg2 = matrix( c(x1seg2, x2seg2, x5seg2), ncol=3)
Y_seg2 = matrix(yseg2,nseg2,1)
subset_data_seg2 <- data[24:32, ]

model_seg3 <- lm(yseg3 ~ x1seg3 + x2seg3 + x5seg3)
X_seg3 = matrix( c(x1seg3, x2seg3, x3seg3), ncol=3)
Y_seg3 = matrix(yseg3,nseg3,1)
subset_data_seg3 <- data[25:39, ]

# Choix du modèle
model = model_seg3
Xmodel = X_seg3
Ymodel = Y_seg3
datamodel = subset_data_seg3
nmodel = nseg3

k = ncol(Xmodel)
K=k+1


xc = cbind(1,Xmodel) 
bhat = model$coefficients 
yf = xc %*% bhat # Fitted values
res = residuals(model) # Residus
scr = t(res) %*% res # Somme des carres des résidus

recrest = recresid(model) # Residus recurrents
recrest2 = recrest^2
scrr = sum(recrest2) # Somme des residus des carres recurrents


### Affichage des resultats 

summary(model)


#########
# Tests #
#########

### Test de Chow

for(i in 4:14){
print(i)
print(sctest(Ymodel ~ Xmodel, type = "Chow", point = i) )}


### Test Cusum

Wr <- efp(Ymodel ~ Xmodel, type = "Rec-CUSUM")
plot(Wr)

cumrr <- cumsum(recrest)/scrr

# Valeurs seuil de la distribution Cusum

c0 = 0.18915
Kp1=K+1

t2 <- ts(Kp1:nmodel)

smin <-((t2-K)/(nmodel-K))-c0
smax <- ((t2-K)/(nmodel-K))+c0

vec2 <- c(smin, cumrr, smax)
cusum2 <- matrix(vec2, ncol = 3); 
matplot(t2, cusum2, type ="l")


### Multicolinéarité: VIF

vif_values <- vif(model)  # Calcul des VIF
vif_table <- datamodel.frame(
  Variable = names(vif_values),
  VIF = vif_values
)
print("Valeurs du test VIF :")
print(vif_table)



#########################
# Tests sur les erreurs #
#########################

### Autocorrélation: Durbin-Watson

dw_test <- dwtest(model)
print("Test Durbin-Watson (Autocorrélation) :")
print(dw_test)

### Autocorrélation: Ljung-Box

lb_test <- Box.test(res, lag = 4, type = "Ljung-Box", fitdf = length(model$coefficients))
print("Test Ljung-Box :")
print(lb_test)

### Autocorrélation: Breusch-Godfrey

bg <- bgtest(model, order=3)
print("Test Breusch-Godfrey :")
print(bg)

### Autocorrélation: Histogramme des résidus

dev.new()
hist(res)

### Autocorrélation: ACF et PACF

pacf(res)
dev.print(device= jpeg, file="pacf.jpeg", width=600)
acf(res)
dev.print(device= jpeg, file="acf.jpeg", width=600)

### Hétéroscédasticité: White

e2=res*res
xcarre=Xmodel*Xmodel

WAUX=lm(formula = e2 ~ Xmodel+xcarre)
summary(WAUX)

### Hétéroscédasticité: Breusch-Pagan

bp <- bptest(model, ~ fitted(model) + I(fitted(model)^2))
print(bp)

### Hétéroscédasticité: Goldfeld et Quandt

gq <- gqtest(model)
print(gq)

# Normalité: Jarque-Bera

# Normalité: Anderson-Darling

# Normalité: Kolmogorov-Smirnov

ks <- ks.test(res, "pnorm")
print(ks)

# Normalité: graphique des résidus

qqnorm(y=res)

# Bootstrap: prévision