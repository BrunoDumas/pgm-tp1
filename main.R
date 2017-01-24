# Rendu : maxime.gasse@gmail.com

#install.packages("bnlearn")
#source("http://bioconductor.org/biocLite.R")
#biocLite(c("graph", "Rgraphviz"))

library("bnlearn")

# charge la base de données alarm
data("alarm")

# infos sur les colonnes (nos variables)
ncol(alarm)
colnames(alarm)

# infos sur les lignes (nos observations)
nrow(alarm)
rownames(alarm)

# dix premières lignes, 5 premières colonnes
alarm[1:10, 1:5]

# lignes 3, 5, 1, colonnes "ANES", "HIST" et "MINV"
alarm[c(3, 5, 1), c("ANES", "HIST", "MINV")]

# Indépendances conditionnelles
res = ci.test(x = "PAP", y = "SHNT", z = as.character(NULL), data = alarm, test = "mi")
res$statistic
res$p.value

res = ci.test(x = "PAP", y = "SHNT", z = "PMB", data = alarm, test = "mi")
res$statistic
res$p.value

# PAP dépendant de SHNT
# PAP indépendant de SHNT sachant PMB

#Inspectez la relation entre les variables PAP et SHNT:
table(alarm[, "PAP"])
plot(alarm[, "PAP"])
prop.table(table(alarm[, "PAP"]))

table(alarm[, "SHNT"])
plot(alarm[, "SHNT"])
prop.table(table(alarm[, "SHNT"]))

ct = table(alarm[, c("PAP", "SHNT")])
prop.table(ct)
prop.table(ct, margin = 1)
prop.table(ct, margin = 2)


# Tests d'indépendances
dep <- function(x, y, z, seuil){
  res = ci.test(x = x, y = y, z = as.character(z), data = alarm, test = "mi")
  
  # TRUE -> dependants
  # FALSE -> dependants
  return(res$p.value < seuil)
}

# STKV indépendant de HR | ∅;
dep(x = "STKV", y = "HR", z = NULL, seuil=0.05)

# STKV indépendant de HR | CO;
dep(x = "STKV", y = "HR", z = "CO", seuil=0.05)

# HR indépendant de CO | ∅;
dep(x = "HR", y = "CO", z = NULL, seuil=0.05)

# HR indépendant de CO | STKV;
dep(x = "HR", y = "CO", z = "STKV", seuil=0.05)

# CO indépendant de STKV | ∅;
dep(x = "CO", y = "STKV", z = NULL, seuil=0.05)

# CO indépendant de STKV | HR.
dep(x = "CO", y = "STKV", z = "HR", seuil=0.05)

# Les variables STKV, HR et CO forment une v-structure
# STKV -> CO <- HR

# Inspectez la relation entre STKV et HR:
mask = rep(TRUE, nrow(alarm))
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x)")

# Inspectez la relation entre STKV et HR sachant CO (vous pouvez remplacer "HIGH" par "LOW" ou "NORMAL"):LOW
mask = alarm[, "CO"] == "HIGH"
p = prop.table(table(alarm[mask, c("STKV", "HR")]), margin = 1)
plot(p, main="p(y|x,z=HIGH)")


#####################################
# Inférence dans un réseau Bayésien #
#####################################

# Tout d'abord, construisez un réseau Bayésien complet (structure et paramètres) avec les instructions suivantes:

# structure
bn = hc(alarm)
graphviz.plot(bn)

# parametres
bn = bn.fit(bn, data = alarm, method = "bayes")
bn[["CO"]]


# Inférence approchée

# On peut calculer (inférer) n'importe quelle probabilité à partir du réseau bayésien avec la commande cpquery(). Exécutez plusieurs fois les instructions suivantes. Qu'observez-vous?

cpquery(bn, event = (STKV == "HIGH"), evidence = (HR == "LOW"))
cpquery(bn, event = (STKV == "HIGH"), evidence = (HR == "LOW" & CO == "LOW"))
# les valeurs changent car sampling

# Inférence exacte

# Récupérez le fichier includes.R qui contient la fonction exact.dist(), et ajoutez-le à votre projet. Vous pouvez désormais faire de l'inférence exacte comme suit:

source("includes.R")

p = exact.dist(bn, event = c("STKV", "HR", "CO"), evidence = TRUE)

sum(p["HIGH", "LOW", ]) / sum(p[, "LOW", ])
# p(STKV=HIGH | HR=LOW)
sum(p["HIGH", "LOW", "LOW"]) / sum(p[, "LOW", "LOW"])
# p(STKV=HIGH | HR=LOW, CO=LOW)


p = exact.dist(bn, event = c("HYP", "STKV","LVV"), evidence = TRUE)
lvv = exact.dist(bn, event = c("LVV"), evidence = TRUE)


# HYP sachant STKV et LVV
t = prop.table(p, margin = c(2,3))
t[,,"HIGH"]=t[,,"HIGH"]*lvv["HIGH"]
t[,,"LOW"]=t[,,"LOW"]*lvv["LOW"]
t[,,"NORMAL"]=t[,,"NORMAL"]*lvv["NORMAL"]
# Somme par rapport aux dimensions d'une matrice en gardant dim 1 et 2
m = margin.table(t, c(1,2))



# L'algorithme PC
vars = colnames(alarm)
g = empty.graph(vars)
for (x in vars) {
  for (y in setdiff(vars, x)) {
    if(dep(x = x, y = y, z = NULL, seuil=0.01)){
      g = set.edge(g, from = x, to = y)
    }
  }
}
graphviz.plot(g)

for (i in c(1:(length(vars)-1))){
  for (j in c((i+1):length(vars))){
    x = vars[i]
    y = vars[j]
    v = g$nodes[[x]]$nbr
    #for(z in v[,"to"],c(x,y)){
    for(z in setdiff(v,c(x,y))){
      if(!dep(x = x, y = y, z = z, seuil=0.01)){
        g = drop.edge(g, from = x, to = y)
        break
      }
    }
  }
}
graphviz.plot(g)
supprimeVstructure <- function(indice){
  for (i in c(1:(length(vars)-1))){
    for (j in c((i+1):length(vars))){
      x = vars[i]
      y = vars[j]
      v = g$nodes[[x]]$nbr
      v = setdiff(v,c(x,y))
      if(indice>length(v)){
        break
      }
      combinaisons = combn(v,indice)
      for(k in c(1:ncol(combinaisons))){
          z = as.array(combinaisons[,k])
          if(!dep(x = x, y = y, z = z, seuil=0.01)){
            g = drop.edge(g, from = x, to = y)
            break
          }
      }
    }
  }
}
supprimeVstructure(1)
graphviz.plot(g)
supprimeVstructure(2)
graphviz.plot(g)
supprimeVstructure(4)
graphviz.plot(g)
supprimeVstructure(5)
graphviz.plot(g)
supprimeVstructure(6)
graphviz.plot(g)


