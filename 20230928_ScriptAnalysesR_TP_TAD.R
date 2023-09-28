#Script analyses statistiques TP TAD - Master 1 
#Date du cours : 28 septembre 2023



# Chargement des librairies -------------------
library(ape) #pour les arbres phylogénétiques et pcoa
library(vegan)
library(stringr) #str_split, split le nom des taxa


#Dossier de travail ------------------------
setwd("~/Documents/20230928TPmeta")


## Chargement des fichiers ----------------------------
# Charger la table de comptage OTU x echantillon
OTU <- read.table(file = "OTU_sample.tsv", sep="\t")
dim(OTU) #547 55
head(OTU)
colnames(OTU)
OTU_Order <- OTU[,order(colnames(OTU))]

# Charger une table décrivant la taxonomie
TAX <- read.table(file = "TAX_sample.tsv", sep="\t")
head(TAX)

# Importer et visionner l'arbre phylogénétique
phy_tree <- read.tree("OTU_meta_2023_mafft_trimal.fas.treefile")
# Phylogenetic tree with 416 tips and 414 internal nodes.
plot(phy_tree, cex=0.5, type='unrooted')

# Importer les métadonnées
Metadata <- read.table("Metadata.txt", h=T)
row.names(Metadata) <- Metadata$sample
Metadata_Order <- Metadata[order(row.names(Metadata)),]
row.names(Metadata_Order) == colnames(OTU_Order)[-55]




## Nettoyage des données ------------------------------

# Rechercher le témoin et s'il y a des OTUs contaminantes dans le témoin
OTU_Order$Temoin
which(OTU_Order$Temoin>0)
# [1]   2 170
OTU_Order[which(OTU_Order$Temoin>0),]
TAX[which(OTU_Order$Temoin>0),]

par(mar=c(5,8,2,2))
barplot(colSums(OTU_Order),
        ylab= '',
        xlab='Nombre de séquences',
        
        cex.names=0.5,
        las=2,
        horiz = T)

# On filtre les OTUs contaminantes (chloroplastes et archébactéries)

# Les archees (domain)
sum(TAX$domain=="Archaea") # 3 OTUs archées

# Les chloroplasts (order)
sum(TAX$order=="Chloroplast") # 46 OTUs chloroplast

# On supprime ces 49 OTUs des 2 tableaux (otu et tax)
exclude <- c(which(TAX$order=="Chloroplast"), which(TAX$domain=="Archaea"))

OTU_Order_Nettoyee <- OTU_Order[-exclude,] # 498 OTUs
colnames(OTU_Order_Nettoyee)
OTU_Order_Nettoyee <- OTU_Order_Nettoyee[,-55]

TAX_Nettoyee <- TAX[-exclude,] # 498 OTUs


# Raréfaction ------------------------------------
OTU_Order_Nettoyee[1:3,1:3]
colSums(OTU_Order_Nettoyee)
OTU_Order_Nettoyee_Rarefie <- seq(1,nrow(OTU_Order_Nettoyee),1) #vecteur de 1 à 498
for (j in 1:ncol(OTU_Order_Nettoyee)) { #pour j de 1 à 55 => pour chaque échantillon
  Liste_OTU_Order_Nettoyee_PerEch <- "Test"
  for (i in 1:nrow(OTU_Order_Nettoyee)) { # pour i de 1 à 498
    Liste_OTU_Order_Nettoyee_PerEch_Ech1_Temp <- rep(row.names(OTU_Order_Nettoyee)[i],OTU_Order_Nettoyee[i,j]) #row.names => noms des OTUs ; nombre de séquence pour cet OTU i dans cet échantillon j => ex vecteurs de 203 fois le nom de l'OTU
    Liste_OTU_Order_Nettoyee_PerEch <- c(Liste_OTU_Order_Nettoyee_PerEch,Liste_OTU_Order_Nettoyee_PerEch_Ech1_Temp) # vecteur avec "Test" puis 203 fois le nom de l'OTU j
  }
  Liste_OTU_Order_Nettoyee_PerEch #vecteuravec "Test" puis x fois le nom de l'OTU i puis y fois le nom de l'OTU i+1 etc... pour l'échantillon j
  Liste_OTU_Order_Nettoyee_PerEch <- Liste_OTU_Order_Nettoyee_PerEch[-1] #Elmination de "Test"
  Liste_OTU_Order_Nettoyee_PerEch_RarefieMin <- Liste_OTU_Order_Nettoyee_PerEch[sample(1:length(Liste_OTU_Order_Nettoyee_PerEch),min(apply(OTU_Order_Nettoyee,2,sum)),replace=FALSE)] #échantillonnage aléatoire avec sample() pour le nombre minimal de séquences
  
  Liste_OTU_Order_Nettoyee_Couplees <- "Test"
  for (i in 1:nrow(OTU_Order_Nettoyee)) { # pour i de 1 à 498
    Liste_OTU_Order_Nettoyee_Couplees_Temp <- sum(Liste_OTU_Order_Nettoyee_PerEch_RarefieMin==row.names(OTU_Order_Nettoyee)[i]) #row.names => noms des OTUs => pour OTU i, combien de fois on la retrouve pour l'échantillon j
    Liste_OTU_Order_Nettoyee_Couplees <- c(Liste_OTU_Order_Nettoyee_Couplees,Liste_OTU_Order_Nettoyee_Couplees_Temp)
  }
  Liste_OTU_Order_Nettoyee_Couplees # vecteur version numérique avec "Test" puis le nombre de séquences OTU i puis le nombre de séquence OTU i+1 etc... pour l'échantillon j
  Liste_OTU_Order_Nettoyee_Couplees <- as.numeric(Liste_OTU_Order_Nettoyee_Couplees[-1]) #Elmination de "Test" puis mis au format numérique
  
  OTU_Order_Nettoyee_Rarefie <- cbind(OTU_Order_Nettoyee_Rarefie,Liste_OTU_Order_Nettoyee_Couplees) #cbind pour tous les échantillons avec comme 1ère colonne vecteur de 1 à 28780
}
head(OTU_Order_Nettoyee_Rarefie)
colSums(OTU_Order_Nettoyee_Rarefie)
OTU_Order_Nettoyee_Rarefie[1:25,1:2]
row.names(OTU_Order_Nettoyee_Rarefie) <- row.names(OTU_Order_Nettoyee)
OTU_Order_Nettoyee_Rarefie <- OTU_Order_Nettoyee_Rarefie[,-1] #Suprression de la colonne avec numéros de 1 à 498
colnames(OTU_Order_Nettoyee_Rarefie) <- colnames(OTU_Order_Nettoyee)
OTU_Order_Nettoyee_Rarefie[1:25,1:2]
apply(OTU_Order_Nettoyee_Rarefie,2,sum) #260





## Courbe de raréfaction ----------------------
# Creer un vecteur pour figurer les échantillons prélevés dans le sédiment en marron, dans l'eau en bleu
colnames(OTU_Order_Nettoyee_Rarefie) == row.names(Metadata_Order)
Metadata_Order_ColCompartiment <- Metadata_Order$Compartiment
Metadata_Order_ColCompartiment[Metadata_Order_ColCompartiment=="Water_column"] <- "blue"
Metadata_Order_ColCompartiment[Metadata_Order_ColCompartiment=="Sediment"] <- "brown"

par(mfrow=c(1,2))
colSums(OTU_Order_Nettoyee)
rarecurve(t(OTU_Order_Nettoyee),label=FALSE,step=50, cex=0.5,col=Metadata_Order_ColCompartiment)
abline(v=260, lty=2)

colSums(OTU_Order_Nettoyee_Rarefie)
rarecurve(t(OTU_Order_Nettoyee_Rarefie),label=FALSE,step=50, cex=0.5,col=Metadata_Order_ColCompartiment)



## Visualisation taxonomie ----------------------
par(mfrow=c(1,1))
OTU_Order_Nettoyee_Rarefie[1:3,1:3]
TAX_Nettoyee
TAX_Nettoyee$class
TAX_Nettoyee[which(str_detect(row.names(TAX_Nettoyee),"Cluster_10\\b")),3]
TAX_Nettoyee[which(str_detect(row.names(TAX_Nettoyee),str_c(sep="",row.names(OTU_Order_Nettoyee_Rarefie)[9],"\\b"))),3]
ClassName <- TAX_Nettoyee[which(str_detect(row.names(TAX_Nettoyee),str_c(sep="",row.names(OTU_Order_Nettoyee_Rarefie)[1],"\\b"))),3]
for (i in 2:nrow(OTU_Order_Nettoyee_Rarefie)) { #pour j de 1 à 498 => pour chaque cluster
  ClassName_Temp <- TAX_Nettoyee[which(str_detect(row.names(TAX_Nettoyee),str_c(sep="",row.names(OTU_Order_Nettoyee_Rarefie)[i],"\\b"))),3]
  ClassName <- c(ClassName,ClassName_Temp)
  }
  
dim(OTU_Order_Nettoyee_Rarefie) #498 54
OTU_Order_Nettoyee_Rarefie_Class <- apply(OTU_Order_Nettoyee_Rarefie[,1:54],2,tapply,ClassName,sum)
dim(OTU_Order_Nettoyee_Rarefie_Class)

par(mar=c(4,4,10,4),xpd=TRUE)
barplot(OTU_Order_Nettoyee_Rarefie_Class, col=rainbow(36))
legend(3,320, legend=row.names(OTU_Order_Nettoyee_Rarefie_Class),ncol=5,
       fill =rainbow(36),cex=0.7)



# Diversité alpha ---------------------
## Nombre d'OTUs -----------------------------
specnumber(t(OTU_Order_Nettoyee_Rarefie))
CompartimentEau <- which(str_detect(colnames(OTU_Order_Nettoyee_Rarefie),"EAU"))
CompartimentSediment <- which(str_detect(colnames(OTU_Order_Nettoyee_Rarefie),"SED"))
specnumber(t(OTU_Order_Nettoyee_Rarefie[,CompartimentEau]))
specnumber(t(OTU_Order_Nettoyee_Rarefie[,CompartimentSediment]))

boxplot(specnumber(t(OTU_Order_Nettoyee_Rarefie[,CompartimentEau])),
        specnumber(t(OTU_Order_Nettoyee_Rarefie[,CompartimentSediment])),
        col=c("blue","brown"),names = c("Eau","Sediment"),ylab="Nbre d'OTUs")

## Indice de Shannon -----------------------
diversity(t(OTU_Order_Nettoyee_Rarefie),index = "shannon")
diversity(t(OTU_Order_Nettoyee_Rarefie[,CompartimentEau]),index = "shannon")
diversity(t(OTU_Order_Nettoyee_Rarefie[,CompartimentSediment]),index = "shannon")

boxplot(diversity(t(OTU_Order_Nettoyee_Rarefie[,CompartimentEau]),index = "shannon"),
        diversity(t(OTU_Order_Nettoyee_Rarefie[,CompartimentSediment]),index = "shannon"),
        col=c("blue","brown"),names = c("Eau","Sediment"),ylab="Indice de Shannon")


# Diversité beta ---------------------
## Indice Bray-Curtis ----------------------
OTU_Order_Nettoyee_Rarefie_MatriceBray <- vegdist(t(OTU_Order_Nettoyee_Rarefie),method="bray",upper=TRUE)

### Calcul NMDS -----------------------------------
OTU_Order_Nettoyee_Rarefie_NMDS <- metaMDS(t(OTU_Order_Nettoyee_Rarefie),distance="bray")
#*** Best solution repeated 3 times
OTU_Order_Nettoyee_Rarefie_NMDS$stress #0.02269167
stressplot(OTU_Order_Nettoyee_Rarefie_NMDS, main="Shepard plot")

### Plot NMDS ---------------------------------------
plot(OTU_Order_Nettoyee_Rarefie_NMDS, type="n",
     cex.main=0.8,
     main=paste("OTU_Order_Nettoyee_Rarefie_NMDS / Bray-Curtis / NMDS - Stress=",
                round(OTU_Order_Nettoyee_Rarefie_NMDS$stress,3)),xlab=c("NMDS1"),ylab=c("NMDS2"))
points(scores(OTU_Order_Nettoyee_Rarefie_NMDS, display = "sites", choices = c(1,2)),
       cex=1.3, lwd=0.5,pch=21,
       bg=Metadata_Order_ColCompartiment,col=Metadata_Order_ColCompartiment)
legend(-4,3.7, legend=c("Eau","Sediment"),ncol=1,
       pt.bg =unique(Metadata_Order_ColCompartiment), pch=21, cex=1)



## Unifrac -------------------



## PCOA --------------------
OTU_Order_Nettoyee_Rarefie.pcoa <- pcoa(dist(t(OTU_Order_Nettoyee_Rarefie)))
OTU_Order_Nettoyee_Rarefie.pcoa
barplot(1-OTU_Order_Nettoyee_Rarefie.pcoa$values$Cumul_eig)
biplot.pcoa(OTU_Order_Nettoyee_Rarefie.pcoa,t(OTU_Order_Nettoyee_Rarefie),
            xlabs="x")



#http://r.qcbs.ca/fr/workshops/
