'''Segundo taller estadística 2
  John Flórez
  ''' 
.rs.restartR() #Reinicia el entorno de R después de algunos errores de configuración
#Librerías necesarias:
library(tidyverse)
library(MASS)        # negative binomial
library(pscl)        # pseudo-R2
library(caret)       # validación cruzada
library(performance) # chequeos de modelo
install.packages("DHARMa")
library(DHARMa)
install.packages("see")
library(see)
install.packages("ggplot2")
library("ggplot2")

rm(list = ls()) # Para limpiar el environment en Rstudio

## Cargamos los datos
data_env <- read.csv(file.choose(), check.names = FALSE)
data_sp <- read.csv(file.choose(), check.names = FALSE)

##Unimos los datos:
data <- data_env %>%
  inner_join(data_sp, by = "Plot_name")

##Especies con más de 49 individuos:
## Y presentes en más de 9 sitios
species_cols <- setdiff(colnames(data_sp), "Plot_name")

species_summary <- data_sp %>%
  pivot_longer(-Plot_name, names_to = "species", values_to = "abundance") %>%
  group_by(species) %>%
  summarise(
    total_ind = sum(abundance, na.rm = TRUE),
    n_sites = sum(abundance > 0, na.rm = TRUE)
  ) %>%
  filter(total_ind >= 50, n_sites >= 10)

species_summary
##Seleccionamos 5 de las anteriores especies
selected_species <- species_summary$species[1:5]
selected_species

### ESpecie automatizable:
sp <- selected_species[1]

data_sp1 <- data %>%
  dplyr::select(
    Plot_name,
    all_of(sp),
    Elevation = `Elevation (m)`,
    Precip = `Annual precipitation (mm)`,
    DryDef = `Dry Season Moisture Deficit (mm)`,
    pH = `pH (water)`,
    Clay = `Clay (%)`,
    TotalN = `Total N (%)`
  ) %>%
  mutate(
    presence = ifelse(.data[[sp]] > 0, 1, 0),
    abundance = .data[[sp]]
  )

#Modelo GLM binomial
m_bin_full <- glm(
  presence ~ Elevation + Precip + DryDef + pH + Clay + TotalN,
  family = binomial,
  data = data_sp1
)
summary(m_bin_full)

#MOdelo por AIC
m_bin_sel <- step(m_bin_full, direction = "both", trace = FALSE)
summary(m_bin_sel)
AIC(m_bin_sel)

#Evaluación del ajuste:
pR2(m_bin_sel)        
check_model(m_bin_sel)

#Se repite para las 5 especies
results <- list()

for (sp in selected_species) {
  
  df <- data %>%
    dplyr::select(
      Plot_name,
      all_of(sp),
      Elevation = `Elevation (m)`,
      Precip = `Annual precipitation (mm)`,
      DryDef = `Dry Season Moisture Deficit (mm)`,
      pH = `pH (water)`,
      Clay = `Clay (%)`,
      TotalN = `Total N (%)`
    ) %>%
    mutate(
      presence = ifelse(.data[[sp]] > 0, 1, 0),
      abundance = .data[[sp]]
    )
  
  m_bin <- step(
    glm(presence ~ ., family = binomial,
        data = df %>% dplyr::select(-Plot_name, -abundance)),
    trace = FALSE
  )
  
  m_nb <- step(
    glm.nb(abundance ~ .,
           data = df %>% dplyr::select(-Plot_name, -presence)),
    trace = FALSE
  )
  
  results[[sp]] <- list(binomial = m_bin, negbin = m_nb)
}

##Comparación presencia y abundancia

lapply(results, function(x) {
  list(
    presence_vars = names(coef(x$binomial))[-1],
    abundance_vars = names(coef(x$negbin))[-1]
  )
})
## Presencia en abundancia:
vars_bin <- names(coef(results[[1]]$binomial))[-1]

glm.nb(
  abundance ~ .,
  data = data_sp1[, c("abundance", vars_bin), drop = FALSE]
) %>% summary()
## Abundancia en presencia:

data_sp1 <- data_sp1 %>%
  mutate(Clay = as.numeric(as.character(Clay)))

results[[1]]$negbin <- glm.nb(
  abundance ~ Elevation + Precip + DryDef + pH + Clay + TotalN,
  data = data_sp1
)

vars_nb <- all.vars(formula(results[[1]]$negbin))[-1]

glm(
  presence ~ Elevation + Clay,
  family = binomial,
  data = data_sp1
)



summary(m_bin)
