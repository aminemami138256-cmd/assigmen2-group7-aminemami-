##***************************
## Assignment 1 Code
##
## Group 8
##
## 2026-10-27
##
##***************************

rm(list = ls())

## Packages used -------
library(tidyverse)
conflicted::conflicts_prefer(dplyr::filter())
library(viridis)
# + scale_color/fill_viridis_c/d()
theme_set(theme_light())

# Startup ends here

## Comment codes ------
# Coding explanations (#, often after the code, but not exclusively)
# names for data frames (df_name), lists (ls_name), ... (Thanks Jacqueline May)

## Reading in data -----

df_trich_sp = read_csv("../data/trichinella_spiralis_te.csv")
df_trich_sp

df_trich_pap = read_csv("../data/GCA_001447755.1.csv")
df_trich_pap

## Filtering the data -----

ls_non_te_sp = c("Simple_repeat", "Low_complexity", "Unknown", 
              "simple_repeat", "Satellite", "rRNA", "snRNA")

df_actual_te_sp = df_trich_sp |> 
  filter(matching_class %notin% ls_non_te_sp)

#Filtering out all non-TE's from the Trichinella spiralis data frame

ls_retransposons_sp = c("LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", 
                        "SINE/tRNA-Deu-L2", "LTR/Gypsy", "LTR/Pao", "LTR/Ngaro", 
                        "LTR/DIRS")

df_trich_sp_retro = df_actual_te_sp |> 
  filter(matching_class %in% ls_retransposons_sp)

df_trich_sp_retro

#Filtering so that only retrotransposons are present in the Trichinella spiralis
#data frame

ls_non_te_pap = c("Simple_repeat", "Low_complexity", "Unknown", 
                  "simple_repeat", "Satellite", "rRNA", "snRNA")

df_actual_te_pap = df_trich_pap |> 
  filter(matching_class %notin% ls_non_te_pap)

#Filtering out all non-TE's from the Trichinella papuae data frame

ls_retransposons_pap = c("LINE", "LINE/CR1", "LINE/RTE-BovB", "SINE/tRNA-RTE", 
                         "LTR/Gypsy", "LTR/Pao", "LTR/Ngaro")

df_trich_pap_retro = df_actual_te_pap |> 
  filter(matching_class %in% ls_retransposons_pap)

df_trich_pap_retro

#Filtering so that only retrotransposons are present in the Trichinella papuae
#data frame

## Combining the data frames -----

df_trich_sp_final = df_trich_sp_retro |> 
  mutate(species = "spiralis") |> 
  mutate(qry_left = str_remove_all(string = qry_left, pattern = "-")) |> 
  mutate(no_bp_in_complement = 
           str_remove_all(string = no_bp_in_complement, pattern = "-")) |> 
  mutate(in_repeat_end = 
           str_remove_all(string = in_repeat_end, pattern = "-"))

df_trich_pap_final = df_trich_pap_retro |> 
  mutate(species = "papuae") |> 
  mutate(qry_left = str_remove_all(string = qry_left, pattern = "[()]")) |> 
  mutate(no_bp_in_complement = 
           str_remove_all(string = no_bp_in_complement, pattern = "[()]")) |> 
  mutate(in_repeat_end = 
           str_remove_all(string = in_repeat_end, pattern = "[()]"))

#Adding a column indicating which species before combining
#Formatting so that the bind_rows function can run without error.
#Formatting done was removing "()" and "-" from values in those 3 columns

df_trich_combined = bind_rows(df_trich_sp_final, df_trich_pap_final)

#Combined data frames for more analysis options

## Figures -----

df_trich_sp_retro |> 
  ggplot(aes(perc_insert)) +
  geom_histogram() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Spiralis")

df_trich_pap_retro |> 
  ggplot(aes(perc_insert)) +
  geom_histogram() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Papuae")

df_trich_sp_retro |> 
  ggplot(aes(perc_insert)) +
  geom_freqpoly() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Spiralis")

df_trich_pap_retro |> 
  ggplot(aes(perc_insert)) +
  geom_freqpoly() +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Papuae")

df_trich_combined |> 
  ggplot(aes(perc_insert)) +
  geom_freqpoly() +
  facet_grid(~ species) +
  labs(x = "Percent Insertion", y = "Frequency", 
       title = "Percent Insertion for Trichinella Papuae and Spiralis")

## Potential novel code -----
  
summary(df_trich_sp_retro$perc_insert)
summary(df_trich_pap_retro$perc_insert)

df_trich_combined_lm = 
  lm(df_trich_combined$perc_insert~df_trich_combined$species)
summary(df_trich_combined_lm)

t.test(df_trich_pap_retro$perc_insert, df_trich_sp_retro$perc_insert)
