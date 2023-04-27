

#### LOAD LIBS ####
library(grid)
library(gridExtra)
library(raster)  # call this before dplyr so it doesn't mask select
library(dplyr)
library(stringr)
library(ggplot2)
library(gaiah)
library(forcats)
library(tikzDevice)

dir.create("outputs")

#### CHOOSE WHETHER TO USE PRE-STORED VALUES OR RECOMPUTE THINGS####
COMPUTE_ISO_POSTERIORS <- TRUE  # if false then it will just load these up from a cache
REMAKE_ALL_SUPP_MAPS <- TRUE
RECOMPUTE_PMGCD_GRID <- TRUE


#### SOME HOUSEKEEPING VARIABLES ####
# these are cleaner labels for the regions to be used in forcats functions:
clean_region_names <- c("Alaska to Alberta" = "AK.EastBC.AB",
                        "Pacific Northwest" = "Wa.To.NorCalCoast",
                        "Coastal California" = "CentCalCoast",
                        "Sierra" = "CalSierra",
                        "Rocky Mountain" = "Basin.Rockies",
                        "Eastern" = "Eastern")

# these are the colors for the different regions ordered as above
region_colors <- c(
  rgb(152, 78, 163, maxColorValue=255), # purple
  rgb(77, 175, 74, maxColorValue=255), # green
  rgb(255, 255, 51, maxColorValue=255),  # yellow
  rgb(247, 129, 191, maxColorValue=255), # pink
  rgb(255, 127, 0, maxColorValue=255),  # orange
  rgb(228, 26, 28, maxColorValue=255) # red
)



#### DO OR GET THE ISOTOPE POSTERIORS AND CLIP TO BREEDING RANGE ####
# (Just a few seconds if loading in from cache.  Minutes to hours(?) to recompute?)
ref_birds <- breeding_wiwa_isotopes

if(COMPUTE_ISO_POSTERIORS == TRUE) {
  isotope_ref_bird_results <- isotope_posterior_probs(isoscape = isomap_job54152_prediction,
                                                      ref_birds = ref_birds,
                                                      isoscape_pred_column = "predkrig",
                                                      isoscape_sd_column = "stdkrig",
                                                      self_assign = TRUE)

  save(isotope_ref_bird_results, file = "outputs/isotope_ref_bird_results.rda", compress = "xz")
} else {
  load("outputs/isotope_ref_bird_results.rda")
}

# clip the isotope posteriors to the known breeding range and normalize again
# I should rewrite the function to do this if we have the breeding range...maybe...
# also, at this point I am going to get a simple list of rasters
Miso <- lapply(isotope_ref_bird_results$loo_results, function(y) {
  x <- y$posterior_probs
  x * wiwa_breed / raster::cellStats(x * wiwa_breed, sum)
})





#### CREATE THE SUPPLEMENT TABLE OF ISOTOPE DATA SUMMARIES ####
# this is the table that Kristina recommends.
iso_ref_birds <- extract_isopredictions(isoscape = isomap_job54152_prediction,
                              birds = breeding_wiwa_isotopes,
                              pred = "predkrig",
                              sd = "stdkrig")

# now we summarise that by location
iso_supp_tab <- iso_ref_birds %>%
  group_by(Location) %>%
  summarise(state_or_province = State[1],
            country = Country[1],
            num_samples = n(),
            lat = mean(lat),
            long = mean(long),
            dHp = mean(iso_pred),
            sd_dHp = mean(iso_sd),
            mean_dHf = mean(Isotope.Value),
            sd_dHf = sd(Isotope.Value)
  )

# write that out to outputs
readr::write_csv(iso_supp_tab, path = "outputs/supp_table_isotope_means.csv")

#### GET AND PROCESS THE GENETICS AND HABITAT ####
Mgen <- genetic_posteriors2rasters(breeding_wiwa_genetic_posteriors, genetic_regions)  # this takes about 25 seconds

Mhab <- wiwa_habitat_unclipped * wiwa_breed  # make sure to clip it with the known breeding range.
Mhab_norm <- Mhab / raster::cellStats(Mhab, stat = sum)  # and for later analyses, treat these as posterior probs

#### FIDDLY BOOKKEEPING: Select Only birds with both genetics and isotopes. Set Names wisely ####
# I want to retain only the birds that appear both in Mgen and Miso, and I
# want to order them so that they are all together from different locations
keepers <- intersect(names(Miso), names(Mgen))

# here we input the genetic regions for those short pops.  I get these from the repunits
# file from the wiwa-popgen repo.  Ultimately I will want to have these with the data from the get-go
spr <- read.table(textConnection("ShortPop  region
                                 wAKDE AK.EastBC.AB
                                 wAKYA AK.EastBC.AB
                                 wAKUG AK.EastBC.AB
                                 wAKJU AK.EastBC.AB
                                 wABCA AK.EastBC.AB
                                 wBCMH AK.EastBC.AB
                                 wWADA Wa.To.NorCalCoast
                                 wORHA Wa.To.NorCalCoast
                                 wORMB Wa.To.NorCalCoast
                                 wCAEU Wa.To.NorCalCoast
                                 wCAHM CentCalCoast
                                 wCABS CentCalCoast
                                 wCASL CentCalCoast
                                 wCATE CalSierra
                                 wCACL CalSierra
                                 wCAHU CalSierra
                                 wMTHM Basin.Rockies
                                 wOREL Basin.Rockies
                                 wCOPP Basin.Rockies
                                 wCOGM Basin.Rockies
                                 eQCCM Eastern
                                 eONHI Eastern
                                 eNBFR Eastern
                                 "), header = TRUE, stringsAsFactors = FALSE)




kbirds <- ref_birds %>%
  filter(ID %in% keepers) %>%  # the names got mangled so I have to do this
  mutate(ShortPop = str_replace_all(Short_Name, "[0-9]", "")) %>%
  left_join(spr) %>%
  arrange(region, ShortPop, Short_Name) %>%
  mutate(clean_region = forcats::fct_recode(region,  # here are a few lines to make cleaner names for the Regions
                                            !!!clean_region_names
  )) %>%
  mutate(clean_region = forcats::fct_relevel(clean_region, names(clean_region_names))) %>%
  mutate(Region = clean_region) %>%
  dplyr::select(ID, region, Region, ShortPop, Short_Name, lat, long, everything())


# Now, I think that it is just going to be better to use the Short_Names for all of these
# (at least it will be easier to think about and check everything).  So do the following:
Miso <- Miso[kbirds$ID]
names(Miso) <- kbirds$Short_Name

Mgen <- Mgen[kbirds$ID]
names(Mgen) <- kbirds$Short_Name

# and, in fact, it will be nice to make these RasterStacks, now that they are in the right
# order
Miso <- raster::stack(Miso)
Mgen <- raster::stack(Mgen)

#### SUMMARIZE THE GENETICS RESULTS (DISTRIBUTION OF POSTERIORS) FOR PAPER ####
# my thought here is to make some bar plots that show the distribution of things

# make a data frame of all this
assy <- breeding_wiwa_genetic_posteriors %>%
  left_join(kbirds %>% select(Short_Name, Region)) %>%
  mutate(ass_to_region = forcats::fct_recode(region,  # here are a few lines to make cleaner names for the Regions
                                            !!!clean_region_names)) %>%
  mutate(ass_to_region = forcats::fct_relevel(ass_to_region, names(clean_region_names)))

# now sort it so that it is organized first by the true Region
# and within each it is sorted descending on the posterior prob of being
# assigned to the correct region.
assy2 <- assy %>%
  group_by(ID) %>%
  mutate(correct_ass_prob = posterior[Region == ass_to_region]) %>%
  arrange(Region, desc(correct_ass_prob))

# now we want to get an ordering for the indivuals
ords <- assy2 %>%
  filter(Region == ass_to_region) %>%
  ungroup() %>%
  arrange(Region, desc(correct_ass_prob)) %>%
  mutate(ind = 1:n()) %>%
  select(ID, ind)

assy3 <- left_join(assy2, ords) %>%
  rename(`Region Assigned` = ass_to_region)



g <- ggplot(assy3, aes(x = ind, y = posterior, fill = `Region Assigned`)) +
  geom_bar(stat = "identity", position = "stack", width = 1.0) +
  scale_fill_manual(values = region_colors) +
  facet_grid(. ~ Region, drop=TRUE, space="free", scales="free") +
  xlab("Individuals ordered by posterior to correct region") +
  theme(axis.text.x=element_blank()) +
  theme(axis.ticks.x=element_blank())


dir.create("outputs/figures", recursive = TRUE)
ggsave(g, filename = "outputs/figures/gsi_sim_q_values.pdf", width = 14, height = 5)



#### MAKE THE 1,1,1 COMBO ####
Combo <- comboize(Mgen, Miso, Mhab_norm, 1, 1, 1)




#### MAKE THE BIRD-MAP FIGURES FOR THE PAPER ####
wmap <- get_wrld_simpl()
dir.create("outputs/figures")
gList <- list();
for(i in c("eNBFR03", "wOREL03")) {
  tmp <- comboize_and_fortify(Mgen[[i]], Miso[[i]], Mhab, iso_beta_levels = 1, hab_beta_levels = 1)
  tmp$bird <- i;
  latlong <- kbirds %>% filter(Short_Name == i)

  gList[[i]] <- ggplot(mapping = aes(x=long, y = lat)) +
    coord_fixed(1.3, xlim = c(-170, -50), ylim = c(33, 70)) +
    geom_polygon(data = wmap, aes(group = group), fill = NA, color = "black", size = .05) +
    geom_raster(data = tmp, mapping = aes(fill = prob), interpolate = TRUE) +
    scale_fill_gradientn(colours = c("#EBEBEB", rainbow(7)), na.value = NA) +
    theme_bw() +
    facet_wrap( ~ beta_vals, ncol = 2) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    geom_point(data = latlong, mapping = aes(x = long, y = lat), shape = 13, size = 2.5)
}

# now we grid those into a and b panels
afig <- grid.arrange(gList[[1]], top = textGrob("(a)", x = unit(0.1, "npc"), just = "left", gp = gpar(fontsize = 20)))
bfig <- grid.arrange(gList[[2]], top = textGrob("(b)", x = unit(0.1, "npc"), just = "left", gp = gpar(fontsize = 20)))

final_fig <- grid.arrange(afig, bfig, ncol = 1)

ggsave(final_fig, filename = "outputs/figures/figure1.pdf", height = 10, width = 14)


#### MAKE ALL THE INDIVIDUAL BIRD-MAP FIGURES FOR THE SUPPLEMENT ####
if(REMAKE_ALL_SUPP_MAPS == TRUE) {
  wmap <- get_wrld_simpl()
  dir.create("outputs/birdmaps")
  for(i in 1:nlayers(Mgen)) {
    tmp <- comboize_and_fortify(Mgen[[i]], Miso[[i]], Mhab, iso_beta_levels = 1, hab_beta_levels = 1)
    tmp$bird <- names(Mgen)[i];

    latlong <- kbirds %>% filter(Short_Name == names(Mgen)[i])


    g <- ggplot(mapping = aes(x=long, y = lat)) +
      coord_fixed(1.3, xlim = c(-170, -50), ylim = c(33, 70)) +
      geom_polygon(data = wmap, aes(group = group), fill = NA, color = "black", size = .05) +
      geom_raster(data = tmp, mapping = aes(fill = prob), interpolate = TRUE) +
      scale_fill_gradientn(colours = c("#EBEBEB", rainbow(7)), na.value = NA) +
      theme_bw() +
      facet_grid(bird ~ beta_vals) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
      geom_point(data = latlong, mapping = aes(x = long, y = lat), shape = 13, size = 2.5)

    ggsave(filename = paste("outputs/birdmaps/bird_", names(Mgen)[i], ".pdf", sep = ""), height = 3, width = 19)

    print(i);
  }

  # later we will latex them altogether.
}



#### COMPUTE THE MEAN POSTERIOR GREAT CIRCLE DISTANCES ####
# First, we make a rasterStack, GCD, that gives the great circle distance between each reference birds true
# location and every cell in the raster.
GCD <- lapply(1:nrow(kbirds), function(i)
  great_circle_raster(wiwa_breed, kbirds$lat[i], kbirds$long[i])
)
names(GCD) <- kbirds$Short_Name
GCD <- raster::stack(GCD)

# compute the posterior mean values here
genpmgc <- raster::cellStats(Mgen * GCD, stat = sum)
names(genpmgc) <- names(Mgen)
isopmgc <- raster::cellStats(Miso * GCD, stat = sum)
combopmgc <- raster::cellStats(Combo * GCD, stat = sum)
habpmgc <- raster::cellStats(GCD * Mhab_norm, stat = sum)
names(habpmgc) <- names(GCD)

# now, make a data frame with all that in there:
wide_pmGCD <- data_frame(Short_Name = names(genpmgc),
                    genetics = genpmgc,
                    isotopes = isopmgc,
                    habitat = habpmgc,
                    combo = combopmgc) %>%
  left_join(kbirds %>% select(Short_Name, Region))

# make one for a scatterplot
pmGCD_tidy <- wide_pmGCD %>%
  tidyr::gather(., key = data_type, value = pmgcd, genetics:habitat)

# make another for boxplots
pmGCD_super <- wide_pmGCD %>%
  tidyr::gather(., key = data_type, value = pmgcd, genetics:combo)


#### MAKE THE POSTERIOR MEAN GCD SCATTERPLOTS ####
# here is everyone on the same figure
tikz("outputs/figures/pmgcd_altogether_figure.tex", width = 10, height = 4)
gplot <- ggplot(pmGCD_tidy, aes(x = combo, y = pmgcd, colour = Region)) +
  geom_point() +
  facet_wrap(~ data_type) +
  geom_abline(slope = 1, intercept = 0) +
  scale_colour_manual(values = region_colors) +
  xlab("$S^\\mathrm{pm}_i$ (km), data sources combined") +
  ylab("$S^\\mathrm{pm}_i$ (km), single data source") +
  theme(legend.position="top")
print(gplot)
dev.off()
# then you have to latex that tikz file. (and later PDFcrop it)
system("cd outputs/figures/; pdflatex pmgcd_altogether.tex;")

#### MAKE THE POSTERIOR MEAN GCD BOXPLOTS ####

tikz("outputs/figures/pmgcd_boxplots_figure.tex", width = 4.2, height = 5)
bplot <- ggplot(pmGCD_super, aes(y = pmgcd, x = data_type, fill = Region)) +
  facet_wrap(~ Region) +
  scale_fill_manual(values = region_colors) +
  geom_boxplot() +
  xlab("Data type used") +
  ylab("{\\large $S^\\mathrm{pm}_i$ (km)}") +
  guides(fill="none") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  theme(axis.title.y=element_text(margin=margin(0,11,0,0)))

print(bplot)
dev.off()
system("cd outputs/figures/; pdflatex pmgcd_boxplots.tex;")




#### DO CALCULATIONS FOR THE ACTUAL CIBOLA MIGRANTS ####
# these are for all the migrants kristen had in the first paper
MigGenAll <- genetic_posteriors2rasters(migrant_wiwa_genetic_posteriors, genetic_regions)  # this takes a minute or so

# these are just the Arizona birds that Kristina Sampled
MigIsoAll <-  isotope_posterior_probs(isoscape = isomap_job54152_prediction,
                                      ref_birds = ref_birds,
                                      assign_birds = migrant_wiwa_isotopes,
                                      isoscape_pred_column = "predkrig",
                                      isoscape_sd_column = "stdkrig",
                                      self_assign = FALSE)

# and from those we need to get just the posterior probs, and clip them by the breeding range
MigIso <- lapply(MigIsoAll$regular, function(y) {
  x <- y$posterior_probs
  x * wiwa_breed / raster::cellStats(x * wiwa_breed, sum)
})

# now, focus just on the Mig birds that we have data for and get everyone in the correct order
# and make stacks of them
MigBirds <- intersect(names(MigIso), names(MigGenAll))
MigGen <- MigGenAll[MigBirds] %>% raster::stack()
MigIso <- MigIso[MigBirds] %>% raster::stack()

# and then comboize them
MigCombo <- comboize(MigGen, MigIso, Mhab_norm, 1, 1, 1)
names(MigCombo) <- names(MigGen)

# identify where they seem to be headed, on the basis of genetics
mig_gen_assignments <- migrant_wiwa_genetic_posteriors %>%
  filter(ID %in% MigBirds) %>%
  group_by(ID) %>%
  arrange(desc(posterior)) %>%
  summarise(ass_reg = first(region, order_by = desc(posterior)),
            maxpost = first(posterior, order_by = desc(posterior)))

#### FOR THE CIBOLA MIGRANTS, COMPUTE THE POSTERIOR MEAN REMAINING MIGRATION DISTANCE ####
# the lat and long are the same for all those birds:
tmp <- migrant_wiwa_isotopes %>% group_by(lat, long) %>% tally()
ciblat <- tmp$lat[1]
ciblong <- tmp$long[1]
MigDistStack <- lapply(1:nlayers(MigCombo), function(x) great_circle_raster(wiwa_breed, ciblat, ciblong)) %>%
  setNames(names(MigCombo)) %>%
  raster::stack()

MigDistMeans <- raster::cellStats(MigCombo * MigDistStack, stat = sum)
names(MigDistMeans) = names(MigCombo)

# make a data frame of that, and add the assignments and the collection dates on there
mig_dist_df <- data_frame(ID = names(MigDistMeans), dist = MigDistMeans) %>%
  left_join(mig_gen_assignments) %>%
  left_join(migrant_wiwa_genetic_posteriors %>% group_by(ID) %>% summarise(collect_date = min(Collection_Date))) %>%
  mutate(collect_date = lubridate::dmy(collect_date),
         year = lubridate::year(collect_date),
         yearday = lubridate::yday(collect_date)) %>%
  mutate(clean_region = forcats::fct_recode(ass_reg,  # here are a few lines to make cleaner names for the Regions
                                            !!!clean_region_names
  )) %>%
  mutate(clean_region = forcats::fct_relevel(clean_region, names(clean_region_names))) %>%
  mutate(`Genetically-Assigned Region` = clean_region)

# save that to outputs
save(mig_dist_df, file = "outputs/mig_dist_df.rda", compress = "xz")

#### MAKE THE PLOTS OF REMAINING MIGRATION DISTANCE OF THE CIBOLA BIRDS ####
rem_all_plot <- ggplot(mig_dist_df, aes(x = yearday, y = dist, colour = `Genetically-Assigned Region`)) +
  facet_wrap(~ year) +
  scale_colour_manual(values = region_colors) +
  geom_point() +
  theme(legend.position="top") +
  guides(colour = guide_legend(nrow=2)) +
  xlab("Days Since January 1") +
  ylab("Posterior mean remaining migration distance (km)")


wa_green <- region_colors[clean_region_names == "Wa.To.NorCalCoast"]
rem_wa_plot <- ggplot(mig_dist_df %>% filter(ass_reg == "Wa.To.NorCalCoast"), aes(x = yearday, y = dist)) +
facet_wrap(~ year) +
geom_point(colour = wa_green) +
geom_smooth() +
xlab("Days Since January 1") +
ylab("Posterior mean remaining migration distance (km)")

top <- grid.arrange(rem_all_plot, top = textGrob("(a)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 20)))
bottom <- grid.arrange(rem_wa_plot, top = textGrob("(b)", x = unit(0, "npc"), just = "left", gp = gpar(fontsize = 20)))
grid.arrange(top, bottom, ncol = 1)
final_fig <- grid.arrange(top, bottom)
ggsave(final_fig, filename = "outputs/figures/rem_mig_dist_fig.pdf", width = 5.5, height = 11)




#### COMPUTE PMGCD OVER DIFFERENT VALUES OF THE BETAS ####
# note that GCD was calculated above, but we will just recompute it here too:
# First, we make a rasterStack, GCD, that gives the great circle distance between each reference birds true
# location and every cell in the raster.


if(RECOMPUTE_PMGCD_GRID == TRUE) {
  GCD <- lapply(1:nrow(kbirds), function(i)
    great_circle_raster(wiwa_breed, kbirds$lat[i], kbirds$long[i])
  )
  names(GCD) <- kbirds$Short_Name
  GCD <- raster::stack(GCD)

  griddy <- c(0, 0.33, 0.66, 1.0, 1.33, 1.66, 2.0, 2.5)
  names(griddy) <- griddy

  out <- lapply(c("0.0" = 0.0, "1.0" = 1.0), function(bgen) {
    lapply(griddy, function(biso) {
      lapply(griddy, function(bhab) {
        print(c(bgen, biso, bhab))
        combo_tmp <- comboize(Mgen, Miso, Mhab_norm, beta_gen = bgen, beta_iso = biso, beta_hab = bhab)
        combopmgc <- raster::cellStats(combo_tmp * GCD, stat = sum)
        names(combopmgc) <- names(Mgen)
        dplyr::data_frame(ID = names(combopmgc), pmgcd = combopmgc)
      }) %>%
        dplyr::bind_rows(.id = "beta_hab")
    }) %>%
      dplyr::bind_rows(.id = "beta_iso")
  }) %>% dplyr::bind_rows(.id = "beta_gen")

  saveRDS(out, file = "outputs/128_pmgcd_vals.rds", compress = "xz")
} else {
  out <- readRDS(file = "outputs/128_pmgcd_vals.rds")
}

# now that we have that, let's get the means for each region
region_mean_pmgcd <- kbirds %>%
  select(region, Short_Name) %>%
  left_join(., out %>% rename(Short_Name = ID)) %>%
  group_by(region, beta_gen, beta_iso, beta_hab) %>%
  summarize(mean_pmgcd = mean(pmgcd),
            trimmed10_mean_pmgcd = mean(pmgcd, trim = 0.10)) %>%
  ungroup() %>% arrange(region, mean_pmgcd)


#### MAKE THE PLOTS OF MEAN PMGCD AS A FUNCTION OF THE BETAS ####

# put better names on it
rmp2 <- region_mean_pmgcd %>%
  mutate(clean_region = forcats::fct_recode(region,  # here are a few lines to make cleaner names for the Regions
                                            !!!clean_region_names
  )) %>%
  mutate(clean_region = forcats::fct_relevel(clean_region, names(clean_region_names))) %>%
  mutate(Region = clean_region)


# Here for untrimmed means
ggplot(rmp2, aes(x = as.numeric(beta_hab), y = mean_pmgcd, colour = beta_iso, linetype = beta_gen)) +
  geom_point() +
  geom_line() +
  facet_wrap(~ Region)

ggsave(filename = "outputs/figures/mean_pmgcd_by_betas.pdf", width = 14, height = 10)

