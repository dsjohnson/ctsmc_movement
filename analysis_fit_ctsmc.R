
### Load necessary packages
source("load_packages.R")

### Setup parallel processing
registerDoFuture()
plan(multisession, workers=6)

### Source 'helper' code
source("helper_code/fit_ctsmc.R")
source("helper_code/hexify.R")
source("helper_code/vector2scalar.R")
source("helper_code/make_disc_path.R")
source("helper_code/make_glm_data.R")
source("helper_code/pmp.R")

### Load pup data
load("data_products/pup_frame.RData")

### create output directories and lists
if(!dir.exists("results")) dir.create("results")


### Convert environmental data to hex grids and save the file names to the pup data
pup_frame %>% mutate(
  hex_df = pmap(.,~future(save_hex_df(...))) %>% values
) -> pup_frame


### Fit models
pup_frame %>% mutate(
  fit = pmap(.,~fit_ctsmc(...))
) -> pup_frame
pup_frame %>% unnest(fit) -> pup_frame

saveRDS(pup_frame, "results/pup_frame_fitted.rds")

plan('sequential')
