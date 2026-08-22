#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# Packages installation ----
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# Install jags
system("sudo apt-get install jags", wait = TRUE)

# Install separately packages relying on github for installation
install.packages("devtools")
install.packages("remotes")

# Install sinkr and matreex
devtools::install_github("marchtaylor/sinkr")
remotes::install_gitlab("lessem/rpackages/matreex", host = "https://forge.inrae.fr")

# install cran packages if needed
packages.in <- c("dplyr", "ggplot2", "targets", "tidyr", "readxl", "cowplot",
                 "data.table", "factoextra", "terra", "ggmcmc", "R2jags", 
                 "betareg", "car", "scales", "MASS", "broom.mixed", "lme4", 
                 "modi", "ggridges", "purrr", "checkmate", "FD", "sf", 
                 "rnaturalearth", "rnaturalearthdata", "sinkr", "egg", "xtable", 
                 "spData", "ggnewscale", "lmerTest", "ggh4x", "future", "clustermq")
for(i in 1:length(packages.in)){
  if(!(packages.in[i] %in% rownames(installed.packages()))){
    install.packages(packages.in[i])
  }
}  

# Clustermq set up
fx = function(x) x * 2
clustermq::Q(fx, x=1:3, n_jobs=1)

# Automatically retrieve credentials from the Onyxia environment
mybucket <- "jbarrere"
key <- Sys.getenv('AWS_ACCESS_KEY_ID')
secret <- Sys.getenv('AWS_SECRET_ACCESS_KEY')
token <- Sys.getenv('AWS_SESSION_TOKEN')
endpoint <- Sys.getenv('AWS_ENDPOINT_URL')

# Initialize the client
library(paws)
s3 <- paws::s3(
  config = list(
    credentials = list(
      creds = list(
        access_key_id = key,
        secret_access_key = secret,
        session_token = token
      )
    ),
    endpoint = endpoint
  )
)

# Check connection
listobjbucket <- s3$list_objects_v2(
  Bucket = mybucket,
  MaxKeys = 5
)
if (!is.null(listobjbucket$Contents)) {
  cat("Connection OK!\nFiles:\n")
  for (obj in listobjbucket$Contents) {
    cat("-", obj$Key, "\n")
  }
}

# Download data_ipmfuture.zip if needed
if (!file.exists("/home/onyxia/work/data_impfuture.zip")) {
  s3$download_file(
    Bucket = mybucket,
    Key = "data_impfuture.zip",
    Filename = "data_impfuture.zip"
  )
}

# Unzip data.zip if needed
if (!dir.exists("/home/onyxia/work/IPM_future/data")) {
  unzip("data_impfuture.zip", exdir = "/home/onyxia/work/IPM_future")
  cat("Done!\n")
}
