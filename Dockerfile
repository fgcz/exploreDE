FROM fgcz_base_image:latest

# Pin ggplot2 to version 3.5.2 for brush compatibility
RUN Rscript -e 'pak::pkg_install("ggplot2@3.5.2", ask = FALSE, dependencies = NA)'

# Install glmpca for GLM-PCA ordination analysis
RUN Rscript -e 'pak::pkg_install("glmpca", ask = FALSE, dependencies = NA)'

# Pin enrichplot to 1.28.4 for cex.params compatibility in cnetplot
RUN Rscript -e 'pak::pkg_install("bioc::enrichplot@1.28.4", ask = FALSE, dependencies = NA)'

# exploreDE is a special case where it is developed on gitlab but the app is hosted on github
# RUN Rscript -e "devtools::install_github('fgcz/exploreDE', dependencies = FALSE, upgrade = 'never', force = TRUE, quiet = TRUE)"

# Copy the app files -- this line is different from the other apps because the code is a package on its own
COPY /inst/app/ /srv/shiny-server/

# Shiny Server runs as user 'shiny'
RUN chown -R shiny:shiny /srv/shiny-server

# Listening port
EXPOSE 3838
