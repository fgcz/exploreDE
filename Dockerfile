FROM fgcz_base_image:latest

# exploreDE is a special case where it is developed on gitlab but the app is hosted on github
# RUN Rscript -e "devtools::install_github('fgcz/exploreDE', dependencies = FALSE, upgrade = 'never', force = TRUE, quiet = TRUE)"

# Copy the app files -- this line is different from the other apps because the code is a package on its own
COPY /inst/app/ /srv/shiny-server/

# Shiny Server runs as user 'shiny'
RUN chown -R shiny:shiny /srv/shiny-server

# Listening port
EXPOSE 3838
