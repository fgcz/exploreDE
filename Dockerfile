FROM fgcz_base_image:latest

# Copy the app files -- this line is different from the other apps because the code is a package on its own
COPY /inst/app/ /srv/shiny-server/

# Shiny Server runs as user 'shiny'
RUN chown -R shiny:shiny /srv/shiny-server

# Listening port
EXPOSE 3838
