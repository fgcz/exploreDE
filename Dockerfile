FROM fgcz_base_image:latest

# Copy the app files
COPY . /srv/shiny-server/

# Shiny Server runs as user 'shiny'
RUN chown -R shiny:shiny /srv/shiny-server

# Listening port
EXPOSE 3838
