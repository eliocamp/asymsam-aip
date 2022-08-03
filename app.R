# app.R
library(plumber)

port <- Sys.getenv('PORT', unset = 8000)

server <- plumb("plumber.R")

server$run(
  host = '0.0.0.0',
  port = as.numeric(port)
)
