# Load API R packages
library(httr)
library(jsonlite)

# Request COVID tracking data
covidtracking_us <- GET("https://covidtracking.com/api/v1/us/daily.json")
covidtracking_ny <- GET("https://covidtracking.com/api/v1/states/ny/daily.json")
covidtracking_ca <- GET("https://covidtracking.com/api/v1/states/ca/daily.json")

# Convert to R data frame
dailyUS = fromJSON(rawToChar(covidtracking_us$content))
dailyNY = fromJSON(rawToChar(covidtracking_ny$content))
dailyCA = fromJSON(rawToChar(covidtracking_ca$content))

# Convert data to KDPF input format
y_us <- rbind(dailyUS$positiveIncrease, dailyUS$hospitalizedCurrently, 
              dailyUS$inIcuCurrently, dailyUS$onVentilatorCurrently, dailyUS$death)[,dim(dailyUS)[1]:1]
y_ny <- rbind(dailyNY$positiveIncrease, dailyNY$hospitalizedCurrently,
              dailyNY$inIcuCurrently, dailyNY$onVentilatorCurrently, dailyNY$death)[,dim(dailyNY)[1]:1]
y_ca <- rbind(dailyCA$positiveIncrease, dailyCA$hospitalizedCurrently,
              dailyCA$inIcuCurrently, dailyCA$onVentilatorCurrently, dailyCA$death)[,dim(dailyCA)[1]:1]

out_daily <- list(y_us = y_us, y_ny = y_ny, y_ca = y_ca)
save(out_daily, file = paste0("./manuscript version/import/covidtracking-",Sys.Date(),".rdata"))
