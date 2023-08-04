#Dealing with massive tables in R

library(RSQLite)
library(dbplyr)
library(dplyr)

#remember your df should not have quotes
my_db <- dbConnect(SQLite(),"db.sqlite") #either connect to the existing db file or creates a new one

dbWriteTable(my_db,"dt","./Distance_Localities&RegularPoints.csv",header = T) #reading files to db without using RAM memory
dbWriteTable(my_db,"intersection","./Intersection_RegularPoints_shapefiles.csv",header = T)
dbWriteTable(my_db,"locality","./Localities_to_join.csv",header = T,overwrite = T)

#you can retrieve any table from db by 'tbl(db,"name_table")' and it exists as the connection only

dt <- tbl(my_db,"dt") %>% #do whatever you want with the syntax known from dplyr
  left_join(tbl(my_db,"intersection"),by = "ID_POINT") %>%
  left_join(tbl(my_db,"locality"),by = "LOCALITY") %>%
  filter(SPECIES != SPECIES_INFERRED & DISTANCE <= 50000) %>% 
  collect() #collect is the only step when you accually read sth into RAM memory

#at each step you can use show_query() function to see corresponding SQL querry