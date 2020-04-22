is.installed <- function(mypkg){ is.element(mypkg, installed.packages()[,1])}
if (!is.installed("ggplot2")){
print("Installing ggplot2:
"); install.packages("ggplot2", repos = "http://cran.us.r-project.org")
}else {
print("INFO: ggplot2 is already installed
")
}
if (!is.installed("RColorBrewer")){
print("Installing RColorBrewer:
"); install.packages("RColorBrewer", repos = "http://cran.us.r-project.org")
}else {
print("INFO: RColorBrewer is already installed
")
}
if (!is.installed("ggsignif")){
print("Installing ggsignif:
"); install.packages("ggsignif", repos = "http://cran.us.r-project.org")
}else {
print("INFO: ggsignif is already installed
")
}
if (!is.installed("gtools")){
print("Installing gtools:
"); install.packages("gtools", repos = "http://cran.us.r-project.org")
}else {
print("INFO: gtools is already installed
")
}
if (!is.installed("PSCBS")){
print("Installing PSCBS:
"); install.packages("PSCBS", repos = "http://cran.us.r-project.org")
}else {
print("INFO: PSCBS is already installed
")
}
if (!is.installed("egg")){
print("Installing egg:
"); install.packages("egg", repos = "http://cran.us.r-project.org")
}else {
print("INFO: egg is already installed
")
}
if (!is.installed("corrplot")){
print("Installing corrplot:
"); install.packages("corrplot", repos = "http://cran.us.r-project.org")
}else {
print("INFO: corrplot is already installed
")
}
if (!is.installed("gplots")){
print("Installing gplots:
"); install.packages("gplots", repos = "http://cran.us.r-project.org")
}else {
print("INFO: gplots is already installed
")
}
if (!is.installed("gridExtra")){
print("Installing gridExtra:
"); install.packages("gridExtra", repos = "http://cran.us.r-project.org")
}else {
print("INFO: gridExtra is already installed
")
}
