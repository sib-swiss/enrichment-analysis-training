
## R and RStudio

### Previous knowledge / Competencies

As is stated in the course prerequisites on the [announcement web page](https://www.sib.swiss/training/upcoming-training-courses), we expect participants to have previous knowledge in:

* statistics beginner level (T-test, multiple testing methods).
* R beginner level (Rstudio, install a library, data frame manipulation, import data from csv file). Test your R skills with the quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSdIyeuabd_ZOWXgI1MWHapmaOMu20L9ESkLDZiWnpmkpujyOg/viewform)

### Technical

This course will be streamed, you are thus required to have your own computer with an internet connection, and with latest the version of [R](https://cran.r-project.org/)
and the free version of [RStudio](https://www.rstudio.com/products/rstudio/download/) installed.

You can install the necessary packages using:

```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
BiocManager::install("enrichplot")
BiocManager::install("biomaRt")
install.packages("ggplot2")
install.packages("ggrepel")
install.packages("msigdbr")
install.packages("ggnewscale")

install.packages("ggridges")
install.packages("tidyverse")



```

