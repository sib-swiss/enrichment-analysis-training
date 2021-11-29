
## R and RStudio

### Previous knowledge / Competencies

As is stated in the course prerequisites at the [announcement web page](https://www.sib.swiss/training/), we expect participants to have previous knowledge in:

* statistics beginner level (T-test, multiple testing methods).
* R beginner level (Rstudio, install a library, matrix manipulation, read files). Test your R skills with the quiz [here](https://docs.google.com/forms/d/e/1FAIpQLSdIyeuabd_ZOWXgI1MWHapmaOMu20L9ESkLDZiWnpmkpujyOg/viewform)

### Technical

This course will be streamed, you are thus required to have your own computer with an internet connection, and with latest version of [R](https://cran.r-project.org/)
and the free version of [RStudio](https://www.rstudio.com/products/rstudio/download/) installed.

You can install the necessary packages using:

```r
install.packages("BiocManager")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pathview")
BiocManager::install("enrichplot")

```

### Data

You can find the data here:
