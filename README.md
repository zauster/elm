elm
===

Exact tests for linear models, as in Schlag and Gossner (2013)

## How to download this package

### Regular users

There are two ways to download the R-package:

#### First option

1. Download the **whole repository** with `Download ZIP` to your right
2. Go to the folder where you saved the .zip file
3. Extract `elm-master.zip`
4. Start R (or R-Studio or whatever you use)
5. Inside R: Go into the newly generated folder `elm-master`, using either the GUI or the `setwd("path\to\folder\elm-master")` command. (In this folder you will find the *whole* repository (the R-Package, the manual and all program code, if you are interested in it)
6. Optionally: Check with `list.files()` if you are in the right folder: You should see a file named `elm_0.1.tar.gz`.
7. Install the package with `install.packages("elm_0.1.tar.gz", type = "source")` and load it afterwards with `library(elm)`
8. `example(elm)` will provide you with a simple example or see the manual for more information on how to use the package.

#### Second option

1. Click on `elm_0.1.tar.gz` in the files list above
2. Click on `View Raw` in the center
3. A small dialog will ask you where you want to save the file. Remember the folder.
4. Start R
5. Inside R: Go to the folder which contains the downloaded file and continue with 6. from above.

To just download the manual, click on `elm.pdf` and `View Raw`.

### Developers

Clone the repository onto your computer, using the `git clone https://github.com/zauster/elm.git` command.
