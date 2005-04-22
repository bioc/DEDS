.First.lib <- function(libname, pkgname) {
    require("methods")
    library.dynam("DEDS", pkgname, libname)
}
