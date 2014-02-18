'.onAttach' <- function(lib, pkg="mets")
  {    
    desc <- packageDescription(pkg)
    ## packageStartupMessage("Loading '", desc$Package, "' package...\n",
    ##                       "\tVersion: ", desc$Version, "\n",
    ##                       "\tOverview: help(package=", desc$Package, ")");
    packageStartupMessage(desc$Package, " version ",desc$Version);
  }
