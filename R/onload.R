'.onAttach' <- function(lib, pkg="mets")
  {    
    desc <- packageDescription(pkg)
    ## packageStartupMessage("Loading '", desc$Package, "' package...\n",
    ##                       "\tVersion: ", desc$Version, "\n",
    ##                       "\tOverview: help(package=", desc$Package, ")");
    packageStartupMessage("Loading '", desc$Package, "' version ",desc$Version,".");

  }
