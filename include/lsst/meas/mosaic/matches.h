// -*- lsst-c++ -*-
#include "lsst/afw/table.h"

namespace lsst {
  namespace meas {
    namespace mosaic {

      afw::table::Schema copySchema(afw::table::Schema & schema,
				    afw::table::Schema & target,
				    std::string targetPrefix="",
				    std::string sourcePrefix="");

      template<typename CatT>
      CatT copyCatalog(afw::table::BaseCatalog & catalog,
					    CatT & target,
					    //afw::table::Schema & sourceSchema=NULL,
					    std::string targetPrefix="",
					    std::string sourcePrefix="");

      afw::table::ReferenceMatchVector  matchesFromCatalog(afw::table::BaseCatalog & catalog);

    }
  }
}
