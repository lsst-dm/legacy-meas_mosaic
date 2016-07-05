#include "lsst/afw/table.h"

namespace lsst { namespace meas { namespace mosaic {

      struct ProcessSchema {
	template <typename T>
	void operator()(afw::table::SchemaItem<T> const & item) const {
	  std::string keyName = item.field.getName();
	  std::string keyNameFix = "";
	  if (!sourcePrefix.empty()) {
	    if (keyName.find(sourcePrefix) == 0)
	      keyNameFix = keyName.substr(sourcePrefix.size());
	  } else {
	    keyNameFix = keyName;
	  }
	  if (!keyNameFix.empty() && existing.find(keyNameFix) == existing.end()) {
	    if (!targetPrefix.empty()) {
	      keyNameFix = targetPrefix + keyNameFix;
	    }
	    std::string fieldDoc = item.field.getDoc();
	    std::string fieldUnits = item.field.getUnits();
	    std::string typeStr = item.field.getTypeString();
	    //std::cout << keyName << "," << keyNameFix << "," << typeStr << "," << fieldDoc << "," << fieldUnits << std::endl;
	    if (typeStr.find("Array") != 0) {
		target->addField<T>(keyNameFix, fieldDoc, fieldUnits);
	    } else {
		int size = item.field.getElementCount();
		target->addField<T>(keyNameFix, fieldDoc, fieldUnits, afw::table::FieldBase<T>(size));
	    }
	  }
	}

	static void apply(afw::table::Schema & schema,
			  afw::table::Schema & target,
			  std::string targetPrefix,
			  std::string sourcePrefix
			  ) {
	  std::set<std::string> existing = target.getNames();
	  ProcessSchema f = { &target, existing, targetPrefix, sourcePrefix };
	  schema.forEach(f);
	}

	afw::table::Schema *target;
	std::set<std::string> existing;
	std::string targetPrefix;
	std::string sourcePrefix;
      };

      afw::table::Schema copySchema(afw::table::Schema & schema,
				    afw::table::Schema & target,
				    std::string targetPrefix,
				    std::string sourcePrefix)
      {
	ProcessSchema::apply(schema, target, targetPrefix, sourcePrefix);
	return target;
      }

      struct ProcessSimpleCatalog {
	template <typename T>
	void operator()(afw::table::SchemaItem<T> const & item) const {
	  std::string keyName = item.field.getName();
	  std::string keyNameFix = "";
	  if (!sourcePrefix.empty()) {
	    if (keyName.find(sourcePrefix) == 0)
	      keyNameFix = keyName.substr(sourcePrefix.size());
	  } else {
	    keyNameFix = keyName;
	  }
	  if (!targetPrefix.empty()) {
	    keyNameFix = targetPrefix + keyNameFix;
	  }
	  if (!keyNameFix.empty()) {
	    afw::table::Key<T> kFrom = item.key;
	    afw::table::Key<T> kTo   = targetSchema.find<T>(keyNameFix).key;
	    for (unsigned int i = 0; i < catalog.size(); i++) {
	      target->get(i)->set(kTo, catalog[i].get(kFrom));
	    }
	  }
	}

	static void apply(afw::table::BaseCatalog & catalog,
			  afw::table::SimpleCatalog & target,
			  std::string targetPrefix,
			  std::string sourcePrefix
			  ) {
	  afw::table::Schema sourceSchema = catalog.getSchema();

	  afw::table::Schema targetSchema = target.getSchema();
	  target.reserve(catalog.size());
	  for (unsigned int i = target.size(); i < catalog.size(); i++) {
	    target.addNew();
	  }
	  ProcessSimpleCatalog f = { &target, targetSchema, catalog, targetPrefix, sourcePrefix };
	  sourceSchema.forEach(f);
	}

	afw::table::SimpleCatalog *target;
	afw::table::Schema targetSchema;
	afw::table::BaseCatalog catalog;
	std::string targetPrefix;
	std::string sourcePrefix;
      };

      afw::table::SimpleCatalog copySimpleCatalog(afw::table::BaseCatalog & catalog,
						  afw::table::SimpleCatalog & target,
						  std::string targetPrefix,
						  std::string sourcePrefix)
      {
	ProcessSimpleCatalog::apply(catalog, target, targetPrefix, sourcePrefix);
	return target;
      }

      struct ProcessSourceCatalog {
	template <typename T>
	void operator()(afw::table::SchemaItem<T> const & item) const {
	  std::string keyName = item.field.getName();
	  std::string keyNameFix = "";
	  if (!sourcePrefix.empty()) {
	    if (keyName.find(sourcePrefix) == 0)
	      keyNameFix = keyName.substr(sourcePrefix.size());
	  } else {
	    keyNameFix = keyName;
	  }
	  if (!targetPrefix.empty()) {
	    keyNameFix = targetPrefix + keyNameFix;
	  }
	  if (!keyNameFix.empty()) {
	    afw::table::Key<T> kFrom = item.key;
	    afw::table::Key<T> kTo   = targetSchema.find<T>(keyNameFix).key;
	    for (unsigned int i = 0; i < catalog.size(); i++) {
	      target->get(i)->set(kTo, catalog[i].get(kFrom));
	    }
	  }
	}

	static void apply(afw::table::BaseCatalog & catalog,
			  afw::table::SourceCatalog & target,
			  std::string targetPrefix,
			  std::string sourcePrefix
			  ) {
	  afw::table::Schema sourceSchema = catalog.getSchema();

	  afw::table::Schema targetSchema = target.getSchema();
	  target.reserve(catalog.size());
	  for (unsigned int i = target.size(); i < catalog.size(); i++) {
	    target.addNew();
	  }
	  ProcessSourceCatalog f = { &target, targetSchema, catalog, targetPrefix, sourcePrefix };
	  sourceSchema.forEach(f);
	}

	afw::table::SourceCatalog *target;
	afw::table::Schema targetSchema;
	afw::table::BaseCatalog catalog;
	std::string targetPrefix;
	std::string sourcePrefix;
      };

      afw::table::SourceCatalog copySourceCatalog(afw::table::BaseCatalog & catalog,
						  afw::table::SourceCatalog & target,
						  std::string targetPrefix,
						  std::string sourcePrefix)
      {
	ProcessSourceCatalog::apply(catalog, target, targetPrefix, sourcePrefix);
	return target;
      }

      afw::table::ReferenceMatchVector matchesFromCatalog(afw::table::BaseCatalog & catalog) {

	afw::table::Schema catSchema = catalog.getSchema();
	afw::table::Schema simpleMinSchema = afw::table::SimpleTable::makeMinimalSchema();
	afw::table::Schema refSchema = copySchema(catSchema, simpleMinSchema, "", "ref_");
	afw::table::SimpleCatalog refCatalog(refSchema);
	copySimpleCatalog(catalog, refCatalog, "", "ref_");

	afw::table::Schema sourceMinSchema = afw::table::SourceTable::makeMinimalSchema();
	afw::table::Schema srcSchema = copySchema(catSchema, sourceMinSchema, "", "src_");
	afw::table::SourceCatalog srcCatalog(srcSchema);
	copySourceCatalog(catalog, srcCatalog, "", "src_");

	afw::table::ReferenceMatchVector matches;
	afw::table::Key<double> distKey = catalog.getSchema().find<double>("distance").key;
	afw::table::BaseCatalog::iterator i = catalog.begin();
	afw::table::SimpleCatalog::iterator j = refCatalog.begin();
	afw::table::SourceCatalog::iterator k = srcCatalog.begin();
	for (; i != catalog.end(); i++, j++, k++) {
	  matches.push_back(afw::table::ReferenceMatch(j, k, i->get(distKey)));
	}

	return matches;
      }

}}} // namespace lsst::meas::mosaic
