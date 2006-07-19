#ifndef DATAMANFORWARD_H
#define DATAMANFORWARD_H


/* used by maps to compare strings case-insensitive */
struct ltstrcase {
	/* case-insensitive string comparison */
	bool operator()(const std::string& s1, const std::string& s2) const {
		return strcasecmp(s1.c_str(), s2.c_str()) < 0;
	};
};

class DataManagerErrors {
public:
	class ErrGeneric {};
	class ErrAssemblyDiverged {};
	class ErrAssemblyMaxIters {};
	class ErrElemNotAllowedInAssembly {};
	class ErrUnknownElem {};
	class ErrUnknownFunction {};
	class ErrUnknownNode {};
	class ErrMissingNodes {};
	class ErrNeedDataManager {};
};

class DataManager;


#endif /* DATAMANFORWARD_H */
