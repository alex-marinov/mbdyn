/*header file for module-FlightGear-Conti_GSOC.cc*/

struct FGNetFDMWordSet : public HighParser::WordSet {
	virtual bool IsWord(const std::string& s) const;
};

struct FGNetCtrlsWordSet : public HighParser::WordSet {
	virtual bool IsWord(const std::string& s) const;
};
extern FGNetFDMWordSet fGNetFDMWordSet;
extern FGNetCtrlsWordSet fGNetCtrlsWordSet;

struct FieldInfo{
	const std::type_info &type; /*data type of the field (double, float, etc)*/
	int typeSize;
	int position; /*position inside of the data structure*/
	size_t offset; /*byte offset from the beginning of the data structure*/
};
typedef struct FieldInfo FieldInfo;

typedef std::map<std::string, FieldInfo> FieldsDescriptionFG; /*key: field name, mapped type: FieldInfo*/

extern FieldsDescriptionFG fieldsDescriptionFGNetFDM;
extern FieldsDescriptionFG fieldsDescriptionFGNetCtrls;
extern int FGNetFDMCurVersion;
extern int FGNetCtrlsCurVersion;

extern "C" int module_init(const char *module_name, void *pdm, void *php);
void buildFieldsDescriptionFG(void);
void readFGStructuresFromFile(const char *filePath, std::string currentVersionLabel, FieldsDescriptionFG &fieldsDescriptionFG, int &currentVersion);
void printOptionsOnTextFile(const char * fileName, FieldsDescriptionFG &fieldsDescriptionFG);
void buildFGBufCast(std::vector<BufCast *>& data, FieldsDescriptionFG &fieldsDescriptionFGInUse);
BufCast *buildOneFGBufCast(size_t& offset, TypeMap_t& swapmap, const std::type_info &fieldType);

/*---------------------------------------------------------------------------------------------------------------------------------------------*/
/*specific definitions for Flight Gear stream output element-----------------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------------------------*/
extern HighParser::WordSet *sendToFGWordSetInUse;

typedef std::map<std::string, ScalarValue *> FlightGearUserChannels;

/*content-type for sending data to Flight Gear*/
struct FlightGearStreamOutputReader : public StreamOutputContentTypeReader{
	virtual StreamContent* Read(DataManager* pDM, MBDynParser& HP);

	void ReadFlightGearScalarValues(DataManager *pDM, MBDynParser& HP,
		FlightGearUserChannels &flightGearUserChannels, FieldsDescriptionFG *dataStructureInUse);
	void setDefaultValues(std::vector<ScalarValue *> &allValues, std::vector<bool> &assignedByUser, FieldsDescriptionFG &fieldsDescriptionFGInUse);
	bool fieldAlreadyUsed(std::string fieldDescription, FlightGearUserChannels &flightGearUserChannels);

	StreamContent::Modifier *buildFGStreamContentModifier(FieldsDescriptionFG &fieldsDescriptionFGInUse);
};

/*---------------------------------------------------------------------------------------------------------------------------------------------*/
/*specific definitions for Flight Gear file drive and file drive caller------------------------------------------------------------------------*/
/*---------------------------------------------------------------------------------------------------------------------------------------------*/

/*FileDrive type for getting data from Flight Gear*/
struct FlightGearFileDriveReader : public FileDriveContentTypeReader{
	virtual StreamDrive::Modifier * Read(std::vector<doublereal> &v0, MBDynParser& HP, int &idrives);

	StreamDrive::Modifier *buildFGStreamDriveModifier(FieldsDescriptionFG *fieldsDescriptionFGInUse);
};

/*FileDriveCaller type for getting data from Flight Gear*/
struct FlightGearFileDriveCallerTypeReader : public FileDriveCallerTypeReader{
	virtual integer Read(const DataManager* pDM, MBDynParser& HP, FileDrive* pDrv);
};

/*--------------modifiers for automatically gettinh which FG data structure is read by the specified FileDrive--------------------*/
class FGNetFDMStreamDriveCopyCast : public StreamDriveCopyCast{
	public:
		FGNetFDMStreamDriveCopyCast(size_t size, const std::vector<BufCast *>& data);
		~FGNetFDMStreamDriveCopyCast(void);
};

class FGNetCtrlsStreamDriveCopyCast : public StreamDriveCopyCast{
	public:
		FGNetCtrlsStreamDriveCopyCast(size_t size, const std::vector<BufCast *>& data);
		~FGNetCtrlsStreamDriveCopyCast(void);
};
/*--------------------------------------------------------------------------------------------------------------------------------*/
