#include"configs.h"

namespace Configs{

const int_t CheckPointMagic=0x53484c486d687045;
const int_t CheckPointVersion=4;
const int_t DataFormatVersion=2;

const int ExportHeaderCount=4;
const char *SaveNameCheckpoint="checkpoint.dat";
const char *SaveNameReadme="readme.txt";
const char *SaveNameIndex="index.dat";
const char *SaveNameDirectory="data/";
const char *SaveRotationalDataExtension=".rot";
const char *SaveOrbitalDataExtension=".pos";
const char *SaveSubstepDataExtension="x";
const char *SaveNameInitialDirectory="system_initial";
const char *VersionString="v0.4.x.dev";
const char AuthorName[]={104, 105, 109, 196, 171, 197, 155, 196, 129, 0};

//deprecated
const char *SaveNameStructure="structure.json";
const char *SaveNameTimestamps="timestamps.dat";
const char *SaveDataExtension=".dat";
}
