#include"physics/CelestialSystem.h"
#include"utils/memio.h"
#include"tests/tests.h"
#include"../resource.h"
#define NOMINMAX
#include<Windows.h>

#define TEST_LIBRARY_DIR "__eptest/"

int publish_resource(int resource_id,const char *filename){
    HINSTANCE hInst=NULL;
    HRSRC hResInfo=FindResource(hInst,MAKEINTRESOURCE(resource_id),RT_HTML);
    if(hResInfo==NULL)return 1;
    HGLOBAL hResData=LoadResource(hInst,hResInfo);
    DWORD dwResSize=SizeofResource(hInst,hResInfo);
    if(hResData==NULL||!dwResSize)return 2;
    LPVOID lpResData=LockResource(hResData);
    if(!lpResData)return 3;
    MFILE *mfile=mopen(lpResData,dwResSize);
    mfile->publish(filename);
    fclose(mfile);
    FreeResource(hResData);
    return 0;
}

int test_prepare(){

#define PUBLISH(IDR,PATH) do{                            \
    int ret=publish_resource(IDR,TEST_LIBRARY_DIR PATH); \
    if(ret)return ret;                                   \
}while(0)
    PUBLISH(IDR_SolarConfig ,"SolarSystem_Config.txt"         );
    PUBLISH(IDR_SolarInitial,"SolarSystem_Initial.txt"        );
    PUBLISH(IDR_SolarExtra  ,"SolarSystem_ExtraParameters.txt");
    PUBLISH(IDR_RingSaturn  ,"Rings/699.txt"                  );
    PUBLISH(IDR_GPMoon      ,"Geopotentials/301.txt"          );
    PUBLISH(IDR_GPEarth     ,"Geopotentials/399.txt"          );
    PUBLISH(IDR_GPMars      ,"Geopotentials/499.txt"          );
    PUBLISH(IDR_GPJupiter   ,"Geopotentials/599.txt"          );
    PUBLISH(IDR_GPSaturn    ,"Geopotentials/699.txt"          );
    PUBLISH(IDR_GPUranus    ,"Geopotentials/799.txt"          );
    PUBLISH(IDR_GPNeptune   ,"Geopotentials/899.txt"          );
#undef PUBLISH
    //msystem ms;
    //ms.load_dir(TEST_LIBRARY_DIR,"SolarSystem_Initial.txt");
    //if(ms.mlist.size())
    //    return 0;
    return 0;
}
