#include"physics/mass.h"
#include"utils/memio.h"
#include"tests/tests.h"
#include"../resource.h"
#include"modules/logger.h"
#define NOMINMAX
#include<Windows.h>

#define TEST_LIBRARY_DIR "__internal/"

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

static msystem ms;
const msystem &get_test_msystem(){
    return ms;
}
msystem get_test_subsystem(std::vector<const char *> sids){
    msystem msdst;
    msdst.copy_params(ms);
    for(const char *ssid:sids){
        int_t mid=ms.get_mid(ssid);
        if(mid<0)continue;
        msdst.push_back(ms[mid]);
    }
    msdst.accel();
    msdst.analyse();
    return msdst;
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
    bool success=ms.load(TEST_LIBRARY_DIR "SolarSystem_Config.txt");
    if(!success)
        return 1;
    
    LogInfo("      Passed, ");
    return 0;
}
