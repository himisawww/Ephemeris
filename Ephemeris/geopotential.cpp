#include"ephemeris.h"
#include"utils.h"

//we support at most order 8 harmonics...
static const int_t Max_N=8;
#include"geopotential.impl"

geopotential *geopotential::load(const char* file,fast_real ref_radius_factor){
    std::string linebuf;
    FILE *fin=fopen(file,"r");
    if(!fin){
        fprintf(stderr,"%s : Error opening harmonics model\n",file);
        return nullptr;
    }
    int_t Nz,Nt,N=0;
    fast_mpvec *c_table;
    fast_real J[(Max_N-1)*(Max_N+3)];
    fast_real *C=J+(Max_N-1);
    fast_real *S=C+(Max_N*(Max_N+1)-2)/2;
    geopotential *ret=nullptr;
    while((linebuf=readline(fin)).size()){
        const char *chbuf=linebuf.c_str();
        if(N==0){
            if(2!=sscanf(chbuf,"%lld%lld",&Nz,&Nt)){
                fprintf(stderr,"%s : Error reading max numbers of zonal/tesseral degrees\n",file);
                break;
            }
            if(Nz<1)Nz=1;
            if(Nz>Max_N)Nz=Max_N;
            if(Nt<1)Nt=1;
            if(Nt>Max_N)Nt=Max_N;
            if(Nz<Nt)Nz=Nt;
            N=std::max(Nz,Nt);
            int_t csize=sizeof(fast_mpvec)*Precompute_Table_size(N);
            memset(J,0,sizeof(J));
            ret=(geopotential*)malloc(sizeof(geopotential)+csize);
            ret->Nz=Nz;
            ret->Nt=Nt;
            ret->N=N;
            c_table=(fast_mpvec*)ret->c_table;
            memset(c_table,0,csize);
            continue;
        }
        int_t n=0,m=0;
        fast_real Cnm=0,Snm=0;
        int_t nr=sscanf(chbuf,"%lld%lld%lf%lf",&n,&m,&Cnm,&Snm);
        if(nr<3||n<2||m<0||m>n||nr==3&&m!=0||m==0&&Snm!=0){
            fprintf(stderr,"%s : %s :\n"
                "\tshould in format of \n"
                "\tdegree(n>=2)  order(m==0)     Jn   [0]\n"
                "\t\tor\n"
                "\tdegree(n>=2)  order(1<=m<=n)  Cnm  Snm\n",file,chbuf);
            free(ret);
            ret=nullptr;
            break;
        }
        double rrffactor=std::pow(ref_radius_factor,n);
        Cnm*=rrffactor;
        Snm*=rrffactor;
        if(m==0){
            if(n<=Nz){
                J[n-2]+=Cnm;
            }
        }
        else{
            if(n<=Nt){
                C[n*(n-1)/2+m-2]+=Cnm;
                S[n*(n-1)/2+m-2]+=Snm;
            }
        }
    }
    fclose(fin);
    if(!ret)return nullptr;

    //precompute coefficient tables
    for(int_t n=2;n<=N;++n){
        fast_mpvec *c=c_table+Precompute_Table_size(n-1);
        if(n==2){
            c[0].x=-3*C[0];c[0].y=-3*S[0];c[0].z=-3*J[0];c[1].x=-6*S[1];c[1].y=6*(C[1]-J[0]);c[1].z=12*S[0];c[2].x=-3*C[0];c[2].y=12*S[0];c[2].z=-15*C[1]+(9*J[0])/2;c[3].x=-6*S[1];c[3].y=-9*C[1]+(3*J[0])/2;c[3].z=-3*S[0];c[4].x=-6*(C[1]+J[0]);c[4].y=-6*S[1];c[4].z=12*C[0];c[5].x=15*S[0];c[5].y=15*C[0];c[5].z=30*S[1];c[6].x=(3*(-14*C[1]+J[0]))/2;c[6].y=24*S[1];c[6].z=-3*C[0];c[7].x=12*C[0];c[7].y=-3*S[0];c[7].z=15*C[1]+(9*J[0])/2;c[8].x=24*S[1];c[8].y=(3*(14*C[1]+J[0]))/2;c[8].z=-3*S[0];c[9].x=(3*(6*C[1]+J[0]))/2;c[9].y=-6*S[1];c[9].z=-3*C[0];
        }
        else if(n==3){
            c[0].x=-6*C[2];c[0].y=-6*S[2];c[0].z=-4*J[1];c[1].x=-30*S[3];c[1].y=30*C[3]-10*J[1];c[1].z=30*S[2];c[2].x=(-9*C[2])/2+45*C[4];c[2].y=(81*S[2])/2+45*S[4];c[2].z=-90*C[3]+12*J[1];c[3].x=-30*S[3];c[3].y=(15*(-10*C[3]+J[1]))/2;c[3].z=(-15*(3*S[2]+14*S[4]))/2;c[4].x=(3*(C[2]+30*C[4]))/2;c[4].y=-6*(S[2]+10*S[4]);c[4].z=15*C[3]-(3*J[1])/2;c[5].x=-10*(3*C[3]+J[1]);c[5].y=-30*S[3];c[5].z=30*C[2];c[6].x=45*(S[2]-2*S[4]);c[6].y=45*(C[2]+2*C[4]);c[6].z=180*S[3];c[7].x=(15*(-18*C[3]+J[1]))/2;c[7].y=180*S[3];c[7].z=(-45*(C[2]+14*C[4]))/2;c[8].x=(-15*(S[2]+26*S[4]))/2;c[8].y=(-15*(C[2]+30*C[4]))/2;c[8].z=-30*S[3];c[9].x=(81*C[2])/2-45*C[4];c[9].y=(-9*(S[2]+10*S[4]))/2;c[9].z=90*C[3]+12*J[1];c[10].x=180*S[3];c[10].y=(15*(18*C[3]+J[1]))/2;c[10].z=(-45*(S[2]-14*S[4]))/2;c[11].x=(-9*(C[2]+70*C[4]))/2;c[11].y=(-9*(S[2]-70*S[4]))/2;c[11].z=-3*J[1];c[12].x=(15*(10*C[3]+J[1]))/2;c[12].y=-30*S[3];c[12].z=(-45*C[2])/2+105*C[4];c[13].x=(-15*(S[2]-30*S[4]))/2;c[13].y=(-15*(C[2]-26*C[4]))/2;c[13].z=-30*S[3];c[14].x=-6*C[2]+60*C[4];c[14].y=(3*(S[2]-30*S[4]))/2;c[14].z=(-3*(10*C[3]+J[1]))/2;
        }
        else if(n==4){
            c[0].x=-10*C[5];c[0].y=-10*S[5];c[0].z=-5*J[2];c[1].x=-90*S[6];c[1].y=90*C[6]-15*J[2];c[1].z=60*S[5];c[2].x=(-5*C[5])/2+315*C[7];c[2].y=(205*S[5])/2+315*S[7];c[2].z=-315*C[6]+25*J[2];c[3].x=-75*S[6]+420*S[8];c[3].y=(-15*(46*C[6]+56*C[8]-3*J[2]))/2;c[3].z=-30*(3*S[5]+28*S[7]);c[4].x=(15*(C[5]+42*C[7]))/2;c[4].y=-45*(S[5]+14*S[7]);c[4].z=(15*(84*C[6]+504*C[8]-5*J[2]))/8;c[5].x=15*(S[6]+28*S[8]);c[5].y=(5*(60*C[6]+840*C[8]-3*J[2]))/8;c[5].z=(15*(S[5]+14*S[7]))/2;c[6].x=-15*(6*C[6]+J[2]);c[6].y=-90*S[6];c[6].z=60*C[5];c[7].x=105*(S[5]-6*S[7]);c[7].y=105*(C[5]+6*C[7]);c[7].z=630*S[6];c[8].x=(45*(-22*C[6]+56*C[8]+J[2]))/2;c[8].y=45*(17*S[6]+28*S[8]);c[8].z=-90*(C[5]+28*C[7]);c[9].x=(-105*(S[5]+30*S[7]))/2;c[9].y=(-105*(C[5]+42*C[7]))/2;c[9].z=-315*(S[6]+12*S[8]);c[10].x=(15*(36*C[6]+1176*C[8]-J[2]))/8;c[10].y=-90*(S[6]+28*S[8]);c[10].z=(15*(C[5]+42*C[7]))/2;c[11].x=(205*C[5])/2-315*C[7];c[11].y=(-5*(S[5]+126*S[7]))/2;c[11].z=5*(63*C[6]+5*J[2]);c[12].x=45*(17*S[6]-28*S[8]);c[12].y=(45*(22*C[6]+56*C[8]+J[2]))/2;c[12].z=-90*(S[5]-28*S[7]);c[13].x=(-15*(5*C[5]+378*C[7]))/2;c[13].y=(-75*S[5])/2+2835*S[7];c[13].z=-5670*C[8]-(75*J[2])/4;c[14].x=-15*(5*S[6]+308*S[8]);c[14].y=(-15*(8*C[6]+1288*C[8]+J[2]))/4;c[14].z=15*(S[5]-14*S[7]);c[15].x=(15*(46*C[6]-56*C[8]+3*J[2]))/2;c[15].y=-15*(5*S[6]+28*S[8]);c[15].z=-90*C[5]+840*C[7];c[16].x=(-105*(S[5]-42*S[7]))/2;c[16].y=(-105*(C[5]-30*C[7]))/2;c[16].z=-315*(S[6]-12*S[8]);c[17].x=30*C[6]-(15*(1288*C[8]+J[2]))/4;c[17].y=-75*S[6]+4620*S[8];c[17].z=15*(C[5]+14*C[7]);c[18].x=-45*(C[5]-14*C[7]);c[18].y=(15*(S[5]-42*S[7]))/2;c[18].z=(-15*(84*C[6]-504*C[8]+5*J[2]))/8;c[19].x=-90*(S[6]-28*S[8]);c[19].y=(-15*(36*C[6]-1176*C[8]+J[2]))/8;c[19].z=(15*(S[5]-42*S[7]))/2;c[20].x=(-15*(20*C[6]-280*C[8]+J[2]))/8;c[20].y=15*(S[6]-28*S[8]);c[20].z=(15*(C[5]-14*C[7]))/2;
        }
        else if(n==5){
            c[0].x=-15*C[9];c[0].y=-15*S[9];c[0].z=-6*J[3];c[1].x=-210*S[10];c[1].y=210*C[10]-21*J[3];c[1].z=105*S[9];c[2].x=(15*(C[9]+168*C[11]))/2;c[2].y=(435*S[9])/2+1260*S[11];c[2].z=-840*C[10]+45*J[3];c[3].x=-105*(S[10]-36*S[12]);c[3].y=(-105*(22*C[10]+72*C[12]-J[3]))/2;c[3].z=(-105*(5*S[9]+72*S[11]))/2;c[4].x=(15*(11*C[9]+84*(7*C[11]-30*C[13])))/8;c[4].y=(-15*(101*S[9]+84*(23*S[11]+30*S[13])))/8;c[4].z=840*C[10]+9450*C[12]-(135*J[3])/4;c[5].x=105*(S[10]+36*S[12]);c[5].y=(105*(28*C[10]+504*C[12]-J[3]))/8;c[5].z=(105*(5*S[9]+108*S[11]+792*S[13]))/8;c[6].x=(-15*(C[9]+84*(C[11]+30*C[13])))/8;c[6].y=(45*(S[9]+28*(S[11]+18*S[13])))/4;c[6].z=(-15*(28*C[10]+504*C[12]-J[3]))/8;c[7].x=-21*(10*C[10]+J[3]);c[7].y=-210*S[10];c[7].z=105*C[9];c[8].x=210*(S[9]-12*S[11]);c[8].y=210*(C[9]+12*C[11]);c[8].z=1680*S[10];c[9].x=(105*(-26*C[10]+216*C[12]+J[3]))/2;c[9].y=105*(23*S[10]+108*S[12]);c[9].z=(-105*(5*C[9]+216*C[11]))/2;c[10].x=-210*(S[9]+33*S[11]-90*S[13]);c[10].y=-210*(C[9]+57*C[11]+90*C[13]);c[10].z=-840*(2*S[10]+45*S[12]);c[11].x=(105*(44*C[10]+1656*C[12]-J[3]))/8;c[11].y=-840*(S[10]+36*S[12]);c[11].z=(105*(5*C[9]+324*C[11]+3960*C[13]))/8;c[12].x=(105*(S[9]+60*S[11]+2232*S[13]))/8;c[12].y=(105*(C[9]+84*(C[11]+30*C[13])))/8;c[12].z=105*(S[10]+36*S[12]);c[13].x=(435*C[9])/2-1260*C[11];c[13].y=(15*(S[9]-168*S[11]))/2;c[13].z=840*C[10]+45*J[3];c[14].x=105*(23*S[10]-108*S[12]);c[14].y=(105*(26*C[10]+216*C[12]+J[3]))/2;c[14].z=(-525*S[9])/2+11340*S[11];c[15].x=(-675*(C[9]+84*(C[11]-2*C[13])))/4;c[15].y=(-675*(S[9]-84*(S[11]+2*S[13])))/4;c[15].z=(-135*(840*C[12]+J[3]))/2;c[16].x=-105*(7*S[10]+468*S[12]);c[16].y=(-105*(8*C[10]+2088*C[12]+J[3]))/4;c[16].z=(525*S[9])/4-945*(3*S[11]+110*S[13]);c[17].x=(15*(C[9]+168*(C[11]+60*C[13])))/2;c[17].y=(15*(11*S[9]-84*(7*S[11]+510*S[13])))/8;c[17].z=(-15*(28*C[10]-3*(840*C[12]+J[3])))/8;c[18].x=(105*(22*C[10]-72*C[12]+J[3]))/2;c[18].y=-105*(S[10]+36*S[12]);c[18].z=(105*(-5*C[9]+72*C[11]))/2;c[19].x=-210*(S[9]-57*S[11]+90*S[13]);c[19].y=-210*(C[9]-33*C[11]-90*C[13]);c[19].z=-840*(2*S[10]-45*S[12]);c[20].x=(105*(8*C[10]-2088*C[12]-J[3]))/4;c[20].y=-735*S[10]+49140*S[12];c[20].z=(105*(5*C[9]+108*C[11]-3960*C[13]))/4;c[21].x=(105*(S[9]-12*(S[11]+330*S[13])))/4;c[21].y=(105*(C[9]+12*(C[11]-330*C[13])))/4;c[21].z=210*S[10];c[22].x=(-15*(101*C[9]+84*(-23*C[11]+30*C[13])))/8;c[22].y=(15*(11*S[9]-84*(7*S[11]+30*S[13])))/8;c[22].z=-840*C[10]+9450*C[12]-(135*J[3])/4;c[23].x=-840*(S[10]-36*S[12]);c[23].y=(-105*(44*C[10]-1656*C[12]+J[3]))/8;c[23].z=(105*(5*S[9]-324*S[11]+3960*S[13]))/8;c[24].x=(15*(11*C[9]+84*(7*C[11]-510*C[13])))/8;c[24].y=(15*(S[9]-168*(S[11]-60*S[13])))/2;c[24].z=(15*(28*C[10]+3*(840*C[12]+J[3])))/8;c[25].x=(-105*(28*C[10]-504*C[12]+J[3]))/8;c[25].y=105*(S[10]-36*S[12]);c[25].z=(105*(5*C[9]-108*C[11]+792*C[13]))/8;c[26].x=(105*(S[9]-84*(S[11]-30*S[13])))/8;c[26].y=(105*(C[9]-60*C[11]+2232*C[13]))/8;c[26].z=105*(S[10]-36*S[12]);c[27].x=(45*(C[9]-28*C[11]+504*C[13]))/4;c[27].y=(-15*(S[9]-84*(S[11]-30*S[13])))/8;c[27].z=(15*(28*C[10]-504*C[12]+J[3]))/8;
        }
        else if(n==6){
            c[0].x=-21*C[14];c[0].y=-21*S[14];c[0].z=-7*J[4];c[1].x=-420*S[15];c[1].y=-28*(-15*C[15]+J[4]);c[1].z=168*S[14];c[2].x=(63*(C[14]+120*C[16]))/2;c[2].y=(819*S[14])/2+3780*S[16];c[2].z=-1890*C[15]+(147*J[4])/2;c[3].x=18900*S[17];c[3].y=-105*(30*C[15]+180*C[17]-J[4]);c[3].z=-630*(S[14]+20*S[16]);c[4].x=(315*(C[14]+60*(C[16]-22*C[18])))/8;c[4].y=(-1575*(3*S[14]+76*S[16]+264*S[18]))/8;c[4].z=(105*(240*C[15]+3960*C[17]-7*J[4]))/8;c[5].x=(315*(5*S[15]+216*S[17]-792*S[19]))/4;c[5].y=(105*(75*C[15]+1728*C[17]+2376*C[19]-2*J[4]))/4;c[5].z=315*(S[14]+30*S[16]+396*S[18]);c[6].x=(-105*(C[14]+108*C[16]+3960*C[18]))/8;c[6].y=105*(S[14]+36*(S[16]+22*S[18]));c[6].z=(-35*(270*C[15]+7128*C[17]+61776*C[19]-7*J[4]))/16;c[7].x=(-105*(S[15]+72*(S[17]+33*S[19])))/4;c[7].y=(-35*(42*C[15]+1512*C[17]+33264*C[19]-J[4]))/16;c[7].z=(-105*(S[14]+36*(S[16]+22*S[18])))/8;c[8].x=-28*(15*C[15]+J[4]);c[8].y=-420*S[15];c[8].z=168*C[14];c[9].x=378*(S[14]-20*S[16]);c[9].y=378*(C[14]+20*C[16]);c[9].z=3780*S[15];c[10].x=-105*(30*C[15]-540*C[17]-J[4]);c[10].y=6300*(S[15]+9*S[17]);c[10].z=-630*(C[14]+60*C[16]);c[11].x=-630*(S[14]+35*S[16]-330*S[18]);c[11].y=-630*(C[14]+75*C[16]+330*C[18]);c[11].z=-6300*(S[15]+33*S[17]);c[12].x=(105*(105*C[15]+4320*C[17]-2*(5940*C[19]+J[4])))/4;c[12].y=(-1575*(11*S[15]+504*S[17]+792*S[19]))/4;c[12].z=315*(C[14]+90*(C[16]+22*C[18]));c[13].x=(945*(S[14]+68*S[16]+2904*S[18]))/8;c[13].y=(945*(C[14]+108*C[16]+3960*C[18]))/8;c[13].z=(945*(5*S[15]+264*(S[17]+13*S[19])))/4;c[14].x=(-35*(66*C[15]+4968*C[17]+204336*C[19]-J[4]))/16;c[14].y=210*(S[15]+72*(S[17]+33*S[19]));c[14].z=(-105*(C[14]+108*C[16]+3960*C[18]))/8;c[15].x=(819*C[14])/2-3780*C[16];c[15].y=(63*(S[14]-120*S[16]))/2;c[15].z=(21*(180*C[15]+7*J[4]))/2;c[16].x=6300*(S[15]-9*S[17]);c[16].y=105*(30*C[15]+540*C[17]+J[4]);c[16].z=-630*(S[14]-60*S[16]);c[17].x=(-315*(7*C[14]+660*(C[16]-6*C[18])))/4;c[17].y=(-315*(7*S[14]-660*(S[16]+6*S[18])))/4;c[17].z=(-105*(11880*C[17]+7*J[4]))/4;c[18].x=(-1575*(5*S[15]+360*S[17]-792*S[19]))/2;c[18].y=(-105*(15*C[15]+2*(3240*C[17]+5940*C[19]+J[4])))/2;c[18].z=630*(S[14]-30*(S[16]+66*S[18]));c[19].x=(315*(C[14]+180*(C[16]+66*C[18])))/4;c[19].y=(1575*(S[14]-60*S[16]-5544*S[18]))/8;c[19].z=(-105*(90*C[15]-11880*C[17]-308880*C[19]-7*J[4]))/16;c[20].x=(315*(S[15]+144*S[17]+8712*S[19]))/2;c[20].y=(-105*(6*C[15]-3672*C[17]-223344*C[19]-J[4]))/16;c[20].z=(-315*(S[14]-12*(S[16]+198*S[18])))/8;c[21].x=105*(30*C[15]-180*C[17]+J[4]);c[21].y=-18900*S[17];c[21].z=-630*(C[14]-20*C[16]);c[22].x=-630*(S[14]-75*S[16]+330*S[18]);c[22].y=-630*(C[14]-5*(7*C[16]+66*C[18]));c[22].z=-6300*(S[15]-33*S[17]);c[23].x=(105*(15*C[15]-2*(3240*C[17]-5940*C[19]+J[4])))/2;c[23].y=(1575*(-5*S[15]+360*S[17]+792*S[19]))/2;c[23].z=630*(C[14]+30*(C[16]-66*C[18]));c[24].x=(945*(S[14]-20*(S[16]+286*S[18])))/4;c[24].y=(945*(C[14]+20*(C[16]-286*C[18])))/4;c[24].z=(4725*(S[15]-1144*S[19]))/2;c[25].x=(-105*(30*C[15]-2520*C[17]-356400*C[19]-J[4]))/16;c[25].y=(1575*(S[15]-24*(S[17]+253*S[19])))/4;c[25].z=(-315*(C[14]+60*(C[16]-22*C[18])))/8;c[26].x=(-1575*(3*C[14]-76*C[16]+264*C[18]))/8;c[26].y=(315*(S[14]-60*(S[16]+22*S[18])))/8;c[26].z=(-105*(240*C[15]-3960*C[17]+7*J[4]))/8;c[27].x=(-1575*(11*S[15]-504*S[17]+792*S[19]))/4;c[27].y=(105*(-105*C[15]+4320*C[17]+11880*C[19]-2*J[4]))/4;c[27].z=315*(S[14]-90*(S[16]-22*S[18]));c[28].x=(1575*(C[14]+60*C[16]-5544*C[18]))/8;c[28].y=(315*(S[14]-180*(S[16]-66*S[18])))/4;c[28].z=(105*(90*C[15]+11880*C[17]-308880*C[19]+7*J[4]))/16;c[29].x=(1575*(S[15]+24*(S[17]-253*S[19])))/4;c[29].y=(105*(30*C[15]+2520*C[17]-356400*C[19]+J[4]))/16;c[29].z=(-315*(S[14]-60*(S[16]+22*S[18])))/8;c[30].x=(-105*(75*C[15]+2*(-864*C[17]+1188*C[19]+J[4])))/4;c[30].y=(315*(5*S[15]-72*(3*S[17]+11*S[19])))/4;c[30].z=315*(C[14]-30*C[16]+396*C[18]);c[31].x=(945*(S[14]-108*S[16]+3960*S[18]))/8;c[31].y=(945*(C[14]-68*C[16]+2904*C[18]))/8;c[31].z=(945*(5*S[15]-264*(S[17]-13*S[19])))/4;c[32].x=(105*(6*C[15]+3672*C[17]-223344*C[19]+J[4]))/16;c[32].y=(315*(S[15]-144*S[17]+8712*S[19]))/2;c[32].z=(-315*(C[14]+12*(C[16]-198*C[18])))/8;c[33].x=105*(C[14]-36*C[16]+792*C[18]);c[33].y=(-105*(S[14]-108*S[16]+3960*S[18]))/8;c[33].z=(35*(270*C[15]-7128*C[17]+61776*C[19]+7*J[4]))/16;c[34].x=210*(S[15]-72*(S[17]-33*S[19]));c[34].y=(35*(66*C[15]-4968*C[17]+204336*C[19]+J[4]))/16;c[34].z=(-105*(S[14]-108*S[16]+3960*S[18]))/8;c[35].x=(35*(42*C[15]-1512*C[17]+33264*C[19]+J[4]))/16;c[35].y=(-105*(S[15]-72*(S[17]-33*S[19])))/4;c[35].z=(-105*(C[14]-36*C[16]+792*C[18]))/8;
        }
        else if(n==7){
            c[0].x=-28*C[20];c[0].y=-28*S[20];c[0].z=-8*J[5];c[1].x=-756*S[21];c[1].y=756*C[21]-36*J[5];c[1].z=252*S[20];c[2].x=77*C[20]+9450*C[22];c[2].y=707*S[20]+9450*S[22];c[2].z=-3780*C[21]+112*J[5];c[3].x=252*(2*S[21]+275*S[23]);c[3].y=-63*(118*C[21]+1100*C[23]-3*J[5]);c[3].z=-63*(21*S[20]+550*S[22]);c[4].x=(105*(C[20]+45*(C[22]-132*C[24])))/2;c[4].y=(-105*(29*S[20]+945*S[22]+5940*S[24]))/2;c[4].z=210*(45*C[21]+990*C[23]-J[5]);c[5].x=(315*(13*S[21]+88*(7*S[23]-117*S[25])))/4;c[5].y=(315*(97*C[21]+2816*C[23]+10296*C[25]-2*J[5]))/4;c[5].z=(315*(7*S[20]+275*S[22]+5148*S[24]))/2;c[6].x=(-35*(23*C[20]+54*(57*C[22]+2420*C[24]-8008*C[26])))/16;c[6].y=(35*(247*S[20]+54*(207*S[22]+5588*S[24]+8008*S[26])))/16;c[6].z=(-35*(405*C[21]+14256*C[23]+216216*C[25]-8*J[5]))/4;c[7].x=(-945*(S[21]+88*(S[23]+39*S[25])))/4;c[7].y=(-315*(54*C[21]+2376*C[23]+61776*C[25]-J[5]))/16;c[7].z=(-315*(7*S[20]+66*(5*S[22]+156*(S[24]+10*S[26]))))/16;c[8].x=(35*(C[20]+54*(3*C[22]+220*C[24]+8008*C[26])))/16;c[8].y=(-35*(S[20]+54*(S[22]+44*(S[24]+26*S[26]))))/2;c[8].z=(35*(54*C[21]+2376*C[23]+61776*C[25]-J[5]))/16;c[9].x=-36*(21*C[21]+J[5]);c[9].y=-756*S[21];c[9].z=252*C[20];c[10].x=630*(S[20]-30*S[22]);c[10].y=630*(C[20]+30*C[22]);c[10].z=7560*S[21];c[11].x=-189*(34*C[21]-1100*C[23]-J[5]);c[11].y=756*(19*S[21]+275*S[23]);c[11].z=-189*(7*C[20]+550*C[22]);c[12].x=-1575*(S[20]+36*(S[22]-22*S[24]));c[12].y=-1575*(C[20]+96*C[22]+792*C[24]);c[12].z=-18900*(S[21]+44*S[23]);c[13].x=(315*(123*C[21]+5280*C[23]-2*(25740*C[25]+J[5])))/4;c[13].y=(-2835*(23*S[21]+440*(3*S[23]+13*S[25])))/4;c[13].z=(315*(7*C[20]+825*C[22]+25740*C[24]))/2;c[14].x=(945*(5*S[20]+6*(63*S[22]+2948*S[24]-8008*S[26])))/8;c[14].y=(945*(5*C[20]+678*C[22]+30360*C[24]+48048*C[26]))/8;c[14].z=(2835*(5*S[21]+88*(4*S[23]+91*S[25])))/2;c[15].x=(-315*(78*C[21]+6600*C[23]+308880*C[25]-J[5]))/16;c[15].y=(4725*(S[21]+88*(S[23]+39*S[25])))/2;c[15].z=(-315*(7*C[20]+990*(C[22]+52*(C[24]+14*C[26]))))/16;c[16].x=(-315*(S[20]+102*S[22]+8712*S[24]+391248*S[26]))/16;c[16].y=(-315*(C[20]+54*(3*C[22]+220*C[24]+8008*C[26])))/16;c[16].z=(-945*(S[21]+88*(S[23]+39*S[25])))/4;c[17].x=707*C[20]-9450*C[22];c[17].y=77*S[20]-9450*S[22];c[17].z=28*(135*C[21]+4*J[5]);c[18].x=756*(19*S[21]-275*S[23]);c[18].y=189*(34*C[21]+1100*C[23]+J[5]);c[18].z=-189*(7*S[20]-550*S[22]);c[19].x=-105*(14*C[20]+1485*(C[22]-12*C[24]));c[19].y=-1470*S[20]+155925*(S[22]+12*S[24]);c[19].z=-420*(2970*C[23]+J[5]);c[20].x=(-315*(97*S[21]+7480*S[23]-51480*S[25]))/2;c[20].y=(-315*(13*C[21]+9680*C[23]+51480*C[25]+2*J[5]))/2;c[20].z=315*(7*S[20]-55*(5*S[22]+468*S[24]));c[21].x=(105*(67*C[20]+270*(49*C[22]+3476*C[24]-8008*C[26])))/16;c[21].y=(105*(157*S[20]-270*(39*S[22]+4532*S[24]+8008*S[26])))/16;c[21].z=(-14175*C[21])/4+210*(2970*C[23]+135135*C[25]+J[5]);c[22].x=1890*(S[21]+22*(7*S[23]+468*S[25]));c[22].y=(-945*(10*C[21]-5368*C[23]-391248*C[25]-J[5]))/16;c[22].z=(-945*(7*S[20]-22*(5*S[22]+468*(3*S[24]+70*S[26]))))/16;c[23].x=(-35*(5*C[20]+54*(27*C[22]+44*(85*C[24]+5278*C[26]))))/16;c[23].y=(35*(-23*S[20]+918*S[22]+2376*(97*S[24]+5642*S[26])))/16;c[23].z=(35*(27*C[21]-2376*C[23]-216216*C[25]-J[5]))/4;c[24].x=63*(118*C[21]-1100*C[23]+3*J[5]);c[24].y=252*(2*S[21]-275*S[23]);c[24].z=-63*(21*C[20]-550*C[22]);c[25].x=-1575*(S[20]-96*S[22]+792*S[24]);c[25].y=-1575*(C[20]-36*(C[22]+22*C[24]));c[25].z=-18900*(S[21]-44*S[23]);c[26].x=(315*(13*C[21]-9680*C[23]+51480*C[25]-2*J[5]))/2;c[26].y=(-315*(97*S[21]-440*(17*S[23]+117*S[25])))/2;c[26].z=315*(7*C[20]+55*(5*C[22]-468*C[24]));c[27].x=(4725*(S[20]-30*S[22]-8008*(S[24]-2*S[26])))/4;c[27].y=(4725*(C[20]+30*C[22]-8008*(C[24]+2*C[26])))/4;c[27].z=4725*(3*S[21]-8008*S[25]);c[28].x=(-945*(34*C[21]-3960*C[23]-583440*C[25]-J[5]))/16;c[28].y=(945*(19*S[21]-440*(S[23]+351*S[25])))/4;c[28].z=(-945*(7*C[20]+550*C[22]-17160*(C[24]+70*C[26])))/16;c[29].x=(-945*(S[20]+14*S[22]-88*(97*S[24]+10738*S[26])))/16;c[29].y=(-945*(C[20]+74*C[22]-88*(85*C[24]+11102*C[26])))/16;c[29].z=(-945*(3*S[21]+88*(S[23]-91*S[25])))/4;c[30].x=(-105*(29*C[20]-945*C[22]+5940*C[24]))/2;c[30].y=(105*(S[20]-45*(S[22]+132*S[24])))/2;c[30].z=-210*(45*C[21]-990*C[23]+J[5]);c[31].x=(-2835*(23*S[21]+440*(-3*S[23]+13*S[25])))/4;c[31].y=(-315*(123*C[21]-5280*C[23]-51480*C[25]+2*J[5]))/4;c[31].z=(315*(7*S[20]-825*S[22]+25740*S[24]))/2;c[32].x=(105*(157*C[20]+270*(39*C[22]-4532*C[24]+8008*C[26])))/16;c[32].y=(105*(67*S[20]-270*(49*S[22]-44*(79*S[24]+182*S[26]))))/16;c[32].z=(105*(135*C[21]+8*(2970*C[23]-135135*C[25]+J[5])))/4;c[33].x=(945*(19*S[21]+440*(S[23]-351*S[25])))/4;c[33].y=(945*(34*C[21]+3960*C[23]-583440*C[25]+J[5]))/16;c[33].z=(-945*(7*S[20]-110*(5*S[22]+156*(S[24]-70*S[26]))))/16;c[34].x=(-105*(7*C[20]+990*(C[22]-12*(C[24]+910*C[26]))))/16;c[34].y=(-105*(7*S[20]-990*(S[22]+12*(S[24]-910*S[26]))))/16;c[34].z=(-105*(3960*C[23]+J[5]))/8;c[35].x=(-315*(97*C[21]+2*(-1408*C[23]+5148*C[25]+J[5])))/4;c[35].y=(4095*S[21])/4-6930*(7*S[23]+117*S[25]);c[35].z=(315*(7*C[20]-275*C[22]+5148*C[24]))/2;c[36].x=(945*(5*S[20]-678*S[22]+30360*S[24]-48048*S[26]))/8;c[36].y=(945*(5*C[20]+6*(-63*C[22]+2948*C[24]+8008*C[26])))/8;c[36].z=(2835*(5*S[21]+88*(-4*S[23]+91*S[25])))/2;c[37].x=(945*(10*C[21]+5368*C[23]-391248*C[25]+J[5]))/16;c[37].y=1890*(S[21]+22*(-7*S[23]+468*S[25]));c[37].z=(-945*(7*C[20]+22*(5*C[22]-468*(3*C[24]-70*C[26]))))/16;c[38].x=(-945*(S[20]-74*S[22]-7480*S[24]+976976*S[26]))/16;c[38].y=(-945*(C[20]-14*C[22]-8536*C[24]+944944*C[26]))/16;c[38].z=(-945*(3*S[21]-88*(S[23]+91*S[25])))/4;c[39].x=(35*(247*C[20]-54*(207*C[22]-5588*C[24]+8008*C[26])))/16;c[39].y=(-35*(23*S[20]-3078*S[22]+2376*(55*S[24]+182*S[26])))/16;c[39].z=(35*(405*C[21]+8*(-1782*C[23]+27027*C[25]+J[5])))/4;c[40].x=(4725*(S[21]-88*(S[23]-39*S[25])))/2;c[40].y=(315*(78*C[21]-6600*C[23]+308880*C[25]+J[5]))/16;c[40].z=(-315*(7*S[20]-990*(S[22]-52*S[24]+728*S[26])))/16;c[41].x=(-35*(23*C[20]+918*C[22]+2376*(-97*C[24]+5642*C[26])))/16;c[41].y=(-35*(5*S[20]-54*(27*S[22]+44*(-85*S[24]+5278*S[26]))))/16;c[41].z=(-35*(27*C[21]+2376*C[23]-216216*C[25]+J[5]))/4;c[42].x=(315*(54*C[21]-2376*C[23]+61776*C[25]+J[5]))/16;c[42].y=(-945*(S[21]-88*(S[23]-39*S[25])))/4;c[42].z=(-315*(7*C[20]-330*C[22]+10296*(C[24]-10*C[26])))/16;c[43].x=(-315*(S[20]-54*(3*S[22]-220*S[24]+8008*S[26])))/16;c[43].y=(-315*(C[20]-102*C[22]+8712*C[24]-391248*C[26]))/16;c[43].z=(-945*(S[21]-88*(S[23]-39*S[25])))/4;c[44].x=(-35*(C[20]-54*(C[22]-44*(C[24]-26*C[26]))))/2;c[44].y=(35*(S[20]-54*(3*S[22]-220*S[24]+8008*S[26])))/16;c[44].z=(-35*(54*C[21]-2376*C[23]+61776*C[25]+J[5]))/16;
        }
        else if(n==8){
            c[0].x=-36*C[27];c[0].y=-36*S[27];c[0].z=-9*J[6];c[1].x=-1260*S[28];c[1].y=45*(28*C[28]-J[6]);c[1].z=360*S[27];c[2].x=9*(17*C[27]+2310*C[29]);c[2].y=9*(127*S[27]+2310*S[29]);c[2].z=-6930*C[28]+162*J[6];c[3].x=1890*(S[28]+110*S[30]);c[3].y=-315*(50*C[28]+660*C[30]-J[6]);c[3].z=-2520*(S[27]+33*S[29]);c[4].x=(63*(C[27]-165*(C[29]+260*C[31])))/2;c[4].y=(-63*(109*S[27]+4455*S[29]+42900*S[31]))/2;c[4].z=(63*(1540*C[28]+42900*C[30]-27*J[6]))/4;c[5].x=(315*(25*S[28]+264*(4*S[30]-273*S[32])))/4;c[5].y=(315*(305*C[28]+10956*C[30]+72072*C[32]-5*J[6]))/4;c[5].z=315*(10*S[27]+495*S[29]+12012*S[31]);c[6].x=(-315*(7*C[27]+66*(17*C[29]+780*(C[31]-14*C[33]))))/16;c[6].y=(315*(103*S[27]+66*(87*S[29]+260*(11*S[31]+42*S[33]))))/16;c[6].z=(-315*(385*C[28]+17160*C[30]+360360*C[32]-6*J[6]))/8;c[7].x=(-315*(29*S[28]+132*(23*S[30]+78*(13*S[32]-40*S[34]))))/8;c[7].y=(-315*(338*C[28]+17952*C[30]+555984*C[32]+823680*C[34]-5*J[6]))/16;c[7].z=(-315*(5*S[27]+33*(9*S[29]+364*S[31]+6240*S[33])))/2;c[8].x=(315*(C[27]+66*(3*C[29]+260*(C[31]+42*C[33]))))/16;c[8].y=(-1575*(S[27]+66*(S[29]+52*(S[31]+30*S[33]))))/8;c[8].z=(315*(616*C[28]+34320*C[30]+1235520*C[32]+14002560*C[34]-9*J[6]))/128;c[9].x=(315*(S[28]+132*(S[30]+78*(S[32]+40*S[34]))))/8;c[9].y=(315*(72*C[28]+4752*C[30]+247104*C[32]+7413120*C[34]-J[6]))/128;c[9].z=(315*(S[27]+66*(S[29]+52*(S[31]+30*S[33]))))/16;c[10].x=-45*(28*C[28]+J[6]);c[10].y=-1260*S[28];c[10].z=360*C[27];c[11].x=990*(S[27]-42*S[29]);c[11].y=990*(C[27]+42*C[29]);c[11].z=13860*S[28];c[12].x=-315*(38*C[28]-1980*C[30]-J[6]);c[12].y=630*(47*S[28]+990*S[30]);c[12].z=-2520*(C[27]+99*C[29]);c[13].x=-3465*(S[27]+36*S[29]-1560*S[31]);c[13].y=-3465*(C[27]+120*(C[29]+13*C[31]));c[13].z=-6930*(7*S[28]+390*S[30]);c[14].x=(1575*(71*C[28]+3036*C[30]-72072*C[32]-J[6]))/4;c[14].y=(-1575*(127*S[28]+8976*S[30]+72072*S[32]))/4;c[14].z=1575*(2*C[27]+33*(9*C[29]+364*C[31]));c[15].x=(3465*(5*S[27]+6*(69*S[29]+3380*S[31]-32760*S[33])))/8;c[15].y=(3465*(5*C[27]+834*C[29]+45240*C[31]+196560*C[33]))/8;c[15].z=(17325*(7*S[28]+312*(2*S[30]+63*S[32])))/4;c[16].x=(-315*(454*C[28]+42240*C[30]+2162160*C[32]-5765760*C[34]-5*J[6]))/16;c[16].y=(315*(367*S[28]+660*(59*S[30]+546*(5*S[32]+8*S[34]))))/8;c[16].z=(-315*(5*C[27]+891*C[29]+60060*(C[31]+24*C[33])))/2;c[17].x=(-3465*(S[27]+114*S[29]+1560*(7*S[31]+354*S[33])))/16;c[17].y=(-3465*(C[27]+66*(3*C[29]+260*(C[31]+42*C[33]))))/16;c[17].z=(-3465*(7*S[28]+780*(S[30]+54*S[32]+816*S[34])))/8;c[18].x=(315*(104*C[28]+13200*C[30]+1235520*C[32]+60128640*C[34]-J[6]))/128;c[18].y=(-1575*(S[28]+132*(S[30]+78*(S[32]+40*S[34]))))/4;c[18].z=(315*(C[27]+66*(3*C[29]+260*(C[31]+42*C[33]))))/16;c[19].x=9*(127*C[27]-2310*C[29]);c[19].y=9*(17*S[27]-2310*S[29]);c[19].z=18*(385*C[28]+9*J[6]);c[20].x=630*(47*S[28]-990*S[30]);c[20].y=315*(38*C[28]+1980*C[30]+J[6]);c[20].z=-2520*(S[27]-99*S[29]);c[21].x=-567*(6*C[27]+715*(C[29]-20*C[31]));c[21].y=-567*(6*S[27]-715*(S[29]+20*S[31]));c[21].z=(-567*(14300*C[30]+3*J[6]))/2;c[22].x=(-1575*(61*S[28]+5016*S[30]-72072*S[32]))/2;c[22].y=(-1575*(5*C[28]+6996*C[30]+72072*C[32]+J[6]))/2;c[22].z=3150*(2*S[27]-33*(3*S[29]+364*S[31]));c[23].x=(315*(89*C[27]+330*(59*C[29]+4420*C[31]-32760*C[33])))/16;c[23].y=(315*(199*S[27]-1650*(9*S[29]+1300*S[31]+6552*S[33])))/16;c[23].z=(-315*(385*C[28]-6*(14300*C[30]+900900*C[32]+3*J[6])))/8;c[24].x=(945*(103*S[28]+396*(43*S[30]+3094*S[32]-7280*S[34])))/8;c[24].y=(-945*(74*C[28]-38016*C[30]-3315312*C[32]-5765760*C[34]-5*J[6]))/16;c[24].z=(-945*(5*S[27]-99*(S[29]+364*(S[31]+40*S[33]))))/2;c[25].x=(-315*(7*C[27]+66*(33*C[29]+4940*C[31]+338520*C[33])))/16;c[25].y=(-315*(29*S[27]-66*(19*S[29]+5980*S[31]+404040*S[33])))/16;c[25].z=(315*(308*C[28]-34320*C[30]-9*(480480*C[32]+10890880*C[34]+J[6])))/32;c[26].x=(-315*(7*S[28]+1716*(S[30]+138*S[32]+8880*S[34])))/8;c[26].y=(315*(28*C[28]-6864*C[30]-1111968*C[32]-65070720*C[34]-J[6]))/32;c[26].z=(315*(S[27]-6864*(S[31]+75*S[33])))/4;c[27].x=315*(50*C[28]-660*C[30]+J[6]);c[27].y=1890*(S[28]-110*S[30]);c[27].z=-2520*(C[27]-33*C[29]);c[28].x=-3465*(S[27]-120*(S[29]-13*S[31]));c[28].y=-3465*(C[27]-36*C[29]-1560*C[31]);c[28].z=-6930*(7*S[28]-390*S[30]);c[29].x=(1575*(5*C[28]-6996*C[30]+72072*C[32]-J[6]))/2;c[29].y=(-1575*(61*S[28]-264*(19*S[30]+273*S[32])))/2;c[29].z=3150*(2*C[27]+99*C[29]-12012*C[31]);c[30].x=(17325*(S[27]-42*(S[29]+260*(S[31]-6*S[33]))))/4;c[30].y=(17325*(C[27]+42*(C[29]-260*(C[31]+6*C[33]))))/4;c[30].z=(121275*(S[28]-4680*S[32]))/2;c[31].x=(-4725*(38*C[28]-5984*C[30]-912912*C[32]+1921920*C[34]-J[6]))/16;c[31].y=(4725*(47*S[28]-44*(23*S[30]+546*(21*S[32]+40*S[34]))))/8;c[31].z=(-4725*(C[27]+99*C[29]-4004*(C[31]+120*C[33])))/2;c[32].x=(-10395*(S[27]+10*(S[29]-52*(23*S[31]+2730*S[33]))))/16;c[32].y=(-10395*(C[27]+94*C[29]-520*(19*C[31]+2982*C[33])))/16;c[32].z=(-10395*(7*S[28]+260*(S[30]-42*(3*S[32]+136*S[34]))))/8;c[33].x=(315*(60*C[28]-2640*C[30]-1441440*C[32]-144144000*C[34]-J[6]))/32;c[33].y=(-315*(29*S[28]+660*(S[30]-546*(S[32]+104*S[34]))))/8;c[33].z=(315*(C[27]+132*(C[29]-5460*C[33])))/4;c[34].x=(-63*(109*C[27]-4455*C[29]+42900*C[31]))/2;c[34].y=(63*(S[27]+165*(S[29]-260*S[31])))/2;c[34].z=(-63*(1540*C[28]-42900*C[30]+27*J[6]))/4;c[35].x=(-1575*(127*S[28]-8976*S[30]+72072*S[32]))/4;c[35].y=(-1575*(71*C[28]-3036*C[30]-72072*C[32]+J[6]))/4;c[35].z=1575*(2*S[27]+33*(-9*S[29]+364*S[31]));c[36].x=(315*(199*C[27]+1650*(9*C[29]-1300*C[31]+6552*C[33])))/16;c[36].y=(315*(89*S[27]-19470*S[29]+85800*(17*S[31]+126*S[33])))/16;c[36].z=(315*(385*C[28]+85800*C[30]-5405400*C[32]+18*J[6]))/8;c[37].x=(4725*(47*S[28]+44*(23*S[30]-546*(21*S[32]-40*S[34]))))/8;c[37].y=(4725*(38*C[28]+5984*C[30]-912912*C[32]-1921920*C[34]+J[6]))/16;c[37].z=(-4725*(S[27]-99*S[29]-4004*(S[31]-120*S[33])))/2;c[38].x=(-945*(9*C[27]+1430*(C[29]-20*(C[31]+714*C[33]))))/16;c[38].y=(-945*(9*S[27]-1430*(S[29]+20*(S[31]-714*S[33]))))/16;c[38].z=(-945*(57200*C[30]-163363200*C[34]+9*J[6]))/64;c[39].x=(-945*(9*S[28]+836*S[30]-24024*(3*S[32]+760*S[34])))/8;c[39].y=(-945*(16*C[28]+9328*C[30]-384384*C[32]-147987840*C[34]+J[6]))/64;c[39].z=(945*(S[27]-66*S[29]-8008*(S[31]-30*S[33])))/8;c[40].x=(-315*(305*C[28]-10956*C[30]+72072*C[32]+5*J[6]))/4;c[40].y=(7875*S[28])/4-20790*(4*S[30]+273*S[32]);c[40].z=315*(10*C[27]-495*C[29]+12012*C[31]);c[41].x=(3465*(5*S[27]-834*S[29]+45240*S[31]-196560*S[33]))/8;c[41].y=(3465*(5*C[27]+6*(-69*C[29]+3380*C[31]+32760*C[33])))/8;c[41].z=(17325*(7*S[28]+312*(-2*S[30]+63*S[32])))/4;c[42].x=(945*(74*C[28]+38016*C[30]-3315312*C[32]+5765760*C[34]+5*J[6]))/16;c[42].y=(945*(103*S[28]-17028*S[30]+72072*(17*S[32]+40*S[34])))/8;c[42].z=(-945*(5*C[27]+99*(C[29]-364*(C[31]-40*C[33]))))/2;c[43].x=(-10395*(S[27]-94*S[29]+520*(-19*S[31]+2982*S[33])))/16;c[43].y=(-10395*(C[27]-10*(C[29]+52*(23*C[31]-2730*C[33]))))/16;c[43].z=(-10395*(7*S[28]-260*(S[30]+126*S[32]-5712*S[34])))/8;c[44].x=(945*(16*C[28]-9328*C[30]-384384*C[32]+147987840*C[34]-J[6]))/64;c[44].y=(-945*(9*S[28]-44*(19*S[30]+546*(3*S[32]-760*S[34]))))/8;c[44].z=(945*(C[27]+66*C[29]-8008*(C[31]+30*C[33])))/8;c[45].x=(315*(103*C[27]-66*(87*C[29]+260*(-11*C[31]+42*C[33]))))/16;c[45].y=(-315*(7*S[27]+66*(-17*S[29]+780*(S[31]+14*S[33]))))/16;c[45].z=(315*(385*C[28]+6*(-2860*C[30]+60060*C[32]+J[6])))/8;c[46].x=(315*(367*S[28]-660*(59*S[30]-546*(5*S[32]-8*S[34]))))/8;c[46].y=(315*(454*C[28]+5*(-8448*C[30]+432432*C[32]+1153152*C[34]+J[6])))/16;c[46].z=(-315*(5*S[27]-891*S[29]+60060*(S[31]-24*S[33])))/2;c[47].x=(-315*(29*C[27]+66*(19*C[29]-5980*C[31]+404040*C[33])))/16;c[47].y=(-315*(7*S[27]-66*(33*S[29]-4940*S[31]+338520*S[33])))/16;c[47].z=(-315*(308*C[28]+34320*C[30]+9*(-480480*C[32]+10890880*C[34]+J[6])))/32;c[48].x=(-315*(29*S[28]-660*(S[30]+546*(S[32]-104*S[34]))))/8;c[48].y=(-315*(60*C[28]+2640*C[30]-1441440*C[32]+144144000*C[34]+J[6]))/32;c[48].z=(315*(S[27]-132*(S[29]-5460*S[33])))/4;c[49].x=(315*(338*C[28]-17952*C[30]+555984*C[32]-823680*C[34]+5*J[6]))/16;c[49].y=(-315*(29*S[28]+132*(-23*S[30]+78*(13*S[32]+40*S[34]))))/8;c[49].z=(315*(-5*C[27]+33*(9*C[29]-364*C[31]+6240*C[33])))/2;c[50].x=(-3465*(S[27]-198*S[29]+17160*(S[31]-42*S[33])))/16;c[50].y=(-3465*(C[27]-6*(19*C[29]+260*(-7*C[31]+354*C[33]))))/16;c[50].z=(-3465*(7*S[28]-780*(S[30]-54*S[32]+816*S[34])))/8;c[51].x=(-315*(28*C[28]+6864*C[30]-1111968*C[32]+65070720*C[34]+J[6]))/32;c[51].y=(-315*(7*S[28]-1716*(S[30]-138*S[32]+8880*S[34])))/8;c[51].z=(315*(C[27]-6864*(C[31]-75*C[33])))/4;c[52].x=(-1575*(C[27]-66*(C[29]-52*(C[31]-30*C[33]))))/8;c[52].y=(315*(S[27]-198*S[29]+17160*(S[31]-42*S[33])))/16;c[52].z=(-315*(616*C[28]-34320*C[30]+9*(137280*C[32]-1555840*C[34]+J[6])))/128;c[53].x=(-1575*(S[28]-132*(S[30]-78*(S[32]-40*S[34]))))/4;c[53].y=(-315*(104*C[28]-13200*C[30]+1235520*C[32]-60128640*C[34]+J[6]))/128;c[53].z=(315*(S[27]-198*S[29]+17160*(S[31]-42*S[33])))/16;c[54].x=(-315*(72*C[28]-4752*C[30]+247104*C[32]-7413120*C[34]+J[6]))/128;c[54].y=(315*(S[28]-132*(S[30]-78*(S[32]-40*S[34]))))/8;c[54].z=(315*(C[27]-66*(C[29]-52*(C[31]-30*C[33]))))/16;
        }
    }
    return ret;
}

void geopotential::unload(geopotential *gp){
    free(gp);
}

int_t geopotential::size() const{
    return sizeof(geopotential)+sizeof(fast_mpvec)*Precompute_Table_size(N);
}

