#include<cmath>
#include<vector>
#include<random>
#include<algorithm>

/*
    partition data to n subsets, s.t. sum of each subset are as close as possible to each other
    store index lists of subsets into res
    return max difference between sum(subset) and sum(data)/n

    Note:
        val_t should be signed type
        sum(data) should not overflow val_t
        elements in data should be all positive
        n>=1
*/
template<typename val_t>
val_t distribute(std::vector<std::vector<size_t>> &res,const std::vector<val_t> &data,size_t n){
    typedef std::vector<val_t> val_lt;
    typedef std::vector<size_t> index_lt;
    const size_t ndat=data.size();
    const size_t npos=-1;

    std::random_device rd;
    std::mt19937_64 g(rd());
    res=std::vector<index_lt>(n);

    val_t dsum=(val_t)0;
    val_lt psum(n,(val_t)0);
    val_t perfect;
    for(size_t i=0;i<ndat;++i){
        size_t ii=g()%n;
        res[ii].push_back(i);
        psum[ii]+=data[i];
        dsum+=data[i];
    }
    for(size_t i=0;i<n;++i){
        res[i].push_back(npos);
    }
    perfect=dsum/n;

label1:
    size_t maxdsi;
    val_t maxds=(val_t)0;
    index_lt tryj;
    for(size_t i=0;i<n;++i){
        val_t ds=std::abs(psum[i]-perfect);
        tryj.push_back(i);
        if(ds>=maxds){
            maxds=ds;
            maxdsi=i;
        }
    }

    std::shuffle(tryj.begin(),tryj.end(),g);

    for(auto j:tryj)if(j!=maxdsi){
        index_lt &idxi=res[maxdsi],&idxj=res[j];

        std::shuffle(idxi.begin(),idxi.end(),g);
        std::shuffle(idxj.begin(),idxj.end(),g);

        for(auto &ii:idxi)for(auto &ij:idxj){
            val_t resi=psum[maxdsi],resj=psum[j];
            if(ii==npos&&ij==npos)continue;
            if(ii!=npos){
                resi-=data[ii];
                resj+=data[ii];
            }
            if(ij!=npos){
                resi+=data[ij];
                resj-=data[ij];
            }
            if(std::abs(resi-perfect)>=maxds||std::abs(resj-perfect)>=maxds)continue;
            psum[maxdsi]=resi;
            psum[j]=resj;
            if(ii!=npos&&ij!=npos){
                auto temp=ii;
                ii=ij;
                ij=temp;
            }
            else if(ii==npos){
                idxi.push_back(ij);
                idxj.erase(idxj.begin()+(&ij-&idxj[0]));
            }
            else{
                idxj.push_back(ii);
                idxi.erase(idxi.begin()+(&ii-&idxi[0]));
            }
            goto label1;
        }
    }

    for(auto &l:res){
        std::sort(l.begin(),l.end());
        l.erase(l.end()-1);
    }
    return maxds;
}

