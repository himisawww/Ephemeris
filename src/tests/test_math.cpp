#include<algorithm>
#include"definitions.h"
#include"math/random.h"
#include"utils/logger.h"
static double max_relative_error=0;

#define TEST_EPSILON_QUATERNION     2e-14
#define TEST_N                      1048576

int test_quaternion(){
    double max_err=0;
    double max_err2=0;
    quat q;
    int_t n_test=0;
    do{
        const double w=randomreal()*2-1;
        if(w==0)
            continue;
        quat q1(w);
        q=w;
        checked_maximize(max_err2,(q-q1).normsqr());
        checked_maximize(max_err,std::abs(q.normsqr()-w*w));
        checked_maximize(max_err,std::abs(q.norm()-std::abs(w)));

        const vec v=randomdirection()*randomreal();
        quat q2(v);
        quat q3(v.x,v.y,v.z);
        q=v;
        checked_maximize(max_err2,(q-q2).normsqr());
        checked_maximize(max_err2,(q-q3).normsqr());
        checked_maximize(max_err,std::abs(q.normsqr()-v.normsqr()));
        checked_maximize(max_err,std::abs(q.norm()-v.norm()));

        quat q4(v,w);
        quat q5(v.x,v.y,v.z,w);
        quat q6(q5);
        q=q4;
        checked_maximize(max_err2,(q-q4).normsqr());
        checked_maximize(max_err2,(q-q5).normsqr());
        checked_maximize(max_err2,(q-q6).normsqr());

        vec vq(q);
        checked_maximize(max_err2,(vq-v).normsqr());

        mat m=randommatrix();
        quat qm(m);
        mat d=mat(qm)-m;
        checked_maximize(max_err2,d.x.normsqr()+d.y.normsqr()+d.z.normsqr());

        double nq2=q.normsqr();
        double nq=q.norm();
        double nqz=q.normalize();
        double nq1=q.norm();
        checked_maximize(max_err2,vec(nq1-1,nq2-nq*nq,nqz-nq).normsqr());

        q3=mat(q);
        q1=q-q3;
        q2=q+q3;
        checked_maximize(max_err2,checked_min(q1.normsqr(),q2.normsqr()));

        std::swap(q,qm);
        qm*=1+w;
        q*=2-w*w;

        q1=q+w;
        q2=q-w;
        q3=q;
        q4=q;
        q3+=w;
        q4-=w;
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr());
        checked_maximize(max_err2,((w+q)-q1).normsqr());
        checked_maximize(max_err2,(q2+(w-q)).normsqr());
        checked_maximize(max_err2,(q1+q2-q*2).normsqr());
        checked_maximize(max_err2,(q1-q2-w*2).normsqr());

        checked_maximize(max_err2,(q1-(q+quat(w))).normsqr());
        checked_maximize(max_err2,(q2-(q-quat(w))).normsqr());
        checked_maximize(max_err2,((quat(w)+q)-q1).normsqr());
        checked_maximize(max_err2,(q2+(quat(w)-q)).normsqr());
        q3=q;
        q4=q;
        q3+=quat(w);
        q4-=quat(w);
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr());

        q1=q*w;
        q2=q/w;
        q3=q;
        q4=q;
        q3*=w;
        q4/=w;
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr()*(w*w));
        checked_maximize(max_err2,((w*q)-q1).normsqr());
        checked_maximize(max_err2,(q2-(w/q).inverse()).normsqr()*(w*w));
        checked_maximize(max_err2,(q1*q2-q*q).normsqr());
        checked_maximize(max_err2,(q1/q2-w*w).normsqr());

        checked_maximize(max_err2,(q1-(q*quat(w))).normsqr());
        checked_maximize(max_err2,(q2-(q/quat(w))).normsqr()*(w*w));
        checked_maximize(max_err2,((quat(w)*q)-q1).normsqr());
        checked_maximize(max_err2,(q2-(quat(w)/q).inverse()).normsqr()*(w*w));
        q3=q;
        q4=q;
        q3*=quat(w);
        q4/=quat(w);
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr()*(w*w));

        q1=q+v;
        q2=q-v;
        q3=q;
        q4=q;
        q3+=v;
        q4-=v;
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr());
        checked_maximize(max_err2,((v+q)-q1).normsqr());
        checked_maximize(max_err2,(q2+(v-q)).normsqr());
        checked_maximize(max_err2,(q1+q2-q*2).normsqr());
        checked_maximize(max_err2,(q1-q2-v*2).normsqr());

        checked_maximize(max_err2,(q1-(q+quat(v))).normsqr());
        checked_maximize(max_err2,(q2-(q-quat(v))).normsqr());
        checked_maximize(max_err2,((quat(v)+q)-q1).normsqr());
        checked_maximize(max_err2,(q2+(quat(v)-q)).normsqr());
        q3=q;
        q4=q;
        q3+=quat(v);
        q4-=quat(v);
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr());

        q1=q*v;
        q2=q/v;
        q3=q;
        q4=q;
        q3*=v;
        q4/=v;
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr()*v.normsqr());
        checked_maximize(max_err2,(v*q.conjugate()+q1.conjugate()).normsqr());
        checked_maximize(max_err2,(q2-(v/q).inverse()).normsqr()*v.normsqr());
        checked_maximize(max_err2,(q2*q1.conjugate()+q.normsqr()).normsqr());
        checked_maximize(max_err2,(q1/q2+v.normsqr()).normsqr());

        checked_maximize(max_err2,(q1-(q*quat(v))).normsqr());
        checked_maximize(max_err2,(q2-(q/quat(v))).normsqr()*v.normsqr());
        checked_maximize(max_err2,(quat(v)*q.conjugate()+q1.conjugate()).normsqr());
        checked_maximize(max_err2,(q2-(quat(v)/q).inverse()).normsqr()*v.normsqr());
        q3=q;
        q4=q;
        q3*=quat(v);
        q4/=quat(v);
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr()*v.normsqr());

        q1=q+qm;
        q2=q-qm;
        q3=q;
        q4=q;
        q3+=qm;
        q4-=qm;
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr());
        checked_maximize(max_err2,(q1+q2-q*2).normsqr());
        checked_maximize(max_err2,(q1-q2-qm*2).normsqr());

        q1=q*qm;
        q2=q/qm;
        q3=q;
        q4=q;
        q3*=qm;
        q4/=qm;
        checked_maximize(max_err2,(q3-q1).normsqr());
        checked_maximize(max_err2,(q2-q4).normsqr());
        checked_maximize(max_err2,(q1*qm.inverse()-q).normsqr());
        checked_maximize(max_err2,(q2/qm.inverse()-q).normsqr());

        q5=qm*q;
        checked_maximize(max_err2,(q5-q1-2*(vec(qm)*vec(q))).normsqr());
        checked_maximize(max_err,std::abs(q5.normsqr()-qm.normsqr()*q.normsqr()));

        checked_maximize(max_err2,(q-(+q)).normsqr());
        checked_maximize(max_err2,(q+(-q)).normsqr());

        quat eq=exp(q);
        checked_maximize(max_err2,(q-log(eq)).normsqr());

        quat lq=log(q);
        checked_maximize(max_err2,(q-exp(lq)).normsqr());

        quat ew=exp(quat(w));
        checked_maximize(max_err2,(w-log(ew)).normsqr());

        quat lw=log(quat(w));
        checked_maximize(max_err2,(w-exp(lw)).normsqr());

        quat ev=exp(quat(v));
        checked_maximize(max_err2,(v-log(ev)).normsqr());

        quat lv=log(quat(v));
        checked_maximize(max_err2,(v-exp(lv)).normsqr());

        checked_maximize(max_err,std::abs(m.normsqr()-3));
        m+=randommatrix()*randomreal();
        checked_maximize(max_err,std::abs((m.transpose()%m).tr()-m.normsqr()));

        double mdet=m.det();
        if(mdet!=0){
            checked_maximize(max_err2,(m.adjoint()/mdet-m.inverse()).normsqr());
            checked_maximize(max_err2,(w*m.inverse()-w/m).normsqr()*(mdet*mdet));
        }

        double t=randomreal()*Constants::pi_mul2;
        quat theta(v*(w*t/2),0);
        m=1+v.rotation_matrix(w*t);
        qm=m;
        double vn=v.norm();
        t*=w*vn/2;
        vec nv=v/vn;
        quat qv(nv*std::sin(t),std::cos(t));
        checked_maximize(max_err2,checked_min((qv-qm).normsqr(),(qv+qm).normsqr()));
        checked_maximize(max_err2,(qv-exp(theta)).normsqr());

        checked_maximize(max_err,std::sqrt(max_err2));
        if(!(max_err<TEST_EPSILON_QUATERNION)){
            LogError(
                "\nMax Rel. Error %.16le Too Large",
                max_err);
            return 1;
        }
        checked_maximize(max_relative_error,max_err/TEST_EPSILON_QUATERNION);
    } while(++n_test<TEST_N);
    return 0;
}

int test_math(){
    int ret;
    ret=test_quaternion();
    if(ret)return ret;

    LogInfo("\n      Passed(%f), ",max_relative_error);
    return 0;
}
