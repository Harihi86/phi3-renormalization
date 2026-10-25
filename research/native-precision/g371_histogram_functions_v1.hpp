#include "g371_soft_functions_v1.hpp"
struct HistogramMappedPoint{std::array<Real,7> point;Real jacobian;};
const std::array<std::array<long,16>,7> histogram_soft_pair={{{870,870,871,872,880,902,966,1085,1323,1671,2228,3211,5031,8320,14059,22377},{870,871,872,880,917,1037,1357,2123,3248,4782,6133,7145,7784,8911,9087,9519},{3997,3925,4461,3898,3883,3799,4376,3779,4051,4556,4252,4194,4091,4028,4409,3837},{3708,4556,4232,4032,4251,3962,4056,3705,4086,4033,4077,3955,4436,4766,3676,4005},{4265,4222,3582,4326,3772,3771,4134,3696,3793,4693,4227,4286,4349,4642,3589,4189},{4778,4068,3985,3856,4744,4390,4235,3725,4075,4114,4245,3944,3458,3746,3989,4184},{876,942,1345,2359,4198,6698,8133,9169,9180,7384,5964,3839,2212,1394,969,874}}};
const std::array<std::array<long,16>,7> histogram_soft_scale={{{870,870,870,870,871,871,871,871,876,905,1047,1592,3552,7638,15904,27058},{3812,2496,2464,2445,2395,2232,2109,2067,1648,1453,1730,2357,3427,4959,7833,22109},{5384,3219,3249,3237,6482,4301,3514,4129,3332,3382,4009,4682,4821,4021,3916,3858},{3519,3478,3732,3465,4912,5143,4245,5975,3745,4366,4740,3709,3344,3859,3858,3446},{3768,4869,3313,3870,3329,3574,5297,3563,3641,3467,4659,3611,3687,4320,3879,6689},{3858,3602,3709,3749,3761,5078,3551,3336,3422,3952,3559,5515,5539,4428,4606,3871},{875,912,1346,2189,3920,6447,7230,9180,8465,8468,5226,4438,3626,1385,953,876}}};

HistogramMappedPoint histogram_map(bool scale,const std::array<Real,7>&p){
    const auto&counts=scale?histogram_soft_scale:histogram_soft_pair;HistogramMappedPoint result;result.jacobian=1;
    for(int axis=0;axis<7;++axis){long cumulative=0;bool found=false;
        for(int bin=0;bin<16;++bin){long mass=counts[axis][bin];Real upper=Real(cumulative+mass)/65536;
            if(arb_lt(p[axis].v,upper.v)){Real slope=Real(65536)/Real(16*mass);result.point[axis]=Real(bin)/16+(p[axis]-Real(cumulative)/65536)*slope;result.jacobian=result.jacobian*slope;found=true;break;}
            Real distance=p[axis]-upper;if(!arb_is_nonnegative(distance.v))throw std::runtime_error("uncertain exact histogram cell; higher arithmetic precision required");cumulative+=mass;}
        if(!found)throw std::invalid_argument("strict histogram cube required");pos(result.point[axis],"positive mapped coordinate");pos(1-result.point[axis],"mapped coordinate below one");}
    pos(result.jacobian,"positive complete histogram Jacobian");return result;
}
Real histogram_host_value(const DoubleHostGraph&g,bool positive,bool scale,const std::array<Real,7>&p){
    auto mapped=histogram_map(scale,p);return mapped.jacobian*soft_host_value(g,positive,scale,mapped.point);
}

