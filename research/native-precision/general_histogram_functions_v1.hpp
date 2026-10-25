struct GeneralMappedPoint{std::array<Real,7> point;Real jacobian;};
const std::array<std::array<long,16>,7> histogram_L4E3G244_angular_pair={{{870,870,871,873,884,919,1015,1245,1577,2070,2545,3637,5335,8701,14591,19533},{870,871,872,880,912,1018,1297,1921,2934,4280,5442,6932,8250,9036,10073,9948},{3707,3802,4524,4076,4367,4371,3923,4251,3760,3712,3927,4618,3531,4481,3993,4493},{3734,4049,3646,3558,4123,4004,4299,4002,4162,4843,4256,3873,4290,4117,4596,3984},{4319,3593,4004,4132,4255,4489,4494,3869,3970,3498,4387,4089,3762,4009,4341,4325},{3990,3979,3949,4247,3987,4334,4305,4068,4072,4488,3970,3658,4008,3804,4256,4421},{874,956,1297,2863,3715,6190,7930,9209,8582,8105,6205,3970,2520,1301,944,875}}};
const std::array<std::array<long,16>,7> histogram_L4E3G244_angular_scale={{{870,870,870,870,870,870,871,871,874,892,1026,1729,3646,8170,15691,26546},{6107,2513,2074,1929,1916,1908,1993,2217,2517,2718,3006,3465,4288,5040,6827,17018},{3226,3722,5075,7224,4249,3704,3611,3582,3875,3492,4068,3791,4217,3732,4232,3736},{3238,3202,3949,3210,6812,5037,3637,3889,4353,3920,3675,5028,4026,3805,4468,3287},{3823,3642,3142,3645,3693,5272,3606,4220,3657,3570,4385,4370,3727,4099,4620,6065},{4072,4310,3334,3406,3941,4015,3890,5916,3750,5207,4519,3762,4758,3613,3530,3513},{873,936,1205,2354,4104,5622,8622,8809,10759,7471,5558,3572,2341,1486,947,877}}};
const std::array<std::array<long,16>,7> histogram_L4E3G370_angular_pair={{{870,870,871,875,894,954,1125,1460,2121,3141,4143,5470,7244,9442,12204,13852},{870,870,871,879,906,993,1244,1745,2644,3867,4775,6349,7792,8866,10732,12133},{3783,3858,4165,4353,4320,4305,4177,4076,4003,3981,3837,4254,3916,4286,3859,4363},{3367,3912,3643,3572,4224,4255,4509,4083,4161,4503,4424,4062,4315,4193,4530,3783},{4353,3792,4000,4050,4348,4314,4447,3883,3836,3743,4360,4301,4048,4158,3882,4021},{4166,4118,3821,4040,4400,4193,4333,4283,3854,4434,4082,3913,3819,3760,4033,4287},{873,959,1288,2476,3315,5435,7120,9091,9080,8690,6866,4451,2641,1411,963,877}}};
const std::array<std::array<long,16>,7> histogram_L4E3G370_angular_scale={{{870,870,870,870,870,870,871,871,875,914,1169,2208,4530,8919,15383,24576},{8562,4049,3332,2907,2685,2630,2539,2694,2606,2682,2700,2972,3407,3668,4791,13312},{3742,3732,4858,7087,3888,3941,3517,3595,4376,3446,3682,3885,4649,3607,3941,3590},{3882,3178,3556,3358,6677,4957,3717,3700,4238,4218,3772,5141,3980,3756,3990,3416},{3491,4458,3484,3684,3669,4950,3387,4019,3575,3526,4669,4136,4347,4037,4028,6076},{4420,4145,3482,3500,4431,4062,3848,5897,4260,4638,4024,3499,4602,3706,3489,3533},{873,928,1180,1970,4154,5613,8036,8157,10484,8118,6265,3927,2580,1411,963,877}}};
const std::array<std::array<long,16>,7>& general_histogram_counts(const std::string&id,bool scale){
if(id=="L4E3G244")return scale?histogram_L4E3G244_angular_scale:histogram_L4E3G244_angular_pair;
if(id=="L4E3G370")return scale?histogram_L4E3G370_angular_scale:histogram_L4E3G370_angular_pair;
throw std::invalid_argument("untrained histogram graph");}

GeneralMappedPoint general_histogram_map(const std::string&id,bool scale,const std::array<Real,7>&p){
    const auto&counts=general_histogram_counts(id,scale);GeneralMappedPoint result;result.jacobian=1;
    for(int axis=0;axis<7;++axis){long cumulative=0;bool found=false;
        for(int bin=0;bin<16;++bin){long mass=counts[axis][bin];Real upper=Real(cumulative+mass)/65536;
            if(arb_lt(p[axis].v,upper.v)){Real slope=Real(65536)/Real(16*mass);result.point[axis]=Real(bin)/16+(p[axis]-Real(cumulative)/65536)*slope;result.jacobian=result.jacobian*slope;found=true;break;}
            Real distance=p[axis]-upper;if(!arb_is_nonnegative(distance.v))throw std::runtime_error("uncertain exact histogram cell; higher arithmetic precision required");cumulative+=mass;}
        if(!found)throw std::invalid_argument("strict histogram cube required");pos(result.point[axis],"positive mapped coordinate");pos(1-result.point[axis],"mapped coordinate below one");}
    pos(result.jacobian,"positive complete histogram Jacobian");return result;
}
Real general_histogram_value(const DoubleHostGraph&g,bool positive,bool scale,const std::array<Real,7>&p){
    auto mapped=general_histogram_map(g.id,scale,p);return mapped.jacobian*general_angular_value(g,positive,scale,mapped.point);
}

