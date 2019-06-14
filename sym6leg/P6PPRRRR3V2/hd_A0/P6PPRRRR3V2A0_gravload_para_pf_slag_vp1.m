% Calculate Gravitation load for parallel robot
% P6PPRRRR3V2A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [6x1]
%   Generalized platform coordinates
% qJ [3x6]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [6x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% koppelP [6x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,alpha3,alpha4,d3,d4,theta1,theta2]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [6x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-16 22:25
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P6PPRRRR3V2A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,6),zeros(3,1),zeros(6,3),zeros(6,3),zeros(10,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [6 1]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: xP has to be [6x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 6]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: qJ has to be [3x6] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [6 3]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: legFrame has to be [6x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [6 3]), ...
  'P6PPRRRR3V2A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [6x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-16 21:45:29
% EndTime: 2019-05-16 21:45:44
% DurationCPUTime: 15.78s
% Computational Cost: add. (8717->625), mult. (22695->1219), div. (270->13), fcn. (25907->63), ass. (0->458)
t3998 = sin(pkin(10));
t3999 = sin(pkin(9));
t4004 = cos(pkin(9));
t4005 = cos(pkin(5));
t4001 = sin(pkin(5));
t4002 = sin(pkin(4));
t4112 = t4001 * t4002;
t4081 = t3999 * t4112;
t4003 = cos(pkin(10));
t4006 = cos(pkin(4));
t4102 = t4003 * t4006;
t3879 = (t3998 * t4004 + t3999 * t4102) * t4005 - t4081;
t4082 = t4004 * t4112;
t4101 = t4004 * t4006;
t3880 = (t3998 * t3999 - t4003 * t4101) * t4005 + t4082;
t4007 = legFrame(6,3);
t3954 = sin(t4007);
t3966 = cos(t4007);
t4261 = t3879 * t3954 + t3880 * t3966;
t4008 = legFrame(5,3);
t3955 = sin(t4008);
t3967 = cos(t4008);
t4260 = t3879 * t3955 + t3880 * t3967;
t4009 = legFrame(4,3);
t3956 = sin(t4009);
t3968 = cos(t4009);
t4259 = t3879 * t3956 + t3880 * t3968;
t4010 = legFrame(3,3);
t3957 = sin(t4010);
t3969 = cos(t4010);
t4258 = t3879 * t3957 + t3880 * t3969;
t4011 = legFrame(2,3);
t3958 = sin(t4011);
t3970 = cos(t4011);
t4257 = t3879 * t3958 + t3880 * t3970;
t4012 = legFrame(1,3);
t3959 = sin(t4012);
t3971 = cos(t4012);
t4256 = t3879 * t3959 + t3880 * t3971;
t4039 = xP(5);
t3991 = sin(t4039);
t3994 = cos(t4039);
t4044 = koppelP(6,3);
t4038 = xP(6);
t3990 = sin(t4038);
t3993 = cos(t4038);
t4050 = koppelP(6,1);
t4231 = koppelP(6,2);
t4248 = t3990 * t4231 - t3993 * t4050;
t3895 = t3991 * t4248 + t3994 * t4044;
t3930 = t3990 * t4050 + t3993 * t4231;
t4040 = xP(4);
t3992 = sin(t4040);
t3995 = cos(t4040);
t3806 = t3895 * t3992 - t3930 * t3995;
t3807 = t3895 * t3995 + t3930 * t3992;
t4045 = koppelP(5,3);
t4051 = koppelP(5,1);
t4232 = koppelP(5,2);
t4247 = t3990 * t4232 - t3993 * t4051;
t3897 = t3991 * t4247 + t3994 * t4045;
t3931 = t3990 * t4051 + t3993 * t4232;
t3810 = t3897 * t3992 - t3931 * t3995;
t3811 = t3897 * t3995 + t3931 * t3992;
t4046 = koppelP(4,3);
t4052 = koppelP(4,1);
t4233 = koppelP(4,2);
t4246 = t3990 * t4233 - t3993 * t4052;
t3899 = t3991 * t4246 + t3994 * t4046;
t3932 = t3990 * t4052 + t3993 * t4233;
t3814 = t3899 * t3992 - t3932 * t3995;
t3815 = t3899 * t3995 + t3932 * t3992;
t4047 = koppelP(3,3);
t4053 = koppelP(3,1);
t4234 = koppelP(3,2);
t4245 = t3990 * t4234 - t3993 * t4053;
t3901 = t3991 * t4245 + t3994 * t4047;
t3933 = t3990 * t4053 + t3993 * t4234;
t3818 = t3901 * t3992 - t3933 * t3995;
t3819 = t3901 * t3995 + t3933 * t3992;
t4048 = koppelP(2,3);
t4054 = koppelP(2,1);
t4235 = koppelP(2,2);
t4244 = t3990 * t4235 - t3993 * t4054;
t3903 = t3991 * t4244 + t3994 * t4048;
t3934 = t3990 * t4054 + t3993 * t4235;
t3822 = t3903 * t3992 - t3934 * t3995;
t3823 = t3903 * t3995 + t3934 * t3992;
t4049 = koppelP(1,3);
t4055 = koppelP(1,1);
t4236 = koppelP(1,2);
t4243 = t3990 * t4236 - t3993 * t4055;
t3905 = t3991 * t4243 + t3994 * t4049;
t3935 = t3990 * t4055 + t3993 * t4236;
t3826 = t3905 * t3992 - t3935 * t3995;
t3827 = t3905 * t3995 + t3935 * t3992;
t4041 = rSges(4,3);
t4042 = rSges(4,2);
t4043 = rSges(4,1);
t4068 = t3990 * t4042 - t3993 * t4043;
t4255 = -t3991 * t4041 + t3994 * t4068;
t3909 = t4005 * t4102 - t4112;
t4127 = t3998 * t4005;
t3878 = t3909 * t3999 + t4004 * t4127;
t3881 = t3909 * t4004 - t3999 * t4127;
t3787 = t3878 * t3971 + t3881 * t3959;
t3793 = -t3878 * t3959 + t3881 * t3971;
t4018 = legFrame(1,1);
t3977 = cos(t4018);
t3965 = sin(t4018);
t4030 = legFrame(1,2);
t3983 = sin(t4030);
t4146 = t3965 * t3983;
t4254 = t3787 * t3977 + t3793 * t4146;
t3786 = t3878 * t3970 + t3881 * t3958;
t3792 = -t3878 * t3958 + t3881 * t3970;
t4017 = legFrame(2,1);
t3976 = cos(t4017);
t3964 = sin(t4017);
t4029 = legFrame(2,2);
t3982 = sin(t4029);
t4147 = t3964 * t3982;
t4253 = t3786 * t3976 + t3792 * t4147;
t3785 = t3878 * t3969 + t3881 * t3957;
t3791 = -t3878 * t3957 + t3881 * t3969;
t4016 = legFrame(3,1);
t3975 = cos(t4016);
t3963 = sin(t4016);
t4028 = legFrame(3,2);
t3981 = sin(t4028);
t4148 = t3963 * t3981;
t4252 = t3785 * t3975 + t3791 * t4148;
t3784 = t3878 * t3968 + t3881 * t3956;
t3790 = -t3878 * t3956 + t3881 * t3968;
t4015 = legFrame(4,1);
t3974 = cos(t4015);
t3962 = sin(t4015);
t4027 = legFrame(4,2);
t3980 = sin(t4027);
t4149 = t3962 * t3980;
t4251 = t3784 * t3974 + t3790 * t4149;
t3783 = t3878 * t3967 + t3881 * t3955;
t3789 = -t3878 * t3955 + t3881 * t3967;
t4014 = legFrame(5,1);
t3973 = cos(t4014);
t3961 = sin(t4014);
t4026 = legFrame(5,2);
t3979 = sin(t4026);
t4150 = t3961 * t3979;
t4250 = t3783 * t3973 + t3789 * t4150;
t3782 = t3878 * t3966 + t3881 * t3954;
t3788 = -t3878 * t3954 + t3881 * t3966;
t4013 = legFrame(6,1);
t3972 = cos(t4013);
t3960 = sin(t4013);
t4025 = legFrame(6,2);
t3978 = sin(t4025);
t4151 = t3960 * t3978;
t4249 = t3782 * t3972 + t3788 * t4151;
t4145 = t3972 * t3978;
t3984 = cos(t4025);
t4229 = g(1) * t3984;
t3794 = -t3954 * t4229 + (-t3954 * t4151 + t3966 * t3972) * g(2) + (t3954 * t4145 + t3960 * t3966) * g(3);
t3795 = t3966 * t4229 + (t3954 * t3972 + t3966 * t4151) * g(2) + (t3954 * t3960 - t3966 * t4145) * g(3);
t3758 = t3794 * t4004 - t3795 * t3999;
t3889 = t3978 * g(1) + (-g(2) * t3960 + g(3) * t3972) * t3984;
t4242 = t3758 * t4006 + t3889 * t4002;
t4144 = t3973 * t3979;
t3985 = cos(t4026);
t4228 = g(1) * t3985;
t3796 = -t3955 * t4228 + (-t3955 * t4150 + t3967 * t3973) * g(2) + (t3955 * t4144 + t3961 * t3967) * g(3);
t3797 = t3967 * t4228 + (t3955 * t3973 + t3967 * t4150) * g(2) + (t3955 * t3961 - t3967 * t4144) * g(3);
t3759 = t3796 * t4004 - t3797 * t3999;
t3890 = t3979 * g(1) + (-g(2) * t3961 + g(3) * t3973) * t3985;
t4241 = t3759 * t4006 + t3890 * t4002;
t4143 = t3974 * t3980;
t3986 = cos(t4027);
t4227 = g(1) * t3986;
t3798 = -t3956 * t4227 + (-t3956 * t4149 + t3968 * t3974) * g(2) + (t3956 * t4143 + t3962 * t3968) * g(3);
t3799 = t3968 * t4227 + (t3956 * t3974 + t3968 * t4149) * g(2) + (t3956 * t3962 - t3968 * t4143) * g(3);
t3760 = t3798 * t4004 - t3799 * t3999;
t3891 = t3980 * g(1) + (-g(2) * t3962 + g(3) * t3974) * t3986;
t4240 = t3760 * t4006 + t3891 * t4002;
t4142 = t3975 * t3981;
t3987 = cos(t4028);
t4226 = g(1) * t3987;
t3800 = -t3957 * t4226 + (-t3957 * t4148 + t3969 * t3975) * g(2) + (t3957 * t4142 + t3963 * t3969) * g(3);
t3801 = t3969 * t4226 + (t3957 * t3975 + t3969 * t4148) * g(2) + (t3957 * t3963 - t3969 * t4142) * g(3);
t3761 = t3800 * t4004 - t3801 * t3999;
t3892 = t3981 * g(1) + (-g(2) * t3963 + g(3) * t3975) * t3987;
t4239 = t3761 * t4006 + t3892 * t4002;
t4141 = t3976 * t3982;
t3988 = cos(t4029);
t4225 = g(1) * t3988;
t3802 = -t3958 * t4225 + (-t3958 * t4147 + t3970 * t3976) * g(2) + (t3958 * t4141 + t3964 * t3970) * g(3);
t3803 = t3970 * t4225 + (t3958 * t3976 + t3970 * t4147) * g(2) + (t3958 * t3964 - t3970 * t4141) * g(3);
t3762 = t3802 * t4004 - t3803 * t3999;
t3893 = t3982 * g(1) + (-g(2) * t3964 + g(3) * t3976) * t3988;
t4238 = t3762 * t4006 + t3893 * t4002;
t4140 = t3977 * t3983;
t3989 = cos(t4030);
t4224 = g(1) * t3989;
t3804 = -t3959 * t4224 + (-t3959 * t4146 + t3971 * t3977) * g(2) + (t3959 * t4140 + t3965 * t3971) * g(3);
t3805 = t3971 * t4224 + (t3959 * t3977 + t3971 * t4146) * g(2) + (t3959 * t3965 - t3971 * t4140) * g(3);
t3763 = t3804 * t4004 - t3805 * t3999;
t3894 = t3983 * g(1) + (-g(2) * t3965 + g(3) * t3977) * t3989;
t4237 = t3763 * t4006 + t3894 * t4002;
t4037 = m(2) + m(3);
t4230 = pkin(8) * sin(pkin(6));
t4019 = sin(qJ(3,6));
t4022 = cos(qJ(3,6));
t3936 = pkin(3) * t4022 + t4019 * t4230;
t3939 = pkin(3) * t4019 - t4022 * t4230;
t3869 = t3936 * t4127 + t3939 * t4003;
t3866 = 0.1e1 / t3869;
t4103 = t4003 * t4005;
t4105 = t4002 * t4003;
t4061 = (t3758 * t4103 + t3889 * t4001) * t4006 - (-t3889 * t4105 + (t3794 * t3999 + t3795 * t4004) * t3998) * t4005 - t3794 * t4082 + t3795 * t4081;
t4104 = t4003 * t4004;
t4113 = t3999 * t4003;
t4067 = t3794 * t4113 + t3795 * t4104;
t4122 = t3998 * t4022;
t4125 = t3998 * t4019;
t4223 = ((-t4061 * t4019 - t4067 * t4022 - t4122 * t4242) * rSges(3,2) + (-t4067 * t4019 + t4061 * t4022 - t4125 * t4242) * rSges(3,1)) * t3866;
t4020 = sin(qJ(3,5));
t4023 = cos(qJ(3,5));
t3937 = pkin(3) * t4023 + t4020 * t4230;
t3940 = pkin(3) * t4020 - t4023 * t4230;
t3870 = t3937 * t4127 + t3940 * t4003;
t3867 = 0.1e1 / t3870;
t4060 = (t3759 * t4103 + t3890 * t4001) * t4006 - (-t3890 * t4105 + (t3796 * t3999 + t3797 * t4004) * t3998) * t4005 - t3796 * t4082 + t3797 * t4081;
t4066 = t3796 * t4113 + t3797 * t4104;
t4121 = t3998 * t4023;
t4124 = t3998 * t4020;
t4222 = ((-t4060 * t4020 - t4066 * t4023 - t4121 * t4241) * rSges(3,2) + (-t4066 * t4020 + t4060 * t4023 - t4124 * t4241) * rSges(3,1)) * t3867;
t4021 = sin(qJ(3,4));
t4024 = cos(qJ(3,4));
t3938 = pkin(3) * t4024 + t4021 * t4230;
t3941 = pkin(3) * t4021 - t4024 * t4230;
t3871 = t3938 * t4127 + t3941 * t4003;
t3868 = 0.1e1 / t3871;
t4059 = (t3760 * t4103 + t3891 * t4001) * t4006 - (-t3891 * t4105 + (t3798 * t3999 + t3799 * t4004) * t3998) * t4005 - t3798 * t4082 + t3799 * t4081;
t4065 = t3798 * t4113 + t3799 * t4104;
t4120 = t3998 * t4024;
t4123 = t3998 * t4021;
t4221 = ((-t4059 * t4021 - t4065 * t4024 - t4120 * t4240) * rSges(3,2) + (-t4065 * t4021 + t4059 * t4024 - t4123 * t4240) * rSges(3,1)) * t3868;
t4031 = sin(qJ(3,3));
t4034 = cos(qJ(3,3));
t3942 = pkin(3) * t4034 + t4031 * t4230;
t3945 = pkin(3) * t4031 - t4034 * t4230;
t3875 = t3942 * t4127 + t3945 * t4003;
t3872 = 0.1e1 / t3875;
t4058 = (t3761 * t4103 + t3892 * t4001) * t4006 - (-t3892 * t4105 + (t3800 * t3999 + t3801 * t4004) * t3998) * t4005 - t3800 * t4082 + t3801 * t4081;
t4064 = t3800 * t4113 + t3801 * t4104;
t4116 = t3998 * t4034;
t4119 = t3998 * t4031;
t4220 = ((-t4058 * t4031 - t4064 * t4034 - t4116 * t4239) * rSges(3,2) + (-t4064 * t4031 + t4058 * t4034 - t4119 * t4239) * rSges(3,1)) * t3872;
t4032 = sin(qJ(3,2));
t4035 = cos(qJ(3,2));
t3943 = pkin(3) * t4035 + t4032 * t4230;
t3946 = pkin(3) * t4032 - t4035 * t4230;
t3876 = t3943 * t4127 + t3946 * t4003;
t3873 = 0.1e1 / t3876;
t4057 = (t3762 * t4103 + t3893 * t4001) * t4006 - (-t3893 * t4105 + (t3802 * t3999 + t3803 * t4004) * t3998) * t4005 - t3802 * t4082 + t3803 * t4081;
t4063 = t3802 * t4113 + t3803 * t4104;
t4115 = t3998 * t4035;
t4118 = t3998 * t4032;
t4219 = ((-t4057 * t4032 - t4063 * t4035 - t4115 * t4238) * rSges(3,2) + (-t4063 * t4032 + t4057 * t4035 - t4118 * t4238) * rSges(3,1)) * t3873;
t4033 = sin(qJ(3,1));
t4036 = cos(qJ(3,1));
t3944 = pkin(3) * t4036 + t4033 * t4230;
t3947 = pkin(3) * t4033 - t4036 * t4230;
t3877 = t3944 * t4127 + t3947 * t4003;
t3874 = 0.1e1 / t3877;
t4056 = (t3763 * t4103 + t3894 * t4001) * t4006 - (-t3894 * t4105 + (t3804 * t3999 + t3805 * t4004) * t3998) * t4005 - t3804 * t4082 + t3805 * t4081;
t4062 = t3804 * t4113 + t3805 * t4104;
t4114 = t3998 * t4036;
t4117 = t3998 * t4033;
t4218 = ((-t4056 * t4033 - t4062 * t4036 - t4114 * t4237) * rSges(3,2) + (-t4062 * t4033 + t4056 * t4036 - t4117 * t4237) * rSges(3,1)) * t3874;
t4126 = t3998 * t4006;
t3907 = t3999 * t4126 - t4104;
t3910 = t3998 * t4101 + t4113;
t3848 = -t3907 * t3954 + t3910 * t3966;
t4217 = (-(t3848 * t4022 - t4019 * t4261) * t4230 + pkin(3) * (t4019 * t3848 + t4022 * t4261)) * t3984;
t3849 = -t3907 * t3955 + t3910 * t3967;
t4216 = (-(t3849 * t4023 - t4020 * t4260) * t4230 + pkin(3) * (t4020 * t3849 + t4023 * t4260)) * t3985;
t3850 = -t3907 * t3956 + t3910 * t3968;
t4215 = (-(t3850 * t4024 - t4021 * t4259) * t4230 + pkin(3) * (t4021 * t3850 + t4024 * t4259)) * t3986;
t3851 = -t3907 * t3957 + t3910 * t3969;
t4214 = (-(t3851 * t4034 - t4031 * t4258) * t4230 + pkin(3) * (t4031 * t3851 + t4034 * t4258)) * t3987;
t3852 = -t3907 * t3958 + t3910 * t3970;
t4213 = (-(t3852 * t4035 - t4032 * t4257) * t4230 + pkin(3) * (t4032 * t3852 + t4035 * t4257)) * t3988;
t3853 = -t3907 * t3959 + t3910 * t3971;
t4212 = (-(t3853 * t4036 - t4033 * t4256) * t4230 + pkin(3) * (t4033 * t3853 + t4036 * t4256)) * t3989;
t4211 = (t3758 * t4002 - t3889 * t4006) / ((-t4003 * t4022 + t4005 * t4125) * t4230 + pkin(3) * (t4003 * t4019 + t4005 * t4122));
t4210 = (t3759 * t4002 - t3890 * t4006) / ((-t4003 * t4023 + t4005 * t4124) * t4230 + pkin(3) * (t4003 * t4020 + t4005 * t4121));
t4209 = (t3760 * t4002 - t3891 * t4006) / ((-t4003 * t4024 + t4005 * t4123) * t4230 + pkin(3) * (t4003 * t4021 + t4005 * t4120));
t4208 = (t3761 * t4002 - t3892 * t4006) / ((-t4003 * t4034 + t4005 * t4119) * t4230 + pkin(3) * (t4003 * t4031 + t4005 * t4116));
t4207 = (t3762 * t4002 - t3893 * t4006) / ((-t4003 * t4035 + t4005 * t4118) * t4230 + pkin(3) * (t4003 * t4032 + t4005 * t4115));
t4206 = (t3763 * t4002 - t3894 * t4006) / ((-t4003 * t4036 + t4005 * t4117) * t4230 + pkin(3) * (t4003 * t4033 + t4005 * t4114));
t4199 = t3866 * t3889;
t4198 = t3867 * t3890;
t4197 = t3868 * t3891;
t4196 = t3872 * t3892;
t4195 = t3873 * t3893;
t4194 = t3874 * t3894;
t3917 = -t3954 * t3999 + t3966 * t4004;
t4187 = t3917 * t3984;
t4186 = t3917 * t3998;
t4185 = t3917 * t4003;
t3918 = -t3955 * t3999 + t3967 * t4004;
t4184 = t3918 * t3985;
t4183 = t3918 * t3998;
t4182 = t3918 * t4003;
t3919 = -t3956 * t3999 + t3968 * t4004;
t4181 = t3919 * t3986;
t4180 = t3919 * t3998;
t4179 = t3919 * t4003;
t3920 = -t3957 * t3999 + t3969 * t4004;
t4178 = t3920 * t3987;
t4177 = t3920 * t3998;
t4176 = t3920 * t4003;
t3921 = -t3958 * t3999 + t3970 * t4004;
t4175 = t3921 * t3988;
t4174 = t3921 * t3998;
t4173 = t3921 * t4003;
t3922 = -t3959 * t3999 + t3971 * t4004;
t4172 = t3922 * t3989;
t4171 = t3922 * t3998;
t4170 = t3922 * t4003;
t4157 = t3936 * t4005;
t4156 = t3937 * t4005;
t4155 = t3938 * t4005;
t4154 = t3942 * t4005;
t4153 = t3943 * t4005;
t4152 = t3944 * t4005;
t4138 = t3991 * t4042;
t4137 = t3991 * t4043;
t4130 = t3994 * t4041;
t3996 = 0.1e1 / t4002;
t3997 = m(1) + t4037;
t4129 = t3996 * t3997;
t4128 = t3996 * t4037;
t4094 = m(3) * t4223;
t4093 = m(3) * t4222;
t4092 = m(3) * t4221;
t4091 = m(3) * t4220;
t4090 = m(3) * t4219;
t4089 = m(3) * t4218;
t4080 = t4128 * t4211;
t4079 = t4128 * t4210;
t4078 = t4128 * t4209;
t4077 = t4128 * t4208;
t4076 = t4128 * t4207;
t4075 = t4128 * t4206;
t4074 = t4129 * t4199;
t4073 = t4129 * t4198;
t4072 = t4129 * t4197;
t4071 = t4129 * t4196;
t4070 = t4129 * t4195;
t4069 = t4129 * t4194;
t3911 = t3954 * t4004 + t3966 * t3999;
t3912 = t3955 * t4004 + t3967 * t3999;
t3913 = t3956 * t4004 + t3968 * t3999;
t3914 = t3957 * t4004 + t3969 * t3999;
t3915 = t3958 * t4004 + t3970 * t3999;
t3916 = t3959 * t4004 + t3971 * t3999;
t3888 = t3991 * t4049 - t3994 * t4243;
t3887 = t3991 * t4048 - t3994 * t4244;
t3886 = t3991 * t4047 - t3994 * t4245;
t3885 = t3991 * t4046 - t3994 * t4246;
t3884 = t3991 * t4045 - t3994 * t4247;
t3883 = t3991 * t4044 - t3994 * t4248;
t3865 = t3916 * t3983 * t4006 + t4002 * t3989;
t3864 = t3915 * t3982 * t4006 + t4002 * t3988;
t3863 = t3914 * t3981 * t4006 + t4002 * t3987;
t3862 = t3913 * t3980 * t4006 + t3986 * t4002;
t3861 = t3912 * t3979 * t4006 + t4002 * t3985;
t3860 = t3911 * t3978 * t4006 + t4002 * t3984;
t3847 = t3907 * t3971 + t3910 * t3959;
t3846 = t3907 * t3970 + t3910 * t3958;
t3845 = t3907 * t3969 + t3910 * t3957;
t3844 = t3907 * t3968 + t3910 * t3956;
t3843 = t3907 * t3967 + t3910 * t3955;
t3842 = t3907 * t3966 + t3910 * t3954;
t3841 = -t3916 * t3977 - t3922 * t4146;
t3840 = -t3915 * t3976 - t3921 * t4147;
t3839 = -t3914 * t3975 - t3920 * t4148;
t3838 = -t3913 * t3974 - t3919 * t4149;
t3837 = -t3912 * t3973 - t3918 * t4150;
t3836 = -t3911 * t3972 - t3917 * t4151;
t3835 = -t3916 * t3965 + t3922 * t4140;
t3834 = -t3915 * t3964 + t3921 * t4141;
t3833 = -t3914 * t3963 + t3920 * t4142;
t3832 = -t3913 * t3962 + t3919 * t4143;
t3831 = -t3912 * t3961 + t3918 * t4144;
t3830 = -t3911 * t3960 + t3917 * t4145;
t3781 = t3847 * t3977 + t3853 * t4146;
t3780 = t3846 * t3976 + t3852 * t4147;
t3779 = t3845 * t3975 + t3851 * t4148;
t3778 = t3844 * t3974 + t3850 * t4149;
t3777 = t3843 * t3973 + t3849 * t4150;
t3776 = t3842 * t3972 + t3848 * t4151;
t3775 = -t3847 * t3965 + t3853 * t4140;
t3774 = -t3846 * t3964 + t3852 * t4141;
t3773 = -t3845 * t3963 + t3851 * t4142;
t3772 = -t3844 * t3962 + t3850 * t4143;
t3771 = -t3843 * t3961 + t3849 * t4144;
t3770 = -t3842 * t3960 + t3848 * t4145;
t3769 = (t3916 * t4003 + t3922 * t4126) * t4152 + (-t3916 * t3998 + t3922 * t4102) * t3947;
t3768 = (t3915 * t4003 + t3921 * t4126) * t4153 + (-t3915 * t3998 + t3921 * t4102) * t3946;
t3767 = (t3914 * t4003 + t3920 * t4126) * t4154 + (-t3914 * t3998 + t3920 * t4102) * t3945;
t3766 = (t3913 * t4003 + t3919 * t4126) * t4155 + (-t3913 * t3998 + t3919 * t4102) * t3941;
t3765 = (t3912 * t4003 + t3918 * t4126) * t4156 + (-t3912 * t3998 + t3918 * t4102) * t3940;
t3764 = (t3911 * t4003 + t3917 * t4126) * t4157 + (-t3911 * t3998 + t3917 * t4102) * t3939;
t3757 = -t3787 * t3965 + t3793 * t4140;
t3756 = -t3786 * t3964 + t3792 * t4141;
t3755 = -t3785 * t3963 + t3791 * t4142;
t3754 = -t3784 * t3962 + t3790 * t4143;
t3753 = -t3783 * t3961 + t3789 * t4144;
t3752 = -t3782 * t3960 + t3788 * t4145;
t3751 = (t3865 * t3998 - t3983 * t4170) * t4152 + t3947 * (t3865 * t4003 + t3983 * t4171);
t3750 = (t3864 * t3998 - t3982 * t4173) * t4153 + t3946 * (t3864 * t4003 + t3982 * t4174);
t3749 = (t3863 * t3998 - t3981 * t4176) * t4154 + t3945 * (t3863 * t4003 + t3981 * t4177);
t3748 = (t3862 * t3998 - t3980 * t4179) * t4155 + t3941 * (t3862 * t4003 + t3980 * t4180);
t3747 = (t3861 * t3998 - t3979 * t4182) * t4156 + t3940 * (t3861 * t4003 + t3979 * t4183);
t3746 = (t3860 * t3998 - t3978 * t4185) * t4157 + t3939 * (t3860 * t4003 + t3978 * t4186);
t3745 = (-(t3916 * t4126 - t4170) * t4152 - t3947 * (t3916 * t4102 + t4171)) * t3989 + t3983 * t3877 * t4002;
t3744 = (-(t3915 * t4126 - t4173) * t4153 - t3946 * (t3915 * t4102 + t4174)) * t3988 + t3982 * t3876 * t4002;
t3743 = (-(t3914 * t4126 - t4176) * t4154 - t3945 * (t3914 * t4102 + t4177)) * t3987 + t3981 * t3875 * t4002;
t3736 = (-(t3913 * t4126 - t4179) * t4155 - t3941 * (t3913 * t4102 + t4180)) * t3986 + t3980 * t3871 * t4002;
t3735 = (-(t3912 * t4126 - t4182) * t4156 - t3940 * (t3912 * t4102 + t4183)) * t3985 + t3979 * t3870 * t4002;
t3734 = (-(t3911 * t4126 - t4185) * t4157 - t3939 * (t3911 * t4102 + t4186)) * t3984 + t3978 * t3869 * t4002;
t3721 = -t3751 * t3965 + t3769 * t3977;
t3720 = t3751 * t3977 + t3769 * t3965;
t3719 = -t3750 * t3964 + t3768 * t3976;
t3718 = t3750 * t3976 + t3768 * t3964;
t3717 = -t3749 * t3963 + t3767 * t3975;
t3716 = t3749 * t3975 + t3767 * t3963;
t3715 = -t3748 * t3962 + t3766 * t3974;
t3714 = t3748 * t3974 + t3766 * t3962;
t3713 = -t3747 * t3961 + t3765 * t3973;
t3712 = t3747 * t3973 + t3765 * t3961;
t3711 = -t3746 * t3960 + t3764 * t3972;
t3710 = t3746 * t3972 + t3764 * t3960;
t3709 = -(t3781 * t4036 + t4033 * t4254) * t4230 + pkin(3) * (t4033 * t3781 - t4036 * t4254);
t3708 = -(t3780 * t4035 + t4032 * t4253) * t4230 + pkin(3) * (t4032 * t3780 - t4035 * t4253);
t3707 = -(t3779 * t4034 + t4031 * t4252) * t4230 + pkin(3) * (t4031 * t3779 - t4034 * t4252);
t3706 = (t3757 * t4033 + t3775 * t4036) * t4230 - (-t3757 * t4036 + t3775 * t4033) * pkin(3);
t3705 = (t3756 * t4032 + t3774 * t4035) * t4230 - (-t3756 * t4035 + t3774 * t4032) * pkin(3);
t3704 = (t3755 * t4031 + t3773 * t4034) * t4230 - (-t3755 * t4034 + t3773 * t4031) * pkin(3);
t3703 = -(t3778 * t4024 + t4021 * t4251) * t4230 + pkin(3) * (t4021 * t3778 - t4024 * t4251);
t3702 = -(t3777 * t4023 + t4020 * t4250) * t4230 + pkin(3) * (t4020 * t3777 - t4023 * t4250);
t3701 = -(t3776 * t4022 + t4019 * t4249) * t4230 + pkin(3) * (t4019 * t3776 - t4022 * t4249);
t3700 = (t3754 * t4021 + t3772 * t4024) * t4230 - (-t3754 * t4024 + t3772 * t4021) * pkin(3);
t3699 = (t3753 * t4020 + t3771 * t4023) * t4230 - (-t3753 * t4023 + t3771 * t4020) * pkin(3);
t3698 = (t3752 * t4019 + t3770 * t4022) * t4230 - (-t3752 * t4022 + t3770 * t4019) * pkin(3);
t1 = [-m(4) * g(1) + (t4172 * t4218 + t4175 * t4219 + t4178 * t4220 + t4181 * t4221 + t4184 * t4222 + t4187 * t4223) * m(3) + ((t4206 * t4212 + t4207 * t4213 + t4208 * t4214 + t4209 * t4215 + t4210 * t4216 + t4211 * t4217) * t4037 + (-t3734 * t4199 - t3735 * t4198 - t3736 * t4197 - t3743 * t4196 - t3744 * t4195 - t3745 * t4194) * t3997) * t3996; -m(4) * g(2) + (-t3836 * t4223 - t3837 * t4222 - t3838 * t4221 - t3839 * t4220 - t3840 * t4219 - t3841 * t4218) * m(3) + ((t3701 * t4211 + t3702 * t4210 + t3703 * t4209 + t3707 * t4208 + t3708 * t4207 + t3709 * t4206) * t4037 + (-t3711 * t4199 - t3713 * t4198 - t3715 * t4197 - t3717 * t4196 - t3719 * t4195 - t3721 * t4194) * t3997) * t3996; -m(4) * g(3) + (-t3830 * t4223 - t3831 * t4222 - t3832 * t4221 - t3833 * t4220 - t3834 * t4219 - t3835 * t4218) * m(3) + ((t3698 * t4211 + t3699 * t4210 + t3700 * t4209 + t3704 * t4208 + t3705 * t4207 + t3706 * t4206) * t4037 + (-t3710 * t4199 - t3712 * t4198 - t3714 * t4197 - t3716 * t4196 - t3718 * t4195 - t3720 * t4194) * t3997) * t3996; -(-t3720 * t3826 - t3721 * t3827) * t4069 + (-t3706 * t3826 - t3709 * t3827) * t4075 - (-t3826 * t3835 - t3827 * t3841) * t4089 - (-t3718 * t3822 - t3719 * t3823) * t4070 + (-t3705 * t3822 - t3708 * t3823) * t4076 - (-t3822 * t3834 - t3823 * t3840) * t4090 - (-t3716 * t3818 - t3717 * t3819) * t4071 + (-t3704 * t3818 - t3707 * t3819) * t4077 - (-t3818 * t3833 - t3819 * t3839) * t4091 - (-t3714 * t3814 - t3715 * t3815) * t4072 + (-t3700 * t3814 - t3703 * t3815) * t4078 - (-t3814 * t3832 - t3815 * t3838) * t4092 - (-t3712 * t3810 - t3713 * t3811) * t4073 + (-t3699 * t3810 - t3702 * t3811) * t4079 - (-t3810 * t3831 - t3811 * t3837) * t4093 - (-t3710 * t3806 - t3711 * t3807) * t4074 + (-t3698 * t3806 - t3701 * t3807) * t4080 - (-t3806 * t3830 - t3807 * t3836) * t4094 + (((-g(2) * t4137 - g(3) * t4042) * t3993 + g(2) * t4130) * t3995 + ((g(2) * t4042 - g(3) * t4137) * t3993 + g(3) * t4130) * t3992 + ((g(2) * t4138 - g(3) * t4043) * t3995 + (g(2) * t4043 + g(3) * t4138) * t3992) * t3990) * m(4); -(-t3720 * t3888 + t3745 * t3827) * t4069 + (-t3706 * t3888 + t3827 * t4212) * t4075 - (-t3827 * t4172 - t3835 * t3888) * t4089 - (-t3718 * t3887 + t3744 * t3823) * t4070 + (-t3705 * t3887 + t3823 * t4213) * t4076 - (-t3823 * t4175 - t3834 * t3887) * t4090 - (-t3716 * t3886 + t3743 * t3819) * t4071 + (-t3704 * t3886 + t3819 * t4214) * t4077 - (-t3819 * t4178 - t3833 * t3886) * t4091 - (-t3714 * t3885 + t3736 * t3815) * t4072 + (-t3700 * t3885 + t3815 * t4215) * t4078 - (-t3815 * t4181 - t3832 * t3885) * t4092 - (-t3712 * t3884 + t3735 * t3811) * t4073 + (-t3699 * t3884 + t3811 * t4216) * t4079 - (-t3811 * t4184 - t3831 * t3884) * t4093 - (-t3710 * t3883 + t3734 * t3807) * t4074 + (-t3698 * t3883 + t3807 * t4217) * t4080 - (-t3807 * t4187 - t3830 * t3883) * t4094 - (t4255 * g(3) + ((t3991 * t4068 + t4130) * t3995 + (t3990 * t4043 + t3993 * t4042) * t3992) * g(1)) * m(4); -(t3721 * t3888 + t3745 * t3826) * t4069 + (t3709 * t3888 + t3826 * t4212) * t4075 - (-t3826 * t4172 + t3841 * t3888) * t4089 - (t3719 * t3887 + t3744 * t3822) * t4070 + (t3708 * t3887 + t3822 * t4213) * t4076 - (-t3822 * t4175 + t3840 * t3887) * t4090 - (t3717 * t3886 + t3743 * t3818) * t4071 + (t3707 * t3886 + t3818 * t4214) * t4077 - (-t3818 * t4178 + t3839 * t3886) * t4091 - (t3715 * t3885 + t3736 * t3814) * t4072 + (t3703 * t3885 + t3814 * t4215) * t4078 - (-t3814 * t4181 + t3838 * t3885) * t4092 - (t3713 * t3884 + t3735 * t3810) * t4073 + (t3702 * t3884 + t3810 * t4216) * t4079 - (-t3810 * t4184 + t3837 * t3884) * t4093 - (t3711 * t3883 + t3734 * t3806) * t4074 + (t3701 * t3883 + t3806 * t4217) * t4080 - (-t3806 * t4187 + t3836 * t3883) * t4094 + (t4255 * g(2) + (-t3992 * t4130 + (t3992 * t4137 + t3995 * t4042) * t3993 + (-t3992 * t4138 + t3995 * t4043) * t3990) * g(1)) * m(4);];
taugX  = t1;
