% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RPRRR9V1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RPRRR9V1G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RPRRR9V1G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RPRRR9V1G3A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:58:19
% EndTime: 2020-08-06 18:58:22
% DurationCPUTime: 2.98s
% Computational Cost: add. (873->185), mult. (1107->324), div. (153->10), fcn. (1101->26), ass. (0->157)
t3993 = sin(pkin(7));
t3994 = cos(pkin(7));
t4114 = MDP(4) * t3994 - MDP(5) * t3993 + MDP(2);
t3995 = legFrame(3,2);
t3971 = sin(t3995);
t3974 = cos(t3995);
t3942 = t3974 * g(1) - t3971 * g(2);
t3999 = sin(qJ(1,3));
t4005 = cos(qJ(1,3));
t3930 = g(3) * t4005 + t3942 * t3999;
t4100 = pkin(5) + pkin(6);
t3983 = qJ(2,3) + t4100;
t3968 = 0.1e1 / t3983;
t4056 = t3999 * t3968;
t4110 = t3930 * t4056;
t3987 = pkin(7) + qJ(3,3);
t3964 = cos(t3987);
t3958 = 0.1e1 / t3964;
t4073 = t3958 * t3968;
t4109 = t3930 * t4073;
t3996 = legFrame(2,2);
t3972 = sin(t3996);
t3975 = cos(t3996);
t3943 = t3975 * g(1) - t3972 * g(2);
t4001 = sin(qJ(1,2));
t4007 = cos(qJ(1,2));
t3932 = g(3) * t4007 + t3943 * t4001;
t3984 = qJ(2,2) + t4100;
t3969 = 0.1e1 / t3984;
t4055 = t4001 * t3969;
t4108 = t3932 * t4055;
t3988 = pkin(7) + qJ(3,2);
t3965 = cos(t3988);
t3959 = 0.1e1 / t3965;
t4070 = t3959 * t3969;
t4107 = t3932 * t4070;
t3985 = qJ(2,1) + t4100;
t3970 = 0.1e1 / t3985;
t4003 = sin(qJ(1,1));
t4053 = t4003 * t3970;
t3997 = legFrame(1,2);
t3973 = sin(t3997);
t3976 = cos(t3997);
t3944 = t3976 * g(1) - t3973 * g(2);
t4009 = cos(qJ(1,1));
t4103 = g(3) * t4009 + t3944 * t4003;
t4041 = t4103 * t4053;
t3989 = pkin(7) + qJ(3,1);
t3966 = cos(t3989);
t3960 = 0.1e1 / t3966;
t4067 = t3960 * t3970;
t4105 = t4103 * t4067;
t4104 = MDP(3) - MDP(6);
t4102 = 2 * pkin(3);
t4101 = 0.2e1 * t3994 ^ 2;
t4004 = cos(qJ(3,3));
t3990 = t4004 ^ 2;
t4096 = t3990 * pkin(3);
t4006 = cos(qJ(3,2));
t3991 = t4006 ^ 2;
t4095 = t3991 * pkin(3);
t4008 = cos(qJ(3,1));
t3992 = t4008 ^ 2;
t4094 = t3992 * pkin(3);
t4093 = t4004 * pkin(2);
t4092 = t4006 * pkin(2);
t4091 = t4008 * pkin(2);
t3981 = pkin(1) * t4005;
t3906 = t3942 * (t3999 * pkin(1) - t4005 * qJ(2,3)) + g(3) * (t3999 * qJ(2,3) + t3981);
t4088 = t3906 * t3958;
t3982 = pkin(1) * t4007;
t3907 = t3943 * (t4001 * pkin(1) - t4007 * qJ(2,2)) + g(3) * (t4001 * qJ(2,2) + t3982);
t4087 = t3907 * t3959;
t3977 = t4009 * pkin(1);
t3908 = t3944 * (t4003 * pkin(1) - t4009 * qJ(2,1)) + g(3) * (t4003 * qJ(2,1) + t3977);
t4086 = t3908 * t3960;
t4085 = t4103 * t3970;
t3998 = sin(qJ(3,3));
t4058 = t3993 * t3998;
t4084 = t3930 / (t3994 * t4004 - t4058);
t4083 = t3930 * t3968;
t4000 = sin(qJ(3,2));
t4057 = t3993 * t4000;
t4082 = t3932 / (t3994 * t4006 - t4057);
t4081 = t3932 * t3969;
t4002 = sin(qJ(3,1));
t4054 = t4002 * t3993;
t4080 = t4103 / (t4008 * t3994 - t4054);
t4011 = pkin(2) / 0.2e1;
t4076 = (t4004 * pkin(3) + t4011) * t3998;
t4075 = (t4006 * pkin(3) + t4011) * t4000;
t4074 = (t4008 * pkin(3) + t4011) * t4002;
t4072 = t3958 * t3971;
t4071 = t3958 * t3974;
t4069 = t3959 * t3972;
t4068 = t3959 * t3975;
t4066 = t3960 * t3973;
t4065 = t3960 * t3976;
t4064 = t3971 * t4005;
t4063 = t3972 * t4007;
t4062 = t3973 * t4009;
t4061 = t3974 * t4005;
t4060 = t3975 * t4007;
t4059 = t3976 * t4009;
t3967 = pkin(1) * t3993;
t4052 = t4004 * (-t3998 * pkin(3) + t3967);
t4051 = t4006 * (-t4000 * pkin(3) + t3967);
t4050 = t4008 * (-t4002 * pkin(3) + t3967);
t4049 = t4003 * t3985 + t3977;
t4048 = t3999 * t3983 + t3981;
t4047 = t4001 * t3984 + t3982;
t4029 = -g(3) * t3999 + t3942 * t4005;
t4039 = t4029 * t4073;
t4028 = -g(3) * t4001 + t3943 * t4007;
t4037 = t4028 * t4070;
t4027 = -g(3) * t4003 + t3944 * t4009;
t4035 = t4027 * t4067;
t4034 = t4005 * t4058;
t4033 = t4007 * t4057;
t4032 = t4009 * t4054;
t3961 = sin(t3987);
t4024 = t3961 * t4109;
t3962 = sin(t3988);
t4023 = t3962 * t4107;
t3963 = sin(t3989);
t4022 = t3963 * t4105;
t4021 = pkin(2) * t4034 + (t4034 * t4102 - t4048) * t4004;
t4020 = pkin(2) * t4033 + (t4033 * t4102 - t4047) * t4006;
t4019 = pkin(2) * t4032 + (t4032 * t4102 - t4049) * t4008;
t4012 = 1 / pkin(3);
t4010 = -pkin(3) / 0.2e1;
t3957 = t3994 * pkin(2) + pkin(1);
t3947 = t4094 + t4091 / 0.2e1 + t4010;
t3946 = t4095 + t4092 / 0.2e1 + t4010;
t3945 = t4096 + t4093 / 0.2e1 + t4010;
t3941 = t3973 * g(1) + t3976 * g(2);
t3940 = t3972 * g(1) + t3975 * g(2);
t3939 = t3971 * g(1) + t3974 * g(2);
t3920 = t3973 * t3963 + t3966 * t4059;
t3919 = t3976 * t3963 - t3966 * t4062;
t3918 = t3972 * t3962 + t3965 * t4060;
t3917 = t3975 * t3962 - t3965 * t4063;
t3916 = t3971 * t3961 + t3964 * t4061;
t3915 = t3974 * t3961 - t3964 * t4064;
t3914 = pkin(1) * t4002 + (-pkin(3) + t4091 + 0.2e1 * t4094) * t3993;
t3913 = pkin(1) * t4000 + (-pkin(3) + t4092 + 0.2e1 * t4095) * t3993;
t3912 = pkin(1) * t3998 + (-pkin(3) + t4093 + 0.2e1 * t4096) * t3993;
t3911 = t4049 * t4054 + (t3992 - 0.1e1) * t4009 * pkin(3);
t3910 = t4047 * t4057 + (t3991 - 0.1e1) * t4007 * pkin(3);
t3909 = t4048 * t4058 + (t3990 - 0.1e1) * t4005 * pkin(3);
t3905 = t3941 * t3963 + t3966 * t4027;
t3904 = -t3941 * t3966 + t3963 * t4027;
t3903 = t3940 * t3962 + t3965 * t4028;
t3902 = -t3940 * t3965 + t3962 * t4028;
t3901 = t3939 * t3961 + t3964 * t4029;
t3900 = -t3939 * t3964 + t3961 * t4029;
t1 = [((t3920 * t4086 - ((t3947 * t4059 + t3973 * t4074) * t4101 + (t3973 * t3914 - t4019 * t3976) * t3994 - t3911 * t3976 + t3973 * t4050) * t4080) * t3970 + (t3918 * t4087 - ((t3946 * t4060 + t3972 * t4075) * t4101 + (t3972 * t3913 - t4020 * t3975) * t3994 - t3910 * t3975 + t3972 * t4051) * t4082) * t3969 + (t3916 * t4088 - ((t3945 * t4061 + t3971 * t4076) * t4101 + (t3971 * t3912 - t4021 * t3974) * t3994 - t3909 * t3974 + t3971 * t4052) * t4084) * t3968) * MDP(7) + (t3916 * t4083 + t3918 * t4081 + t3920 * t4085 + (t3900 * t4072 + t3902 * t4069 + t3904 * t4066) * t4012) * MDP(13) + (-t3916 * t4024 - t3918 * t4023 - t3920 * t4022 + (t3901 * t4072 + t3903 * t4069 + t3905 * t4066) * t4012) * MDP(14) - g(1) * MDP(15) + t4104 * (t3916 * t4039 + t3918 * t4037 + t3920 * t4035) + t4114 * (t3916 * t4109 + t3918 * t4107 + t3920 * t4105); ((t3919 * t4086 - ((-t3947 * t4062 + t3976 * t4074) * t4101 + (t3976 * t3914 + t4019 * t3973) * t3994 + t3911 * t3973 + t3976 * t4050) * t4080) * t3970 + (t3917 * t4087 - ((-t3946 * t4063 + t3975 * t4075) * t4101 + (t3975 * t3913 + t4020 * t3972) * t3994 + t3910 * t3972 + t3975 * t4051) * t4082) * t3969 + (t3915 * t4088 - ((-t3945 * t4064 + t3974 * t4076) * t4101 + (t3974 * t3912 + t4021 * t3971) * t3994 + t3909 * t3971 + t3974 * t4052) * t4084) * t3968) * MDP(7) + (t3915 * t4083 + t3917 * t4081 + t3919 * t4085 + (t3900 * t4071 + t3902 * t4068 + t3904 * t4065) * t4012) * MDP(13) + (-t3915 * t4024 - t3917 * t4023 - t3919 * t4022 + (t3901 * t4071 + t3903 * t4068 + t3905 * t4065) * t4012) * MDP(14) - g(2) * MDP(15) + t4104 * (t3915 * t4039 + t3917 * t4037 + t3919 * t4035) + t4114 * (t3915 * t4109 + t3917 * t4107 + t3919 * t4105); ((-t3985 * t4009 * t4103 + (-t3908 - (-pkin(3) * t3966 - t3957) * t4103) * t4003) * t3970 + (-t3984 * t4007 * t3932 + (-t3907 - (-pkin(3) * t3965 - t3957) * t3932) * t4001) * t3969 + (-t3983 * t4005 * t3930 + (-t3906 - (-pkin(3) * t3964 - t3957) * t3930) * t3999) * t3968) * MDP(7) + (-t3964 * t4110 - t3965 * t4108 - t3966 * t4041) * MDP(13) + (t3961 * t4110 + t3962 * t4108 + t3963 * t4041) * MDP(14) - g(3) * MDP(15) - t4104 * (t4027 * t4053 + t4028 * t4055 + t4029 * t4056) + t4114 * (-t4108 - t4110 - t4041);];
taugX  = t1;
