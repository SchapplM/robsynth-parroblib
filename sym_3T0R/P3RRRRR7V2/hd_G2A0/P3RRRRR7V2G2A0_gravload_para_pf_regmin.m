% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x18]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:08
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR7V2G2A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G2A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 09:57:40
% EndTime: 2020-08-07 09:57:43
% DurationCPUTime: 3.41s
% Computational Cost: add. (2181->345), mult. (3690->501), div. (210->14), fcn. (2796->66), ass. (0->260)
t4097 = legFrame(1,2);
t4061 = sin(t4097);
t4064 = cos(t4097);
t4031 = t4064 * g(1) - t4061 * g(2);
t4106 = sin(qJ(1,1));
t4115 = cos(qJ(1,1));
t4007 = -g(3) * t4106 + t4031 * t4115;
t4092 = qJ(2,1) + qJ(3,1);
t4052 = cos(t4092);
t4114 = cos(qJ(2,1));
t4025 = 0.1e1 / (t4114 * pkin(2) + pkin(3) * t4052 + pkin(1));
t4221 = t4025 * t4115;
t4158 = t4064 * t4221;
t4141 = t4007 * t4158;
t4159 = t4061 * t4221;
t4142 = t4007 * t4159;
t4222 = t4025 * t4106;
t4165 = t4007 * t4222;
t4095 = legFrame(3,2);
t4059 = sin(t4095);
t4062 = cos(t4095);
t4029 = t4062 * g(1) - t4059 * g(2);
t4100 = sin(qJ(1,3));
t4109 = cos(qJ(1,3));
t4011 = -g(3) * t4100 + t4029 * t4109;
t4090 = qJ(2,3) + qJ(3,3);
t4050 = cos(t4090);
t4108 = cos(qJ(2,3));
t4023 = 0.1e1 / (t4108 * pkin(2) + pkin(3) * t4050 + pkin(1));
t4226 = t4023 * t4100;
t4271 = t4011 * t4226;
t4096 = legFrame(2,2);
t4060 = sin(t4096);
t4063 = cos(t4096);
t4030 = t4063 * g(1) - t4060 * g(2);
t4103 = sin(qJ(1,2));
t4112 = cos(qJ(1,2));
t4012 = -g(3) * t4103 + t4030 * t4112;
t4091 = qJ(2,2) + qJ(3,2);
t4051 = cos(t4091);
t4111 = cos(qJ(2,2));
t4024 = 0.1e1 / (t4111 * pkin(2) + pkin(3) * t4051 + pkin(1));
t4224 = t4024 * t4103;
t4270 = t4012 * t4224;
t4049 = sin(t4092);
t4250 = qJ(3,1) + qJ(1,1);
t4057 = qJ(2,1) + t4250;
t4251 = -qJ(3,1) + qJ(1,1);
t4058 = -qJ(2,1) + t4251;
t4105 = sin(qJ(2,1));
t4125 = 0.2e1 * qJ(3,1);
t4130 = pkin(2) ^ 2;
t4239 = 0.2e1 * pkin(1);
t4244 = qJ(1,1) - 0.2e1 * qJ(2,1);
t4245 = qJ(1,1) + 0.2e1 * qJ(2,1);
t4104 = sin(qJ(3,1));
t4256 = t4104 * pkin(1);
t4080 = (pkin(5) + pkin(6) + pkin(7));
t4265 = 2 * t4080;
t4269 = ((cos(t4058) + cos(t4057)) * t4239 + (sin(t4058) + sin(t4057)) * t4265 + (cos(0.2e1 * qJ(3,1) - t4244) + cos(t4125 + t4245) + 0.2e1 * t4115) * pkin(3) + (cos(qJ(3,1) - t4244) + cos(qJ(3,1) + t4245) + cos(t4251) + cos(t4250)) * pkin(2)) / (-t4130 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t4256 + pkin(2) * t4049 + (sin(t4125 + qJ(2,1)) - t4105) * pkin(3))) / 0.2e1;
t4048 = sin(t4091);
t4248 = qJ(3,2) + qJ(1,2);
t4055 = qJ(2,2) + t4248;
t4249 = -qJ(3,2) + qJ(1,2);
t4056 = -qJ(2,2) + t4249;
t4102 = sin(qJ(2,2));
t4122 = 0.2e1 * qJ(3,2);
t4242 = qJ(1,2) - 0.2e1 * qJ(2,2);
t4243 = qJ(1,2) + 0.2e1 * qJ(2,2);
t4101 = sin(qJ(3,2));
t4258 = t4101 * pkin(1);
t4268 = ((cos(t4056) + cos(t4055)) * t4239 + (sin(t4056) + sin(t4055)) * t4265 + (cos(0.2e1 * qJ(3,2) - t4242) + cos(t4122 + t4243) + 0.2e1 * t4112) * pkin(3) + (cos(qJ(3,2) - t4242) + cos(qJ(3,2) + t4243) + cos(t4249) + cos(t4248)) * pkin(2)) / (-t4130 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t4258 + pkin(2) * t4048 + (sin(t4122 + qJ(2,2)) - t4102) * pkin(3))) / 0.2e1;
t4047 = sin(t4090);
t4246 = qJ(3,3) + qJ(1,3);
t4053 = qJ(2,3) + t4246;
t4247 = -qJ(3,3) + qJ(1,3);
t4054 = -qJ(2,3) + t4247;
t4099 = sin(qJ(2,3));
t4119 = 0.2e1 * qJ(3,3);
t4240 = -0.2e1 * qJ(2,3) + qJ(1,3);
t4241 = 0.2e1 * qJ(2,3) + qJ(1,3);
t4098 = sin(qJ(3,3));
t4260 = t4098 * pkin(1);
t4267 = ((cos(t4054) + cos(t4053)) * t4239 + (sin(t4054) + sin(t4053)) * t4265 + (cos(0.2e1 * qJ(3,3) - t4240) + cos(t4119 + t4241) + 0.2e1 * t4109) * pkin(3) + (cos(qJ(3,3) - t4240) + cos(qJ(3,3) + t4241) + cos(t4247) + cos(t4246)) * pkin(2)) / (-t4130 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t4260 + pkin(2) * t4047 + (sin(t4119 + qJ(2,3)) - t4099) * pkin(3))) / 0.2e1;
t4266 = 0.2e1 * pkin(3);
t4085 = t4108 ^ 2;
t4264 = 0.2e1 * t4085;
t4087 = t4111 ^ 2;
t4263 = 0.2e1 * t4087;
t4089 = t4114 ^ 2;
t4262 = 0.2e1 * t4089;
t4107 = cos(qJ(3,3));
t4084 = t4107 ^ 2;
t4073 = pkin(3) * t4084;
t4110 = cos(qJ(3,2));
t4086 = t4110 ^ 2;
t4074 = pkin(3) * t4086;
t4113 = cos(qJ(3,1));
t4088 = t4113 ^ 2;
t4075 = pkin(3) * t4088;
t4259 = t4098 * pkin(3);
t4257 = t4101 * pkin(3);
t4255 = t4104 * pkin(3);
t4254 = t4107 * pkin(2);
t4065 = t4107 * pkin(3);
t4253 = t4110 * pkin(2);
t4066 = t4110 * pkin(3);
t4252 = t4113 * pkin(2);
t4067 = t4113 * pkin(3);
t4149 = pkin(1) * t4100 - t4109 * t4080;
t4198 = t4098 * t4099;
t3993 = t4149 * t4198 + (t4084 - 0.1e1) * t4100 * pkin(3);
t4235 = t3993 * t4059;
t4234 = t3993 * t4062;
t4148 = pkin(1) * t4103 - t4112 * t4080;
t4197 = t4101 * t4102;
t3994 = t4148 * t4197 + (t4086 - 0.1e1) * t4103 * pkin(3);
t4233 = t3994 * t4060;
t4232 = t3994 * t4063;
t4147 = pkin(1) * t4106 - t4115 * t4080;
t4196 = t4104 * t4105;
t3995 = t4147 * t4196 + (t4088 - 0.1e1) * t4106 * pkin(3);
t4231 = t3995 * t4061;
t4230 = t3995 * t4064;
t4184 = pkin(3) * t4198;
t4044 = t4065 + pkin(2);
t4214 = t4044 * t4108;
t4229 = 0.1e1 / (pkin(1) - t4184 + t4214) / t4098;
t4183 = pkin(3) * t4197;
t4045 = t4066 + pkin(2);
t4211 = t4045 * t4111;
t4228 = 0.1e1 / (pkin(1) - t4183 + t4211) / t4101;
t4182 = pkin(3) * t4196;
t4046 = t4067 + pkin(2);
t4208 = t4046 * t4114;
t4227 = 0.1e1 / (pkin(1) - t4182 + t4208) / t4104;
t4225 = t4023 * t4109;
t4223 = t4024 * t4112;
t4118 = pkin(2) / 0.2e1;
t4219 = (t4065 + t4118) * t4098;
t4218 = (t4066 + t4118) * t4101;
t4217 = (t4067 + t4118) * t4104;
t4216 = t4044 * t4059;
t4215 = t4044 * t4062;
t4213 = t4045 * t4060;
t4212 = t4045 * t4063;
t4210 = t4046 * t4061;
t4209 = t4046 * t4064;
t4207 = t4059 * t4100;
t4206 = t4060 * t4103;
t4205 = t4061 * t4106;
t4204 = t4062 * t4100;
t4203 = t4063 * t4103;
t4202 = t4064 * t4106;
t4128 = pkin(3) ^ 2;
t4201 = t4084 * t4128;
t4200 = t4086 * t4128;
t4199 = t4088 * t4128;
t4038 = pkin(1) * t4099 - t4259;
t4195 = t4107 * t4038;
t4039 = pkin(1) * t4102 - t4257;
t4194 = t4110 * t4039;
t4040 = pkin(1) * t4105 - t4255;
t4193 = t4113 * t4040;
t4129 = 0.1e1 / pkin(3);
t4131 = 0.1e1 / pkin(2);
t4192 = t4129 * t4131;
t4191 = -t4128 / 0.2e1 + t4130 / 0.2e1;
t4190 = pkin(2) * t4065;
t4189 = pkin(2) * t4066;
t4188 = pkin(2) * t4067;
t4187 = t4044 * t4259;
t4186 = t4045 * t4257;
t4185 = t4046 * t4255;
t4026 = t4059 * g(1) + t4062 * g(2);
t4152 = g(3) * t4109 + t4029 * t4100;
t3978 = -t4026 * t4050 + t4047 * t4152;
t4181 = t3978 * t4229;
t3979 = t4026 * t4047 + t4050 * t4152;
t4180 = t3979 * t4229;
t4027 = t4060 * g(1) + t4063 * g(2);
t4151 = g(3) * t4112 + t4030 * t4103;
t3980 = -t4027 * t4051 + t4048 * t4151;
t4179 = t3980 * t4228;
t3981 = t4027 * t4048 + t4051 * t4151;
t4178 = t3981 * t4228;
t4028 = t4061 * g(1) + t4064 * g(2);
t4150 = g(3) * t4115 + t4031 * t4106;
t3982 = -t4028 * t4052 + t4049 * t4150;
t4177 = t3982 * t4227;
t3983 = t4028 * t4049 + t4052 * t4150;
t4176 = t3983 * t4227;
t3984 = -t4026 * t4108 + t4099 * t4152;
t4175 = t3984 * t4229;
t3985 = t4026 * t4099 + t4108 * t4152;
t4174 = t3985 * t4229;
t3986 = -t4027 * t4111 + t4102 * t4151;
t4173 = t3986 * t4228;
t3987 = t4027 * t4102 + t4111 * t4151;
t4172 = t3987 * t4228;
t3988 = -t4028 * t4114 + t4105 * t4150;
t4171 = t3988 * t4227;
t3989 = t4028 * t4105 + t4114 * t4150;
t4170 = t3989 * t4227;
t4168 = t4011 * t4225;
t4166 = t4012 * t4223;
t4164 = t4007 * t4221;
t4163 = t4059 * t4225;
t4162 = t4062 * t4225;
t4161 = t4060 * t4223;
t4160 = t4063 * t4223;
t4156 = t4100 * t4198;
t4154 = t4103 * t4197;
t4153 = t4106 * t4196;
t4146 = t4050 * t4168;
t4145 = t4108 * t4168;
t4144 = t4051 * t4166;
t4143 = t4111 * t4166;
t4140 = t4105 * t4164;
t4139 = t4114 * t4164;
t4138 = t4011 * t4163;
t4137 = t4012 * t4161;
t4136 = t4011 * t4162;
t4135 = t4012 * t4160;
t4008 = t4156 * t4266 - t4149;
t4134 = pkin(2) * t4156 + t4008 * t4107;
t4009 = t4154 * t4266 - t4148;
t4133 = pkin(2) * t4154 + t4009 * t4110;
t4010 = t4153 * t4266 - t4147;
t4132 = pkin(2) * t4153 + t4010 * t4113;
t4117 = -pkin(3) / 0.2e1;
t4079 = -t4128 + t4130;
t4034 = t4075 + t4252 / 0.2e1 + t4117;
t4033 = t4074 + t4253 / 0.2e1 + t4117;
t4032 = t4073 + t4254 / 0.2e1 + t4117;
t4022 = t4188 + t4191 + t4199;
t4021 = t4189 + t4191 + t4200;
t4020 = t4190 + t4191 + t4201;
t4001 = t4256 + (-pkin(3) + t4252 + 0.2e1 * t4075) * t4105;
t4000 = t4258 + (-pkin(3) + t4253 + 0.2e1 * t4074) * t4102;
t3999 = t4260 + (-pkin(3) + t4254 + 0.2e1 * t4073) * t4099;
t3998 = pkin(1) * t4255 + (t4079 + 0.2e1 * t4188 + 0.2e1 * t4199) * t4105;
t3997 = pkin(1) * t4257 + (t4079 + 0.2e1 * t4189 + 0.2e1 * t4200) * t4102;
t3996 = pkin(1) * t4259 + (t4079 + 0.2e1 * t4190 + 0.2e1 * t4201) * t4099;
t3977 = -0.2e1 * t4022 * t4115 * t4089 - ((pkin(1) - 0.2e1 * t4182) * t4115 + t4106 * t4080) * t4208 + pkin(3) * ((pkin(1) * t4196 - pkin(3) + t4075) * t4115 + t4080 * t4153);
t3976 = -0.2e1 * t4021 * t4112 * t4087 - ((pkin(1) - 0.2e1 * t4183) * t4112 + t4103 * t4080) * t4211 + pkin(3) * ((pkin(1) * t4197 - pkin(3) + t4074) * t4112 + t4080 * t4154);
t3975 = -0.2e1 * t4020 * t4109 * t4085 - ((pkin(1) - 0.2e1 * t4184) * t4109 + t4100 * t4080) * t4214 + pkin(3) * ((pkin(1) * t4198 - pkin(3) + t4073) * t4109 + t4080 * t4156);
t3971 = (t4034 * t4202 + t4061 * t4217) * t4262 + (t4061 * t4001 - t4132 * t4064) * t4114 - t4230 + t4061 * t4193;
t3970 = (-t4034 * t4205 + t4064 * t4217) * t4262 + (t4064 * t4001 + t4132 * t4061) * t4114 + t4231 + t4064 * t4193;
t3969 = (t4033 * t4203 + t4060 * t4218) * t4263 + (t4060 * t4000 - t4133 * t4063) * t4111 - t4232 + t4060 * t4194;
t3968 = (-t4033 * t4206 + t4063 * t4218) * t4263 + (t4063 * t4000 + t4133 * t4060) * t4111 + t4233 + t4063 * t4194;
t3967 = (t4032 * t4204 + t4059 * t4219) * t4264 + (t4059 * t3999 - t4134 * t4062) * t4108 - t4234 + t4059 * t4195;
t3966 = (-t4032 * t4207 + t4062 * t4219) * t4264 + (t4062 * t3999 + t4134 * t4059) * t4108 + t4235 + t4062 * t4195;
t3965 = (-t4022 * t4202 - t4061 * t4185) * t4262 + (-t4061 * t3998 + t4010 * t4209) * t4114 + pkin(3) * t4230 - t4040 * t4210;
t3964 = (t4022 * t4205 - t4064 * t4185) * t4262 + (-t3998 * t4064 - t4010 * t4210) * t4114 - pkin(3) * t4231 - t4040 * t4209;
t3963 = (-t4021 * t4203 - t4060 * t4186) * t4263 + (-t4060 * t3997 + t4009 * t4212) * t4111 + pkin(3) * t4232 - t4039 * t4213;
t3962 = (t4021 * t4206 - t4063 * t4186) * t4263 + (-t3997 * t4063 - t4009 * t4213) * t4111 - pkin(3) * t4233 - t4039 * t4212;
t3961 = (-t4020 * t4204 - t4059 * t4187) * t4264 + (-t4059 * t3996 + t4008 * t4215) * t4108 + pkin(3) * t4234 - t4038 * t4216;
t3960 = (t4020 * t4207 - t4062 * t4187) * t4264 + (-t3996 * t4062 - t4008 * t4216) * t4108 - pkin(3) * t4235 - t4038 * t4215;
t1 = [0, -t4135 - t4136 - t4141, t4150 * t4158 + t4151 * t4160 + t4152 * t4162, 0, 0, 0, 0, 0, -t4062 * t4145 - t4063 * t4143 - t4064 * t4139 + (t3967 * t4175 + t3969 * t4173 + t3971 * t4171) * t4131, t4064 * t4140 + t4099 * t4136 + t4102 * t4135 + (t3967 * t4174 + t3969 * t4172 + t3971 * t4170) * t4131, 0, 0, 0, 0, 0, -t4062 * t4146 - t4063 * t4144 - t4052 * t4141 + (t3967 * t4181 + t3969 * t4179 + t3971 * t4177 + (t3961 * t4181 + t3963 * t4179 + t3965 * t4177) * t4129) * t4131, t4049 * t4141 + t4047 * t4136 + t4048 * t4135 + (t3967 * t4180 + t3969 * t4178 + t3971 * t4176 + (t3961 * t4180 + t3963 * t4178 + t3965 * t4176) * t4129) * t4131, -g(1); 0, t4137 + t4138 + t4142, -t4150 * t4159 - t4151 * t4161 - t4152 * t4163, 0, 0, 0, 0, 0, t4059 * t4145 + t4060 * t4143 + t4061 * t4139 + (t3966 * t4175 + t3968 * t4173 + t3970 * t4171) * t4131, -t4061 * t4140 - t4099 * t4138 - t4102 * t4137 + (t3966 * t4174 + t3968 * t4172 + t3970 * t4170) * t4131, 0, 0, 0, 0, 0, t4059 * t4146 + t4060 * t4144 + t4052 * t4142 + (t3966 * t4181 + t3968 * t4179 + t3970 * t4177 + (t3960 * t4181 + t3962 * t4179 + t3964 * t4177) * t4129) * t4131, -t4049 * t4142 - t4047 * t4138 - t4048 * t4137 + (t3966 * t4180 + t3968 * t4178 + t3970 * t4176 + (t3960 * t4180 + t3962 * t4178 + t3964 * t4176) * t4129) * t4131, -g(2); 0, t4270 + t4271 + t4165, -t4150 * t4222 - t4151 * t4224 - t4152 * t4226, 0, 0, 0, 0, 0, t3984 * t4267 + t3986 * t4268 + t3988 * t4269 + t4108 * t4271 + t4111 * t4270 + t4114 * t4165, t3985 * t4267 + t3987 * t4268 + t3989 * t4269 - t4099 * t4271 - t4102 * t4270 - t4105 * t4165, 0, 0, 0, 0, 0, t4050 * t4271 + t4051 * t4270 + t4052 * t4165 + (t3975 * t4181 + t3976 * t4179 + t3977 * t4177) * t4192 + t3978 * t4267 + t3980 * t4268 + t3982 * t4269, -t4049 * t4165 - t4047 * t4271 - t4048 * t4270 + (t3975 * t4180 + t3976 * t4178 + t3977 * t4176) * t4192 + t3979 * t4267 + t3981 * t4268 + t3983 * t4269, -g(3);];
tau_reg  = t1;
