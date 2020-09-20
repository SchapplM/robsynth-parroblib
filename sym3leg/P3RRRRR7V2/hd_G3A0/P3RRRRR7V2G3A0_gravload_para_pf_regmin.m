% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR7V2G3A0
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
% Datum: 2020-08-07 10:47
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR7V2G3A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:36:03
% EndTime: 2020-08-07 10:36:06
% DurationCPUTime: 2.76s
% Computational Cost: add. (2181->345), mult. (3702->495), div. (210->14), fcn. (2808->66), ass. (0->257)
t4145 = qJ(2,1) + qJ(3,1);
t4104 = sin(t4145);
t4303 = qJ(3,1) + qJ(1,1);
t4112 = qJ(2,1) + t4303;
t4304 = -qJ(3,1) + qJ(1,1);
t4113 = -qJ(2,1) + t4304;
t4158 = sin(qJ(2,1));
t4159 = sin(qJ(1,1));
t4178 = 2 * qJ(3,1);
t4183 = pkin(2) ^ 2;
t4292 = 2 * pkin(1);
t4297 = qJ(1,1) - 0.2e1 * qJ(2,1);
t4298 = qJ(1,1) + 0.2e1 * qJ(2,1);
t4157 = sin(qJ(3,1));
t4309 = t4157 * pkin(1);
t4133 = (pkin(5) + pkin(6) + pkin(7));
t4320 = 2 * t4133;
t4327 = ((-sin(t4113) - sin(t4112)) * t4292 + (cos(t4113) + cos(t4112)) * t4320 + (sin((2 * qJ(3,1)) - t4297) - sin(t4178 + t4298) - 0.2e1 * t4159) * pkin(3) + (sin(qJ(3,1) - t4297) - sin(qJ(3,1) + t4298) - sin(t4304) - sin(t4303)) * pkin(2)) / (-t4183 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t4309 + pkin(2) * t4104 + (sin(t4178 + qJ(2,1)) - t4158) * pkin(3))) / 0.2e1;
t4144 = qJ(2,2) + qJ(3,2);
t4103 = sin(t4144);
t4301 = qJ(3,2) + qJ(1,2);
t4110 = qJ(2,2) + t4301;
t4302 = -qJ(3,2) + qJ(1,2);
t4111 = -qJ(2,2) + t4302;
t4155 = sin(qJ(2,2));
t4156 = sin(qJ(1,2));
t4175 = 2 * qJ(3,2);
t4295 = qJ(1,2) - 0.2e1 * qJ(2,2);
t4296 = qJ(1,2) + 0.2e1 * qJ(2,2);
t4154 = sin(qJ(3,2));
t4311 = t4154 * pkin(1);
t4326 = ((-sin(t4111) - sin(t4110)) * t4292 + (cos(t4111) + cos(t4110)) * t4320 + (sin((2 * qJ(3,2)) - t4295) - sin(t4175 + t4296) - 0.2e1 * t4156) * pkin(3) + (sin(qJ(3,2) - t4295) - sin(qJ(3,2) + t4296) - sin(t4302) - sin(t4301)) * pkin(2)) / (-t4183 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t4311 + pkin(2) * t4103 + (sin(t4175 + qJ(2,2)) - t4155) * pkin(3))) / 0.2e1;
t4143 = qJ(2,3) + qJ(3,3);
t4102 = sin(t4143);
t4299 = qJ(3,3) + qJ(1,3);
t4108 = qJ(2,3) + t4299;
t4300 = -qJ(3,3) + qJ(1,3);
t4109 = -qJ(2,3) + t4300;
t4152 = sin(qJ(2,3));
t4153 = sin(qJ(1,3));
t4172 = 2 * qJ(3,3);
t4293 = -0.2e1 * qJ(2,3) + qJ(1,3);
t4294 = qJ(1,3) + 0.2e1 * qJ(2,3);
t4151 = sin(qJ(3,3));
t4313 = t4151 * pkin(1);
t4325 = ((-sin(t4109) - sin(t4108)) * t4292 + (cos(t4109) + cos(t4108)) * t4320 + (sin((2 * qJ(3,3)) - t4293) - sin(t4172 + t4294) - 0.2e1 * t4153) * pkin(3) + (sin(qJ(3,3) - t4293) - sin(qJ(3,3) + t4294) - sin(t4300) - sin(t4299)) * pkin(2)) / (-t4183 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t4313 + pkin(2) * t4102 + (sin(t4172 + qJ(2,3)) - t4152) * pkin(3))) / 0.2e1;
t4150 = legFrame(1,2);
t4116 = sin(t4150);
t4119 = cos(t4150);
t4083 = t4119 * g(1) - t4116 * g(2);
t4168 = cos(qJ(1,1));
t4200 = g(3) * t4168 + t4083 * t4159;
t4107 = cos(t4145);
t4167 = cos(qJ(2,1));
t4077 = 0.1e1 / (t4167 * pkin(2) + pkin(3) * t4107 + pkin(1));
t4275 = t4077 * t4159;
t4206 = t4119 * t4275;
t4190 = t4200 * t4206;
t4207 = t4116 * t4275;
t4191 = t4200 * t4207;
t4274 = t4077 * t4168;
t4212 = t4200 * t4274;
t4149 = legFrame(2,2);
t4115 = sin(t4149);
t4118 = cos(t4149);
t4082 = t4118 * g(1) - t4115 * g(2);
t4165 = cos(qJ(1,2));
t4201 = g(3) * t4165 + t4082 * t4156;
t4106 = cos(t4144);
t4164 = cos(qJ(2,2));
t4076 = 0.1e1 / (t4164 * pkin(2) + pkin(3) * t4106 + pkin(1));
t4277 = t4076 * t4156;
t4208 = t4118 * t4277;
t4194 = t4201 * t4208;
t4209 = t4115 * t4277;
t4195 = t4201 * t4209;
t4276 = t4076 * t4165;
t4214 = t4201 * t4276;
t4148 = legFrame(3,2);
t4114 = sin(t4148);
t4117 = cos(t4148);
t4081 = t4117 * g(1) - t4114 * g(2);
t4162 = cos(qJ(1,3));
t4202 = g(3) * t4162 + t4081 * t4153;
t4105 = cos(t4143);
t4161 = cos(qJ(2,3));
t4075 = 0.1e1 / (t4161 * pkin(2) + pkin(3) * t4105 + pkin(1));
t4279 = t4075 * t4153;
t4210 = t4117 * t4279;
t4198 = t4202 * t4210;
t4211 = t4114 * t4279;
t4199 = t4202 * t4211;
t4278 = t4075 * t4162;
t4216 = t4202 * t4278;
t4324 = -g(3) * t4159 + t4083 * t4168;
t4323 = -g(3) * t4156 + t4082 * t4165;
t4322 = -g(3) * t4153 + t4081 * t4162;
t4321 = 0.2e1 * pkin(3);
t4319 = 0.2e1 * t4161 ^ 2;
t4318 = 0.2e1 * t4164 ^ 2;
t4317 = 0.2e1 * t4167 ^ 2;
t4160 = cos(qJ(3,3));
t4137 = t4160 ^ 2;
t4126 = pkin(3) * t4137;
t4163 = cos(qJ(3,2));
t4139 = t4163 ^ 2;
t4127 = pkin(3) * t4139;
t4166 = cos(qJ(3,1));
t4141 = t4166 ^ 2;
t4128 = pkin(3) * t4141;
t4312 = t4151 * pkin(3);
t4310 = t4154 * pkin(3);
t4308 = t4157 * pkin(3);
t4307 = t4160 * pkin(2);
t4120 = t4160 * pkin(3);
t4306 = t4163 * pkin(2);
t4121 = t4163 * pkin(3);
t4305 = t4166 * pkin(2);
t4122 = t4166 * pkin(3);
t4242 = pkin(1) * t4162 + t4153 * t4133;
t4249 = t4151 * t4152;
t4045 = t4242 * t4249 + (t4137 - 0.1e1) * t4162 * pkin(3);
t4288 = t4045 * t4114;
t4287 = t4045 * t4117;
t4241 = pkin(1) * t4165 + t4156 * t4133;
t4248 = t4154 * t4155;
t4046 = t4241 * t4248 + (t4139 - 0.1e1) * t4165 * pkin(3);
t4286 = t4046 * t4115;
t4285 = t4046 * t4118;
t4240 = pkin(1) * t4168 + t4159 * t4133;
t4247 = t4157 * t4158;
t4047 = t4240 * t4247 + (t4141 - 0.1e1) * t4168 * pkin(3);
t4284 = t4047 * t4116;
t4283 = t4047 * t4119;
t4232 = pkin(3) * t4249;
t4099 = t4120 + pkin(2);
t4265 = t4099 * t4161;
t4282 = 0.1e1 / (pkin(1) - t4232 + t4265) / t4151;
t4231 = pkin(3) * t4248;
t4100 = t4121 + pkin(2);
t4262 = t4100 * t4164;
t4281 = 0.1e1 / (pkin(1) - t4231 + t4262) / t4154;
t4230 = pkin(3) * t4247;
t4101 = t4122 + pkin(2);
t4259 = t4101 * t4167;
t4280 = 0.1e1 / (pkin(1) - t4230 + t4259) / t4157;
t4171 = pkin(2) / 0.2e1;
t4270 = (t4120 + t4171) * t4151;
t4269 = (t4121 + t4171) * t4154;
t4268 = (t4122 + t4171) * t4157;
t4267 = t4099 * t4114;
t4266 = t4099 * t4117;
t4264 = t4100 * t4115;
t4263 = t4100 * t4118;
t4261 = t4101 * t4116;
t4260 = t4101 * t4119;
t4258 = t4114 * t4162;
t4257 = t4115 * t4165;
t4256 = t4116 * t4168;
t4255 = t4117 * t4162;
t4254 = t4118 * t4165;
t4253 = t4119 * t4168;
t4181 = pkin(3) ^ 2;
t4252 = t4137 * t4181;
t4251 = t4139 * t4181;
t4250 = t4141 * t4181;
t4090 = pkin(1) * t4152 - t4312;
t4246 = t4160 * t4090;
t4091 = pkin(1) * t4155 - t4310;
t4245 = t4163 * t4091;
t4092 = pkin(1) * t4158 - t4308;
t4244 = t4166 * t4092;
t4182 = 0.1e1 / pkin(3);
t4184 = 0.1e1 / pkin(2);
t4243 = t4182 * t4184;
t4239 = -t4181 / 0.2e1 + t4183 / 0.2e1;
t4238 = pkin(2) * t4120;
t4237 = pkin(2) * t4121;
t4236 = pkin(2) * t4122;
t4235 = t4099 * t4312;
t4234 = t4100 * t4310;
t4233 = t4101 * t4308;
t4078 = t4114 * g(1) + t4117 * g(2);
t4030 = -t4078 * t4105 + t4102 * t4322;
t4229 = t4030 * t4282;
t4031 = t4078 * t4102 + t4105 * t4322;
t4228 = t4031 * t4282;
t4079 = t4115 * g(1) + t4118 * g(2);
t4032 = -t4079 * t4106 + t4103 * t4323;
t4227 = t4032 * t4281;
t4033 = t4079 * t4103 + t4106 * t4323;
t4226 = t4033 * t4281;
t4080 = t4116 * g(1) + t4119 * g(2);
t4034 = -t4080 * t4107 + t4104 * t4324;
t4225 = t4034 * t4280;
t4035 = t4080 * t4104 + t4107 * t4324;
t4224 = t4035 * t4280;
t4036 = -t4078 * t4161 + t4152 * t4322;
t4223 = t4036 * t4282;
t4037 = t4078 * t4152 + t4161 * t4322;
t4222 = t4037 * t4282;
t4038 = -t4079 * t4164 + t4155 * t4323;
t4221 = t4038 * t4281;
t4039 = t4079 * t4155 + t4164 * t4323;
t4220 = t4039 * t4281;
t4040 = -t4080 * t4167 + t4158 * t4324;
t4219 = t4040 * t4280;
t4041 = t4080 * t4158 + t4167 * t4324;
t4218 = t4041 * t4280;
t4217 = t4202 * t4279;
t4215 = t4201 * t4277;
t4213 = t4200 * t4275;
t4205 = t4162 * t4249;
t4204 = t4165 * t4248;
t4203 = t4168 * t4247;
t4197 = t4152 * t4217;
t4196 = t4161 * t4217;
t4193 = t4155 * t4215;
t4192 = t4164 * t4215;
t4189 = t4158 * t4213;
t4188 = t4167 * t4213;
t4060 = t4205 * t4321 - t4242;
t4187 = pkin(2) * t4205 + t4060 * t4160;
t4061 = t4204 * t4321 - t4241;
t4186 = pkin(2) * t4204 + t4061 * t4163;
t4062 = t4203 * t4321 - t4240;
t4185 = pkin(2) * t4203 + t4062 * t4166;
t4170 = -pkin(3) / 0.2e1;
t4132 = -t4181 + t4183;
t4086 = t4128 + t4305 / 0.2e1 + t4170;
t4085 = t4127 + t4306 / 0.2e1 + t4170;
t4084 = t4126 + t4307 / 0.2e1 + t4170;
t4074 = t4236 + t4239 + t4250;
t4073 = t4237 + t4239 + t4251;
t4072 = t4238 + t4239 + t4252;
t4053 = t4309 + (-pkin(3) + t4305 + 0.2e1 * t4128) * t4158;
t4052 = t4311 + (-pkin(3) + t4306 + 0.2e1 * t4127) * t4155;
t4051 = t4313 + (-pkin(3) + t4307 + 0.2e1 * t4126) * t4152;
t4050 = pkin(1) * t4308 + (t4132 + 0.2e1 * t4236 + 0.2e1 * t4250) * t4158;
t4049 = pkin(1) * t4310 + (t4132 + 0.2e1 * t4237 + 0.2e1 * t4251) * t4155;
t4048 = pkin(1) * t4312 + (t4132 + 0.2e1 * t4238 + 0.2e1 * t4252) * t4152;
t4029 = t4074 * t4159 * t4317 + ((pkin(1) - 0.2e1 * t4230) * t4159 - t4168 * t4133) * t4259 - pkin(3) * ((pkin(1) * t4247 - pkin(3) + t4128) * t4159 - t4133 * t4203);
t4028 = t4073 * t4156 * t4318 + ((pkin(1) - 0.2e1 * t4231) * t4156 - t4165 * t4133) * t4262 - pkin(3) * ((pkin(1) * t4248 - pkin(3) + t4127) * t4156 - t4133 * t4204);
t4027 = t4072 * t4153 * t4319 + ((pkin(1) - 0.2e1 * t4232) * t4153 - t4162 * t4133) * t4265 - pkin(3) * ((pkin(1) * t4249 - pkin(3) + t4126) * t4153 - t4133 * t4205);
t4023 = (t4086 * t4253 + t4116 * t4268) * t4317 + (t4116 * t4053 - t4185 * t4119) * t4167 - t4283 + t4116 * t4244;
t4022 = (-t4086 * t4256 + t4119 * t4268) * t4317 + (t4119 * t4053 + t4185 * t4116) * t4167 + t4284 + t4119 * t4244;
t4021 = (t4085 * t4254 + t4115 * t4269) * t4318 + (t4115 * t4052 - t4186 * t4118) * t4164 - t4285 + t4115 * t4245;
t4020 = (-t4085 * t4257 + t4118 * t4269) * t4318 + (t4118 * t4052 + t4186 * t4115) * t4164 + t4286 + t4118 * t4245;
t4019 = (t4084 * t4255 + t4114 * t4270) * t4319 + (t4114 * t4051 - t4187 * t4117) * t4161 - t4287 + t4114 * t4246;
t4018 = (-t4084 * t4258 + t4117 * t4270) * t4319 + (t4117 * t4051 + t4187 * t4114) * t4161 + t4288 + t4117 * t4246;
t4017 = (-t4074 * t4253 - t4116 * t4233) * t4317 + (-t4116 * t4050 + t4062 * t4260) * t4167 + pkin(3) * t4283 - t4092 * t4261;
t4016 = (t4074 * t4256 - t4119 * t4233) * t4317 + (-t4119 * t4050 - t4062 * t4261) * t4167 - pkin(3) * t4284 - t4092 * t4260;
t4015 = (-t4073 * t4254 - t4115 * t4234) * t4318 + (-t4115 * t4049 + t4061 * t4263) * t4164 + pkin(3) * t4285 - t4091 * t4264;
t4014 = (t4073 * t4257 - t4118 * t4234) * t4318 + (-t4118 * t4049 - t4061 * t4264) * t4164 - pkin(3) * t4286 - t4091 * t4263;
t4013 = (-t4072 * t4255 - t4114 * t4235) * t4319 + (-t4114 * t4048 + t4060 * t4266) * t4161 + pkin(3) * t4287 - t4090 * t4267;
t4012 = (t4072 * t4258 - t4117 * t4235) * t4319 + (-t4117 * t4048 - t4060 * t4267) * t4161 - pkin(3) * t4288 - t4090 * t4266;
t1 = [0, -t4194 - t4198 - t4190, -t4206 * t4324 - t4208 * t4323 - t4210 * t4322, 0, 0, 0, 0, 0, -t4117 * t4196 - t4118 * t4192 - t4119 * t4188 + (t4019 * t4223 + t4021 * t4221 + t4023 * t4219) * t4184, t4117 * t4197 + t4118 * t4193 + t4119 * t4189 + (t4019 * t4222 + t4021 * t4220 + t4023 * t4218) * t4184, 0, 0, 0, 0, 0, -t4105 * t4198 - t4106 * t4194 - t4107 * t4190 + (t4019 * t4229 + t4021 * t4227 + t4023 * t4225 + (t4013 * t4229 + t4015 * t4227 + t4017 * t4225) * t4182) * t4184, t4102 * t4198 + t4103 * t4194 + t4104 * t4190 + (t4019 * t4228 + t4021 * t4226 + t4023 * t4224 + (t4013 * t4228 + t4015 * t4226 + t4017 * t4224) * t4182) * t4184, -g(1); 0, t4195 + t4199 + t4191, t4207 * t4324 + t4209 * t4323 + t4211 * t4322, 0, 0, 0, 0, 0, t4114 * t4196 + t4115 * t4192 + t4116 * t4188 + (t4018 * t4223 + t4020 * t4221 + t4022 * t4219) * t4184, -t4114 * t4197 - t4115 * t4193 - t4116 * t4189 + (t4018 * t4222 + t4020 * t4220 + t4022 * t4218) * t4184, 0, 0, 0, 0, 0, t4105 * t4199 + t4106 * t4195 + t4107 * t4191 + (t4018 * t4229 + t4020 * t4227 + t4022 * t4225 + (t4012 * t4229 + t4014 * t4227 + t4016 * t4225) * t4182) * t4184, -t4102 * t4199 - t4103 * t4195 - t4104 * t4191 + (t4018 * t4228 + t4020 * t4226 + t4022 * t4224 + (t4012 * t4228 + t4014 * t4226 + t4016 * t4224) * t4182) * t4184, -g(2); 0, -t4214 - t4216 - t4212, -t4274 * t4324 - t4276 * t4323 - t4278 * t4322, 0, 0, 0, 0, 0, t4036 * t4325 + t4038 * t4326 + t4040 * t4327 - t4161 * t4216 - t4164 * t4214 - t4167 * t4212, t4037 * t4325 + t4039 * t4326 + t4041 * t4327 + t4152 * t4216 + t4155 * t4214 + t4158 * t4212, 0, 0, 0, 0, 0, -t4105 * t4216 - t4106 * t4214 - t4107 * t4212 + (t4027 * t4229 + t4028 * t4227 + t4029 * t4225) * t4243 + t4030 * t4325 + t4032 * t4326 + t4034 * t4327, t4102 * t4216 + t4103 * t4214 + t4104 * t4212 + (t4027 * t4228 + t4028 * t4226 + t4029 * t4224) * t4243 + t4031 * t4325 + t4033 * t4326 + t4035 * t4327, -g(3);];
tau_reg  = t1;
