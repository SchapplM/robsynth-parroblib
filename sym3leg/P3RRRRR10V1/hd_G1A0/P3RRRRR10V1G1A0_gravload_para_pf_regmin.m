% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRRRR10V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d1,d2,d4]';
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
% Datum: 2020-08-06 22:09
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRRRR10V1G1A0_gravload_para_pf_regmin(xP, qJ, g, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_regmin: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V1G1A0_gravload_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 21:55:52
% EndTime: 2020-08-06 21:56:00
% DurationCPUTime: 7.37s
% Computational Cost: add. (2781->429), mult. (7185->764), div. (171->13), fcn. (6960->26), ass. (0->351)
t4252 = cos(pkin(3));
t4241 = t4252 ^ 2;
t4458 = t4241 * pkin(6);
t4251 = sin(pkin(3));
t4415 = t4251 * t4252;
t4271 = cos(qJ(3,1));
t4272 = cos(qJ(2,1));
t4263 = sin(qJ(2,1));
t4262 = sin(qJ(3,1));
t4371 = t4272 * t4262;
t4317 = t4263 * t4371;
t4374 = t4271 * t4251;
t4367 = pkin(2) * t4374;
t4248 = t4271 ^ 2;
t4372 = t4272 * t4248;
t4391 = t4262 * t4251;
t4231 = t4263 * pkin(6);
t4479 = t4231 + pkin(1);
t4489 = (pkin(2) * t4372 + t4271 * t4479) * t4252 - t4317 * t4367 + pkin(6) * (t4272 - 0.1e1) * (t4272 + 0.1e1) * t4391;
t4268 = cos(qJ(3,2));
t4269 = cos(qJ(2,2));
t4259 = sin(qJ(3,2));
t4260 = sin(qJ(2,2));
t4395 = t4259 * t4260;
t4319 = t4251 * t4395;
t4379 = t4268 * t4269;
t4364 = pkin(2) * t4379;
t4245 = t4268 ^ 2;
t4378 = t4269 * t4245;
t4396 = t4259 * t4251;
t4230 = t4260 * pkin(6);
t4480 = t4230 + pkin(1);
t4488 = (pkin(2) * t4378 + t4268 * t4480) * t4252 - t4319 * t4364 + pkin(6) * (t4269 - 0.1e1) * (t4269 + 0.1e1) * t4396;
t4265 = cos(qJ(3,3));
t4266 = cos(qJ(2,3));
t4256 = sin(qJ(3,3));
t4257 = sin(qJ(2,3));
t4400 = t4256 * t4257;
t4321 = t4251 * t4400;
t4385 = t4265 * t4266;
t4366 = pkin(2) * t4385;
t4242 = t4265 ^ 2;
t4384 = t4266 * t4242;
t4401 = t4256 * t4251;
t4229 = t4257 * pkin(6);
t4481 = t4229 + pkin(1);
t4487 = (pkin(2) * t4384 + t4265 * t4481) * t4252 - t4321 * t4366 + pkin(6) * (t4266 - 0.1e1) * (t4266 + 0.1e1) * t4401;
t4486 = -0.2e1 * pkin(6);
t4485 = 0.2e1 * pkin(6);
t4455 = t4262 * pkin(2);
t4219 = pkin(1) * t4455;
t4373 = t4271 * t4272;
t4389 = t4263 * t4271;
t4184 = pkin(2) * t4389 - t4272 * pkin(6);
t4422 = t4184 * t4252;
t4127 = 0.1e1 / (pkin(1) * t4422 + (-t4219 + (pkin(2) * t4373 + t4231) * pkin(5)) * t4251);
t4456 = t4259 * pkin(2);
t4218 = pkin(1) * t4456;
t4394 = t4260 * t4268;
t4183 = pkin(2) * t4394 - t4269 * pkin(6);
t4423 = t4183 * t4252;
t4126 = 0.1e1 / (pkin(1) * t4423 + (-t4218 + (t4364 + t4230) * pkin(5)) * t4251);
t4457 = t4256 * pkin(2);
t4217 = pkin(1) * t4457;
t4399 = t4257 * t4265;
t4182 = pkin(2) * t4399 - t4266 * pkin(6);
t4424 = t4182 * t4252;
t4125 = 0.1e1 / (pkin(1) * t4424 + (-t4217 + (t4366 + t4229) * pkin(5)) * t4251);
t4484 = pkin(1) * t4257;
t4483 = pkin(1) * t4260;
t4482 = pkin(1) * t4263;
t4255 = legFrame(1,3);
t4223 = sin(t4255);
t4226 = cos(t4255);
t4264 = sin(qJ(1,1));
t4273 = cos(qJ(1,1));
t4157 = t4223 * t4264 - t4273 * t4226;
t4274 = pkin(6) ^ 2;
t4275 = pkin(2) ^ 2;
t4307 = t4248 * t4275 - t4274;
t4478 = t4307 * t4157;
t4254 = legFrame(2,3);
t4222 = sin(t4254);
t4225 = cos(t4254);
t4261 = sin(qJ(1,2));
t4270 = cos(qJ(1,2));
t4156 = t4222 * t4261 - t4270 * t4225;
t4308 = t4245 * t4275 - t4274;
t4477 = t4308 * t4156;
t4253 = legFrame(3,3);
t4221 = sin(t4253);
t4224 = cos(t4253);
t4258 = sin(qJ(1,3));
t4267 = cos(qJ(1,3));
t4155 = t4221 * t4258 - t4267 * t4224;
t4309 = t4242 * t4275 - t4274;
t4476 = t4309 * t4155;
t4244 = t4266 ^ 2;
t4188 = (t4244 - 0.2e1) * t4457 - pkin(5);
t4204 = t4256 * pkin(5) + pkin(2);
t4464 = pkin(2) * t4242;
t4312 = t4204 - 0.2e1 * t4464;
t4383 = t4266 * t4256;
t4361 = pkin(6) * t4383;
t4469 = t4385 * t4458 - t4188 * t4265 * t4415 + (t4312 * t4241 - t4361 * t4415 - t4204 + t4464) * t4257;
t4247 = t4269 ^ 2;
t4189 = (t4247 - 0.2e1) * t4456 - pkin(5);
t4207 = t4259 * pkin(5) + pkin(2);
t4463 = pkin(2) * t4245;
t4311 = t4207 - 0.2e1 * t4463;
t4377 = t4269 * t4259;
t4360 = pkin(6) * t4377;
t4468 = t4379 * t4458 - t4189 * t4268 * t4415 + (t4311 * t4241 - t4360 * t4415 - t4207 + t4463) * t4260;
t4250 = t4272 ^ 2;
t4190 = (t4250 - 0.2e1) * t4455 - pkin(5);
t4210 = t4262 * pkin(5) + pkin(2);
t4462 = pkin(2) * t4248;
t4310 = t4210 - 0.2e1 * t4462;
t4467 = (-pkin(6) * t4317 - t4190 * t4271) * t4415 + t4373 * t4458 + (t4310 * t4241 - t4210 + t4462) * t4263;
t4466 = (t4241 - 0.1e1) * t4485;
t4465 = pkin(1) * pkin(2);
t4461 = pkin(2) * t4265;
t4460 = pkin(2) * t4268;
t4459 = pkin(2) * t4271;
t4220 = t4251 * g(3);
t4359 = pkin(6) * t4415;
t4095 = (t4266 * t4359 + t4188 * t4241 + pkin(5) + (-t4244 + 0.1e1) * t4457) * t4265 - t4481 * t4383 + (t4241 * t4361 + t4312 * t4415) * t4257;
t4216 = pkin(1) * t4252 * pkin(6);
t4409 = t4252 * t4257;
t4414 = t4251 * t4266;
t4418 = t4251 * (-pkin(5) * t4229 + t4217);
t4123 = 0.1e1 / ((pkin(1) * t4409 + pkin(5) * t4414) * t4461 - t4266 * t4216 - t4418);
t4451 = t4095 * t4123;
t4096 = (t4269 * t4359 + t4189 * t4241 + pkin(5) + (-t4247 + 0.1e1) * t4456) * t4268 - t4480 * t4377 + (t4241 * t4360 + t4311 * t4415) * t4260;
t4407 = t4252 * t4260;
t4413 = t4251 * t4269;
t4417 = t4251 * (-pkin(5) * t4230 + t4218);
t4124 = 0.1e1 / ((pkin(1) * t4407 + pkin(5) * t4413) * t4460 - t4216 * t4269 - t4417);
t4450 = t4096 * t4124;
t4097 = (t4272 * t4359 + t4190 * t4241 + pkin(5) + (-t4250 + 0.1e1) * t4455) * t4271 - t4479 * t4371 + (t4310 * t4415 + t4371 * t4458) * t4263;
t4405 = t4252 * t4263;
t4412 = t4251 * t4272;
t4416 = t4251 * (-pkin(5) * t4231 + t4219);
t4122 = 0.1e1 / ((pkin(1) * t4405 + pkin(5) * t4412) * t4459 - t4216 * t4272 - t4416);
t4449 = t4097 * t4122;
t4381 = t4267 * t4266;
t4398 = t4258 * t4257;
t4167 = t4252 * t4381 - t4398;
t4382 = t4267 * t4257;
t4397 = t4258 * t4266;
t4173 = t4252 * t4397 + t4382;
t4176 = t4221 * g(1) - t4224 * g(2);
t4179 = t4224 * g(1) + t4221 * g(2);
t4107 = -g(3) * t4414 + t4176 * t4167 + t4179 * t4173;
t4243 = 0.1e1 / t4265;
t4448 = t4107 * t4243;
t4164 = t4252 * t4398 - t4381;
t4168 = t4252 * t4382 + t4397;
t4108 = -t4179 * t4164 - t4176 * t4168 + t4257 * t4220;
t4447 = t4108 * t4243;
t4375 = t4270 * t4269;
t4393 = t4261 * t4260;
t4169 = t4252 * t4375 - t4393;
t4376 = t4270 * t4260;
t4392 = t4261 * t4269;
t4174 = t4252 * t4392 + t4376;
t4177 = t4222 * g(1) - t4225 * g(2);
t4180 = t4225 * g(1) + t4222 * g(2);
t4109 = -g(3) * t4413 + t4177 * t4169 + t4180 * t4174;
t4246 = 0.1e1 / t4268;
t4446 = t4109 * t4246;
t4165 = t4252 * t4393 - t4375;
t4170 = t4252 * t4376 + t4392;
t4110 = -t4180 * t4165 - t4177 * t4170 + t4260 * t4220;
t4445 = t4110 * t4246;
t4369 = t4273 * t4272;
t4388 = t4264 * t4263;
t4171 = t4252 * t4369 - t4388;
t4370 = t4273 * t4263;
t4387 = t4264 * t4272;
t4175 = t4252 * t4387 + t4370;
t4178 = t4223 * g(1) - t4226 * g(2);
t4181 = t4226 * g(1) + t4223 * g(2);
t4111 = -g(3) * t4412 + t4178 * t4171 + t4181 * t4175;
t4249 = 0.1e1 / t4271;
t4444 = t4111 * t4249;
t4166 = t4252 * t4388 - t4369;
t4172 = t4252 * t4370 + t4387;
t4112 = -t4181 * t4166 - t4178 * t4172 + t4263 * t4220;
t4443 = t4112 * t4249;
t4315 = t4252 * t4399;
t4368 = pkin(2) * pkin(5) * t4251;
t4442 = 0.1e1 / ((-t4265 * t4368 + t4216) * t4266 - t4315 * t4465 + t4418) * t4243;
t4314 = t4252 * t4394;
t4441 = 0.1e1 / ((-t4268 * t4368 + t4216) * t4269 - t4314 * t4465 + t4417) * t4246;
t4313 = t4252 * t4389;
t4440 = 0.1e1 / ((-pkin(5) * t4367 + t4216) * t4272 - t4313 * t4465 + t4416) * t4249;
t4439 = t4122 * t4249;
t4438 = t4123 * t4243;
t4437 = t4124 * t4246;
t4131 = t4176 * t4267 + t4179 * t4258;
t4436 = t4125 * t4131;
t4134 = -t4176 * t4258 + t4179 * t4267;
t4435 = t4125 * t4134;
t4132 = t4177 * t4270 + t4180 * t4261;
t4434 = t4126 * t4132;
t4135 = -t4177 * t4261 + t4180 * t4270;
t4433 = t4126 * t4135;
t4133 = t4178 * t4273 + t4181 * t4264;
t4432 = t4127 * t4133;
t4136 = -t4178 * t4264 + t4181 * t4273;
t4431 = t4127 * t4136;
t4410 = t4252 * t4256;
t4430 = t4125 * (pkin(2) * t4410 + t4182 * t4251);
t4408 = t4252 * t4259;
t4429 = t4126 * (pkin(2) * t4408 + t4183 * t4251);
t4406 = t4252 * t4262;
t4428 = t4127 * (pkin(2) * t4406 + t4184 * t4251);
t4205 = pkin(5) + t4457;
t4421 = t4205 * t4251;
t4208 = pkin(5) + t4456;
t4420 = t4208 * t4251;
t4211 = pkin(5) + t4455;
t4419 = t4211 * t4251;
t4276 = 0.1e1 / pkin(2);
t4411 = t4251 * t4276;
t4404 = t4252 * t4265;
t4403 = t4252 * t4268;
t4402 = t4252 * t4271;
t4390 = t4262 * t4263;
t4386 = t4265 * t4251;
t4380 = t4268 * t4251;
t4358 = pkin(2) * t4404;
t4357 = pkin(2) * t4403;
t4356 = pkin(2) * t4402;
t4213 = pkin(6) + t4484;
t4214 = pkin(6) + t4483;
t4215 = pkin(6) + t4482;
t4288 = t4164 * t4256 + t4258 * t4386;
t4080 = -t4179 * t4288 - ((t4252 * t4400 + t4386) * t4267 + t4258 * t4383) * t4176 - g(3) * (-t4321 + t4404);
t4349 = t4080 * t4442;
t4161 = t4315 - t4401;
t4137 = -t4161 * t4258 + t4265 * t4381;
t4081 = t4137 * t4179 - (t4161 * t4267 + t4258 * t4385) * t4176 + g(3) * (t4257 * t4386 + t4410);
t4348 = t4081 * t4442;
t4162 = t4314 - t4396;
t4138 = -t4162 * t4261 + t4268 * t4375;
t4082 = t4138 * t4180 - t4177 * (t4162 * t4270 + t4261 * t4379) + g(3) * (t4260 * t4380 + t4408);
t4347 = t4082 * t4441;
t4163 = t4313 - t4391;
t4139 = -t4163 * t4264 + t4271 * t4369;
t4083 = t4139 * t4181 - t4178 * (t4163 * t4273 + t4264 * t4373) + g(3) * (t4263 * t4374 + t4406);
t4346 = t4083 * t4440;
t4287 = t4165 * t4259 + t4261 * t4380;
t4084 = -t4287 * t4180 - t4177 * ((t4252 * t4395 + t4380) * t4270 + t4261 * t4377) - g(3) * (-t4319 + t4403);
t4345 = t4084 * t4441;
t4286 = t4166 * t4262 + t4264 * t4374;
t4085 = -t4286 * t4181 - t4178 * ((t4252 * t4390 + t4374) * t4273 + t4264 * t4371) - g(3) * (-t4251 * t4390 + t4402);
t4344 = t4085 * t4440;
t4324 = t4205 * t4415;
t4328 = (t4252 + 0.1e1) * (t4252 - 0.1e1) * t4275;
t4343 = (-t4257 * t4328 * t4384 + (t4244 * t4466 + t4266 * t4324 + t4213 - t4458) * t4461 - ((-t4241 * t4229 + t4481) * t4266 - t4257 * t4324) * pkin(6)) * t4438;
t4323 = t4208 * t4415;
t4342 = (-t4260 * t4328 * t4378 + (t4247 * t4466 + t4269 * t4323 + t4214 - t4458) * t4460 - ((-t4241 * t4230 + t4480) * t4269 - t4260 * t4323) * pkin(6)) * t4437;
t4322 = t4211 * t4415;
t4341 = (-t4263 * t4328 * t4372 + (t4250 * t4466 + t4272 * t4322 + t4215 - t4458) * t4459 - ((-t4241 * t4231 + t4479) * t4272 - t4263 * t4322) * pkin(6)) * t4439;
t4340 = t4095 * t4438;
t4339 = t4096 * t4437;
t4338 = t4097 * t4439;
t4191 = t4258 * g(1) - t4267 * g(2);
t4192 = t4267 * g(1) + t4258 * g(2);
t4104 = (t4220 + (-t4191 * t4224 - t4192 * t4221) * t4252) * t4266 + t4257 * (t4191 * t4221 - t4192 * t4224);
t4337 = t4104 * t4243 * t4256;
t4193 = t4261 * g(1) - t4270 * g(2);
t4194 = t4270 * g(1) + t4261 * g(2);
t4105 = (t4220 + (-t4193 * t4225 - t4194 * t4222) * t4252) * t4269 + t4260 * (t4193 * t4222 - t4194 * t4225);
t4336 = t4105 * t4246 * t4259;
t4195 = g(1) * t4264 - g(2) * t4273;
t4196 = g(1) * t4273 + g(2) * t4264;
t4106 = (t4220 + (-t4195 * t4226 - t4196 * t4223) * t4252) * t4272 + t4263 * (t4195 * t4223 - t4196 * t4226);
t4335 = t4106 * t4249 * t4262;
t4334 = t4155 * t4401;
t4333 = t4156 * t4396;
t4332 = t4157 * t4391;
t4158 = t4221 * t4267 + t4258 * t4224;
t4331 = t4158 * t4401;
t4159 = t4222 * t4270 + t4261 * t4225;
t4330 = t4159 * t4396;
t4160 = t4223 * t4273 + t4264 * t4226;
t4329 = t4160 * t4391;
t4327 = t4267 * t4421;
t4326 = t4270 * t4420;
t4325 = t4273 * t4419;
t4320 = t4258 * t4421;
t4318 = t4261 * t4420;
t4316 = t4264 * t4419;
t4306 = t4155 * t4358;
t4305 = t4156 * t4357;
t4304 = t4157 * t4356;
t4303 = t4158 * t4358;
t4302 = t4159 * t4357;
t4301 = t4160 * t4356;
t4297 = t4309 * t4158;
t4296 = t4308 * t4159;
t4295 = t4307 * t4160;
t4206 = 0.2e1 * t4229 + pkin(1);
t4285 = t4206 * t4267 + t4320;
t4209 = 0.2e1 * t4230 + pkin(1);
t4284 = t4209 * t4270 + t4318;
t4212 = 0.2e1 * t4231 + pkin(1);
t4283 = t4212 * t4273 + t4316;
t4282 = t4213 * t4267 + t4257 * t4320;
t4281 = t4214 * t4270 + t4260 * t4318;
t4280 = t4215 * t4273 + t4263 * t4316;
t4151 = t4264 * t4212 - t4325;
t4150 = t4261 * t4209 - t4326;
t4149 = t4258 * t4206 - t4327;
t4142 = t4264 * t4215 - t4263 * t4325;
t4141 = t4261 * t4214 - t4260 * t4326;
t4140 = t4258 * t4213 - t4257 * t4327;
t4118 = -t4178 * t4166 + t4181 * t4172;
t4117 = t4181 * t4171 - t4178 * t4175;
t4116 = -t4177 * t4165 + t4180 * t4170;
t4115 = t4180 * t4169 - t4177 * t4174;
t4114 = -t4176 * t4164 + t4179 * t4168;
t4113 = t4179 * t4167 - t4176 * t4173;
t4103 = t4181 * (t4172 * t4271 - t4273 * t4391) + t4178 * t4139;
t4102 = -t4181 * (t4172 * t4262 + t4273 * t4374) + t4178 * t4286;
t4101 = t4180 * (t4170 * t4268 - t4270 * t4396) + t4138 * t4177;
t4100 = -t4180 * (t4170 * t4259 + t4270 * t4380) + t4177 * t4287;
t4099 = t4179 * (t4168 * t4265 - t4267 * t4401) + t4137 * t4176;
t4098 = -t4179 * (t4168 * t4256 + t4267 * t4386) + t4176 * t4288;
t4091 = t4157 * t4231 + t4160 * t4422 + (t4157 * t4373 - t4329) * pkin(2);
t4090 = t4156 * t4230 + t4159 * t4423 + (t4156 * t4379 - t4330) * pkin(2);
t4089 = t4155 * t4229 + t4158 * t4424 + (t4155 * t4385 - t4331) * pkin(2);
t4088 = -t4160 * t4231 + t4157 * t4422 + (-t4160 * t4373 - t4332) * pkin(2);
t4087 = -t4159 * t4230 + t4156 * t4423 + (-t4159 * t4379 - t4333) * pkin(2);
t4086 = -t4158 * t4229 + t4155 * t4424 + (-t4158 * t4385 - t4334) * pkin(2);
t4079 = (t4304 * t4485 + t4295) * t4250 + ((t4151 * t4226 + t4223 * t4283) * t4459 - t4405 * t4478) * t4272 + (t4142 * t4226 + t4223 * t4280 - t4304) * pkin(6);
t4078 = (t4305 * t4485 + t4296) * t4247 + ((t4150 * t4225 + t4222 * t4284) * t4460 - t4407 * t4477) * t4269 + (t4141 * t4225 + t4222 * t4281 - t4305) * pkin(6);
t4077 = (t4306 * t4485 + t4297) * t4244 + ((t4149 * t4224 + t4221 * t4285) * t4461 - t4409 * t4476) * t4266 + (t4140 * t4224 + t4221 * t4282 - t4306) * pkin(6);
t4076 = (t4301 * t4486 + t4478) * t4250 + ((t4223 * t4151 - t4283 * t4226) * t4459 + t4295 * t4405) * t4272 + (t4223 * t4142 - t4280 * t4226 + t4301) * pkin(6);
t4075 = (t4302 * t4486 + t4477) * t4247 + ((t4222 * t4150 - t4284 * t4225) * t4460 + t4296 * t4407) * t4269 + (t4222 * t4141 - t4281 * t4225 + t4302) * pkin(6);
t4074 = (t4303 * t4486 + t4476) * t4244 + ((t4221 * t4149 - t4285 * t4224) * t4461 + t4297 * t4409) * t4266 + (t4221 * t4140 - t4282 * t4224 + t4303) * pkin(6);
t4073 = -t4489 * t4157 + t4467 * t4160 + t4332 * t4482;
t4072 = -t4488 * t4156 + t4468 * t4159 + t4333 * t4483;
t4071 = -t4487 * t4155 + t4469 * t4158 + t4334 * t4484;
t4070 = t4467 * t4157 + t4489 * t4160 - t4329 * t4482;
t4069 = t4468 * t4156 + t4488 * t4159 - t4330 * t4483;
t4068 = t4469 * t4155 + t4487 * t4158 - t4331 * t4484;
t1 = [0, -t4089 * t4436 - t4090 * t4434 - t4091 * t4432, -t4089 * t4435 - t4090 * t4433 - t4091 * t4431, 0, 0, 0, 0, 0, -(t4073 * t4444 + t4091 * t4118) * t4127 - (t4072 * t4446 + t4090 * t4116) * t4126 - (t4071 * t4448 + t4089 * t4114) * t4125, -(t4073 * t4443 + t4091 * t4117) * t4127 - (t4072 * t4445 + t4090 * t4115) * t4126 - (t4071 * t4447 + t4089 * t4113) * t4125, 0, 0, 0, 0, 0, -(-t4073 * t4106 + t4091 * t4103) * t4127 - (-t4072 * t4105 + t4090 * t4101) * t4126 - (-t4071 * t4104 + t4089 * t4099) * t4125 + (t4074 * t4349 + t4075 * t4345 + t4076 * t4344) * t4411, -(t4073 * t4335 + t4091 * t4102) * t4127 - (t4072 * t4336 + t4090 * t4100) * t4126 - (t4071 * t4337 + t4089 * t4098) * t4125 + (t4074 * t4348 + t4075 * t4347 + t4076 * t4346) * t4411, -g(1); 0, -t4086 * t4436 - t4087 * t4434 - t4088 * t4432, -t4086 * t4435 - t4087 * t4433 - t4088 * t4431, 0, 0, 0, 0, 0, -(t4070 * t4444 + t4088 * t4118) * t4127 - (t4069 * t4446 + t4087 * t4116) * t4126 - (t4068 * t4448 + t4086 * t4114) * t4125, -(t4070 * t4443 + t4088 * t4117) * t4127 - (t4069 * t4445 + t4087 * t4115) * t4126 - (t4068 * t4447 + t4086 * t4113) * t4125, 0, 0, 0, 0, 0, -(-t4070 * t4106 + t4088 * t4103) * t4127 - (-t4069 * t4105 + t4087 * t4101) * t4126 - (-t4068 * t4104 + t4086 * t4099) * t4125 + (-t4077 * t4349 - t4078 * t4345 - t4079 * t4344) * t4411, -(t4070 * t4335 + t4088 * t4102) * t4127 - (t4069 * t4336 + t4087 * t4100) * t4126 - (t4068 * t4337 + t4086 * t4098) * t4125 + (-t4077 * t4348 - t4078 * t4347 - t4079 * t4346) * t4411, -g(2); 0, t4131 * t4430 + t4132 * t4429 + t4133 * t4428, t4134 * t4430 + t4135 * t4429 + t4136 * t4428, 0, 0, 0, 0, 0, t4107 * t4340 + t4109 * t4339 + t4111 * t4338 + t4114 * t4430 + t4116 * t4429 + t4118 * t4428, t4108 * t4340 + t4110 * t4339 + t4112 * t4338 + t4113 * t4430 + t4115 * t4429 + t4117 * t4428, 0, 0, 0, 0, 0, -t4104 * t4451 - t4105 * t4450 - t4106 * t4449 + t4099 * t4430 + t4101 * t4429 + t4103 * t4428 + (t4080 * t4343 + t4084 * t4342 + t4085 * t4341) * t4276, t4337 * t4451 + t4336 * t4450 + t4335 * t4449 + t4098 * t4430 + t4100 * t4429 + t4102 * t4428 + (t4081 * t4343 + t4082 * t4342 + t4083 * t4341) * t4276, -g(3);];
tau_reg  = t1;
