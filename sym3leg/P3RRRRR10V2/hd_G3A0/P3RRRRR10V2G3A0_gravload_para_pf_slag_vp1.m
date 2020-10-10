% Calculate Gravitation load for parallel robot
% P3RRRRR10V2G3A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m [4x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
%
% Output:
% taugX [3x1]
%   forces required to compensate gravitation load
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 03:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1(xP, qJ, g, legFrame, ...
  koppelP, pkin, m, rSges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: g has to be [3x1] (double)');
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR10V2G3A0_gravload_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From gravvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 02:22:16
% EndTime: 2020-08-07 02:22:29
% DurationCPUTime: 13.39s
% Computational Cost: add. (2379->590), mult. (4518->1008), div. (36->13), fcn. (3144->26), ass. (0->384)
t4255 = legFrame(3,2);
t4223 = sin(t4255);
t4226 = cos(t4255);
t4168 = g(1) * t4226 - g(2) * t4223;
t4260 = sin(qJ(1,3));
t4269 = cos(qJ(1,3));
t4137 = -g(3) * t4260 + t4168 * t4269;
t4259 = sin(qJ(2,3));
t4268 = cos(qJ(2,3));
t4140 = g(3) * t4269 + t4168 * t4260;
t4165 = g(1) * t4223 + g(2) * t4226;
t4253 = sin(pkin(4));
t4254 = cos(pkin(4));
t4568 = -t4140 * t4254 + t4165 * t4253;
t4288 = t4268 * t4137 + t4568 * t4259;
t4256 = legFrame(2,2);
t4224 = sin(t4256);
t4227 = cos(t4256);
t4169 = g(1) * t4227 - g(2) * t4224;
t4263 = sin(qJ(1,2));
t4272 = cos(qJ(1,2));
t4138 = -g(3) * t4263 + t4169 * t4272;
t4262 = sin(qJ(2,2));
t4271 = cos(qJ(2,2));
t4141 = g(3) * t4272 + t4169 * t4263;
t4166 = g(1) * t4224 + g(2) * t4227;
t4569 = -t4141 * t4254 + t4166 * t4253;
t4287 = t4271 * t4138 + t4569 * t4262;
t4257 = legFrame(1,2);
t4225 = sin(t4257);
t4228 = cos(t4257);
t4170 = g(1) * t4228 - g(2) * t4225;
t4266 = sin(qJ(1,1));
t4275 = cos(qJ(1,1));
t4139 = -g(3) * t4266 + t4170 * t4275;
t4265 = sin(qJ(2,1));
t4274 = cos(qJ(2,1));
t4142 = g(3) * t4275 + t4170 * t4266;
t4167 = g(1) * t4225 + g(2) * t4228;
t4570 = -t4142 * t4254 + t4167 * t4253;
t4286 = t4274 * t4139 + t4570 * t4265;
t4267 = cos(qJ(3,3));
t4232 = t4267 * pkin(3);
t4205 = t4232 + pkin(2);
t4576 = t4205 * t4259;
t4270 = cos(qJ(3,2));
t4233 = t4270 * pkin(3);
t4206 = t4233 + pkin(2);
t4575 = t4206 * t4262;
t4273 = cos(qJ(3,1));
t4234 = t4273 * pkin(3);
t4207 = t4234 + pkin(2);
t4574 = t4207 * t4265;
t4405 = m(1) * rSges(1,1) + (m(2) + m(3)) * pkin(1);
t4248 = t4268 ^ 2;
t4277 = pkin(8) + pkin(7);
t4573 = t4248 * t4277;
t4250 = t4271 ^ 2;
t4572 = t4250 * t4277;
t4252 = t4274 ^ 2;
t4571 = t4252 * t4277;
t4216 = t4277 * t4274;
t4176 = -pkin(2) * t4265 + t4216;
t4310 = t4253 * t4176;
t4215 = t4277 * t4271;
t4175 = -pkin(2) * t4262 + t4215;
t4311 = t4253 * t4175;
t4214 = t4277 * t4268;
t4174 = -pkin(2) * t4259 + t4214;
t4312 = t4253 * t4174;
t4567 = t4140 * t4253 + t4165 * t4254;
t4566 = t4141 * t4253 + t4166 * t4254;
t4565 = t4142 * t4253 + t4167 * t4254;
t4258 = sin(qJ(3,3));
t4404 = pkin(6) * t4258 + pkin(3);
t4247 = t4267 ^ 2;
t4544 = t4247 * pkin(3);
t4547 = pkin(2) * t4267;
t4564 = (-t4404 + 0.2e1 * t4544 + t4547) * t4259 - t4267 * t4214;
t4261 = sin(qJ(3,2));
t4403 = pkin(6) * t4261 + pkin(3);
t4249 = t4270 ^ 2;
t4543 = t4249 * pkin(3);
t4546 = pkin(2) * t4270;
t4563 = (-t4403 + 0.2e1 * t4543 + t4546) * t4262 - t4270 * t4215;
t4264 = sin(qJ(3,1));
t4402 = pkin(6) * t4264 + pkin(3);
t4251 = t4273 ^ 2;
t4542 = t4251 * pkin(3);
t4545 = pkin(2) * t4273;
t4562 = (-t4402 + 0.2e1 * t4542 + t4545) * t4265 - t4273 * t4216;
t4213 = t4277 * t4265;
t4237 = pkin(2) * t4274;
t4430 = t4213 + t4237;
t4173 = pkin(1) + t4430;
t4231 = t4264 * pkin(3);
t4349 = (t4252 - 0.2e1) * t4231 - pkin(6);
t4480 = t4253 * t4266;
t4561 = -t4173 * t4275 + t4349 * t4480;
t4212 = t4277 * t4262;
t4236 = pkin(2) * t4271;
t4431 = t4212 + t4236;
t4172 = pkin(1) + t4431;
t4230 = t4261 * pkin(3);
t4350 = (t4250 - 0.2e1) * t4230 - pkin(6);
t4482 = t4253 * t4263;
t4560 = -t4172 * t4272 + t4350 * t4482;
t4211 = t4277 * t4259;
t4235 = pkin(2) * t4268;
t4432 = t4211 + t4235;
t4171 = pkin(1) + t4432;
t4229 = t4258 * pkin(3);
t4351 = (t4248 - 0.2e1) * t4229 - pkin(6);
t4484 = t4253 * t4260;
t4559 = -t4171 * t4269 + t4351 * t4484;
t4557 = m(1) * rSges(1,2);
t4556 = t4254 / 0.2e1;
t4555 = m(3) / pkin(3);
t4554 = pkin(1) * t4254;
t4238 = pkin(1) * t4259;
t4239 = pkin(1) * t4262;
t4240 = pkin(1) * t4265;
t4553 = pkin(2) * t4248;
t4552 = pkin(2) * t4250;
t4551 = pkin(2) * t4252;
t4220 = t4254 * pkin(2);
t4541 = t4258 * pkin(2);
t4540 = t4260 * pkin(2);
t4539 = t4261 * pkin(2);
t4538 = t4263 * pkin(2);
t4537 = t4264 * pkin(2);
t4536 = t4266 * pkin(2);
t4245 = pkin(2) - t4277;
t4244 = pkin(2) + t4277;
t4200 = t4229 - pkin(6);
t4427 = m(2) * rSges(2,1) + pkin(2) * m(3);
t4149 = (rSges(3,1) * t4267 - rSges(3,2) * t4258) * m(3) + t4427;
t4189 = (-pkin(7) - rSges(3,3)) * m(3) + m(2) * rSges(2,2);
t4429 = m(2) * (rSges(2,3) + pkin(6)) + pkin(6) * m(3);
t4303 = (t4149 * t4259 + t4189 * t4268) * t4254 - ((rSges(3,1) * t4258 + rSges(3,2) * t4267) * m(3) + t4429) * t4253 + t4557;
t4315 = t4149 * t4268 - t4189 * t4259;
t4408 = t4171 * t4541;
t4417 = t4254 * t4238;
t4425 = pkin(1) * t4229;
t4433 = t4405 * g(3);
t4479 = t4253 * t4268;
t4535 = (t4433 * t4269 + (-t4303 * t4260 + t4315 * t4269) * g(3) + (t4303 * t4269 + (t4315 + t4405) * t4260) * t4168) / (-(pkin(6) * t4479 + t4417) * t4544 + ((-pkin(6) * t4211 + t4200 * t4235 + t4425) * t4253 + t4174 * t4554) * t4267 + t4253 * t4408);
t4202 = t4230 - pkin(6);
t4150 = (rSges(3,1) * t4270 - rSges(3,2) * t4261) * m(3) + t4427;
t4302 = (t4150 * t4262 + t4189 * t4271) * t4254 - ((rSges(3,1) * t4261 + rSges(3,2) * t4270) * m(3) + t4429) * t4253 + t4557;
t4314 = t4150 * t4271 - t4189 * t4262;
t4407 = t4172 * t4539;
t4416 = t4254 * t4239;
t4424 = pkin(1) * t4230;
t4477 = t4253 * t4271;
t4534 = (t4433 * t4272 + (-t4302 * t4263 + t4314 * t4272) * g(3) + (t4302 * t4272 + (t4314 + t4405) * t4263) * t4169) / (-(pkin(6) * t4477 + t4416) * t4543 + ((-pkin(6) * t4212 + t4202 * t4236 + t4424) * t4253 + t4175 * t4554) * t4270 + t4253 * t4407);
t4204 = t4231 - pkin(6);
t4151 = (rSges(3,1) * t4273 - rSges(3,2) * t4264) * m(3) + t4427;
t4301 = (t4151 * t4265 + t4189 * t4274) * t4254 - ((rSges(3,1) * t4264 + rSges(3,2) * t4273) * m(3) + t4429) * t4253 + t4557;
t4313 = t4151 * t4274 - t4189 * t4265;
t4406 = t4173 * t4537;
t4415 = t4254 * t4240;
t4423 = pkin(1) * t4231;
t4475 = t4253 * t4274;
t4533 = (t4433 * t4275 + (-t4301 * t4266 + t4313 * t4275) * g(3) + (t4301 * t4275 + (t4313 + t4405) * t4266) * t4170) / (-(pkin(6) * t4475 + t4415) * t4542 + ((-pkin(6) * t4213 + t4204 * t4237 + t4423) * t4253 + t4176 * t4554) * t4273 + t4253 * t4406);
t4208 = t4235 + pkin(1);
t4473 = t4254 * t4267;
t4414 = pkin(1) * t4473;
t4422 = pkin(6) * t4544;
t4532 = (t4288 * t4189 + (t4137 * t4259 - t4268 * t4568) * t4149) / ((-t4268 * t4422 + (-pkin(6) * t4432 + t4208 * t4229) * t4267 + t4408) * t4253 + (t4214 - t4576) * t4414);
t4209 = t4236 + pkin(1);
t4471 = t4254 * t4270;
t4413 = pkin(1) * t4471;
t4421 = pkin(6) * t4543;
t4531 = (t4287 * t4189 + (t4138 * t4262 - t4271 * t4569) * t4150) / ((-t4271 * t4421 + (-pkin(6) * t4431 + t4209 * t4230) * t4270 + t4407) * t4253 + (t4215 - t4575) * t4413);
t4210 = t4237 + pkin(1);
t4469 = t4254 * t4273;
t4412 = pkin(1) * t4469;
t4420 = pkin(6) * t4542;
t4530 = (t4286 * t4189 + (t4139 * t4265 - t4274 * t4570) * t4151) / ((-t4274 * t4420 + (-pkin(6) * t4430 + t4210 * t4231) * t4273 + t4406) * t4253 + (t4216 - t4574) * t4412);
t4517 = (t4232 + t4244) * (t4232 + t4245);
t4516 = (t4233 + t4244) * (t4233 + t4245);
t4515 = (t4234 + t4244) * (t4234 + t4245);
t4467 = t4254 * t4277;
t4194 = pkin(1) * t4467;
t4514 = t4194 * t4267;
t4513 = t4194 * t4270;
t4512 = t4194 * t4273;
t4199 = t4229 + pkin(6);
t4511 = t4199 * t4253;
t4201 = t4230 + pkin(6);
t4510 = t4201 * t4253;
t4203 = t4231 + pkin(6);
t4509 = t4203 * t4253;
t4508 = t4205 * t4223;
t4507 = t4205 * t4226;
t4506 = t4206 * t4224;
t4505 = t4206 * t4227;
t4504 = t4207 * t4225;
t4503 = t4207 * t4228;
t4502 = (t4254 + 0.1e1) * (t4254 - 0.1e1);
t4195 = t4238 + t4277;
t4501 = t4223 * t4195;
t4500 = t4223 * t4260;
t4196 = t4239 + t4277;
t4499 = t4224 * t4196;
t4498 = t4224 * t4263;
t4197 = t4240 + t4277;
t4497 = t4225 * t4197;
t4496 = t4225 * t4266;
t4495 = t4226 * t4259;
t4494 = t4226 * t4260;
t4493 = t4227 * t4262;
t4492 = t4227 * t4263;
t4491 = t4228 * t4265;
t4490 = t4228 * t4266;
t4489 = t4244 * t4245;
t4280 = pkin(3) ^ 2;
t4488 = t4247 * t4280;
t4487 = t4249 * t4280;
t4486 = t4251 * t4280;
t4485 = t4253 * t4259;
t4483 = t4253 * t4262;
t4481 = t4253 * t4265;
t4478 = t4253 * t4269;
t4476 = t4253 * t4272;
t4474 = t4253 * t4275;
t4472 = t4254 * t4269;
t4470 = t4254 * t4272;
t4468 = t4254 * t4275;
t4466 = t4258 * t4260;
t4465 = t4259 * t4268;
t4177 = pkin(1) + 0.2e1 * t4211;
t4464 = t4260 * t4177;
t4463 = t4260 * t4259;
t4462 = t4260 * t4268;
t4461 = t4260 * t4277;
t4460 = t4261 * t4263;
t4459 = t4262 * t4271;
t4178 = pkin(1) + 0.2e1 * t4212;
t4458 = t4263 * t4178;
t4457 = t4263 * t4262;
t4456 = t4263 * t4271;
t4455 = t4263 * t4277;
t4454 = t4264 * t4266;
t4453 = t4265 * t4274;
t4179 = pkin(1) + 0.2e1 * t4213;
t4452 = t4266 * t4179;
t4451 = t4266 * t4265;
t4450 = t4266 * t4274;
t4449 = t4266 * t4277;
t4447 = t4268 * t4269;
t4446 = t4269 * t4226;
t4445 = t4269 * t4259;
t4443 = t4271 * t4272;
t4442 = t4272 * t4227;
t4441 = t4272 * t4262;
t4439 = t4274 * t4275;
t4438 = t4275 * t4228;
t4437 = t4275 * t4265;
t4436 = t4277 * t4269;
t4435 = t4277 * t4272;
t4434 = t4277 * t4275;
t4243 = t4254 ^ 2;
t4428 = t4243 - 0.1e1 / 0.2e1;
t4426 = 0.2e1 * t4220;
t4419 = -0.2e1 * t4502;
t4411 = pkin(3) * t4466;
t4410 = pkin(3) * t4460;
t4409 = pkin(3) * t4454;
t4372 = t4254 * t4461;
t4401 = (0.2e1 * t4205 * t4372 + t4269 * t4517) * t4253 * t4248;
t4370 = t4254 * t4455;
t4400 = (0.2e1 * t4206 * t4370 + t4272 * t4516) * t4253 * t4250;
t4368 = t4254 * t4449;
t4399 = (0.2e1 * t4207 * t4368 + t4275 * t4515) * t4253 * t4252;
t4398 = t4177 * t4205 * t4269;
t4397 = t4178 * t4206 * t4272;
t4396 = t4179 * t4207 * t4275;
t4395 = t4259 * t4517;
t4394 = t4262 * t4516;
t4393 = t4265 * t4515;
t4392 = t4205 * t4500;
t4391 = t4205 * t4494;
t4390 = t4206 * t4498;
t4389 = t4206 * t4492;
t4388 = t4207 * t4496;
t4387 = t4207 * t4490;
t4386 = t4223 * t4485;
t4385 = t4223 * t4463;
t4384 = t4224 * t4483;
t4383 = t4224 * t4457;
t4382 = t4225 * t4481;
t4381 = t4225 * t4451;
t4380 = t4226 * t4463;
t4379 = t4227 * t4457;
t4378 = t4228 * t4451;
t4377 = t4254 * t4489;
t4376 = t4258 * t4478;
t4375 = t4261 * t4476;
t4374 = t4264 * t4474;
t4373 = t4254 * t4445;
t4371 = t4254 * t4441;
t4369 = t4254 * t4437;
t4367 = t4226 * t4485;
t4366 = t4227 * t4483;
t4365 = t4228 * t4481;
t4361 = -0.2e1 * pkin(2) * t4467;
t4113 = (-t4567 * rSges(3,1) + t4288 * rSges(3,2)) * t4267 + t4258 * (t4288 * rSges(3,1) + t4567 * rSges(3,2));
t4180 = pkin(6) * t4267 - t4541;
t4360 = t4113 / ((-t4180 * t4205 * t4253 + t4514) * t4268 + (pkin(1) * t4205 * t4258 - t4180 * t4211) * t4253 - t4414 * t4576) * t4555;
t4114 = (-t4566 * rSges(3,1) + t4287 * rSges(3,2)) * t4270 + t4261 * (t4287 * rSges(3,1) + t4566 * rSges(3,2));
t4181 = pkin(6) * t4270 - t4539;
t4359 = t4114 / ((-t4181 * t4206 * t4253 + t4513) * t4271 + (pkin(1) * t4206 * t4261 - t4181 * t4212) * t4253 - t4413 * t4575) * t4555;
t4115 = (-t4565 * rSges(3,1) + t4286 * rSges(3,2)) * t4273 + t4264 * (t4286 * rSges(3,1) + t4565 * rSges(3,2));
t4182 = pkin(6) * t4273 - t4537;
t4358 = t4115 / ((-t4182 * t4207 * t4253 + t4512) * t4274 + (pkin(1) * t4207 * t4264 - t4182 * t4213) * t4253 - t4412 * t4574) * t4555;
t4357 = pkin(3) * t4385;
t4356 = pkin(3) * t4380;
t4355 = pkin(3) * t4383;
t4354 = pkin(3) * t4379;
t4353 = pkin(3) * t4381;
t4352 = pkin(3) * t4378;
t4348 = t4260 * t4395;
t4347 = t4263 * t4394;
t4346 = t4266 * t4393;
t4345 = t4199 * t4385;
t4344 = t4199 * t4380;
t4343 = t4201 * t4383;
t4342 = t4201 * t4379;
t4341 = t4203 * t4381;
t4340 = t4203 * t4378;
t4339 = t4445 * t4479;
t4338 = t4441 * t4477;
t4337 = t4437 * t4475;
t4336 = t4199 * t4392;
t4335 = t4199 * t4391;
t4334 = t4201 * t4390;
t4333 = t4201 * t4389;
t4332 = t4203 * t4388;
t4331 = t4203 * t4387;
t4330 = t4259 * t4404;
t4329 = t4262 * t4403;
t4328 = t4265 * t4402;
t4327 = t4205 * t4248 * t4419;
t4326 = t4206 * t4250 * t4419;
t4325 = t4207 * t4252 * t4419;
t4190 = t4211 + pkin(1);
t4324 = t4190 * t4268 + t4553;
t4191 = t4212 + pkin(1);
t4323 = t4191 * t4271 + t4552;
t4192 = t4213 + pkin(1);
t4322 = t4192 * t4274 + t4551;
t4309 = t4259 * t4214 - pkin(2) + t4553;
t4308 = t4262 * t4215 - pkin(2) + t4552;
t4307 = t4265 * t4216 - pkin(2) + t4551;
t4306 = pkin(6) * t4463 + (pkin(2) * t4465 + t4195 - t4573) * t4478;
t4305 = pkin(6) * t4457 + (pkin(2) * t4459 + t4196 - t4572) * t4476;
t4304 = pkin(6) * t4451 + (pkin(2) * t4453 + t4197 - t4571) * t4474;
t4297 = t4309 * t4466;
t4296 = t4308 * t4460;
t4295 = t4307 * t4454;
t4294 = (-t4254 * t4540 + t4436) * t4259 + (pkin(2) * t4269 + t4372) * t4268 + t4253 * t4411;
t4293 = (-t4254 * t4538 + t4435) * t4262 + (pkin(2) * t4272 + t4370) * t4271 + t4253 * t4410;
t4292 = (-t4254 * t4536 + t4434) * t4265 + (pkin(2) * t4275 + t4368) * t4274 + t4253 * t4409;
t4291 = -t4254 * t4229 + t4312;
t4290 = -t4254 * t4230 + t4311;
t4289 = -t4254 * t4231 + t4310;
t4285 = -t4309 * t4258 - t4351 * t4267;
t4284 = -t4308 * t4261 - t4350 * t4270;
t4283 = -t4307 * t4264 - t4349 * t4273;
t4282 = pkin(2) ^ 2;
t4242 = pkin(1) * t4277;
t4164 = -pkin(6) * t4253 * t4277 - pkin(1) * t4220;
t4163 = t4254 * t4434 - t4536;
t4162 = t4254 * t4435 - t4538;
t4161 = t4254 * t4436 - t4540;
t4154 = -t4254 * t4451 + t4439;
t4153 = -t4254 * t4457 + t4443;
t4152 = -t4254 * t4463 + t4447;
t4136 = t4242 - t4393;
t4135 = t4242 - t4394;
t4134 = t4242 - t4395;
t1 = [(-(t4154 * t4228 + t4382) * t4542 + (t4289 * t4225 - t4292 * t4228) * t4273 - (t4225 * t4254 + t4228 * t4480) * t4537) * t4533 + ((t4283 * t4225 - t4562 * t4490) * t4243 + ((t4274 * t4438 + 0.2e1 * t4382) * t4542 + (-t4225 * t4310 - t4561 * t4228) * t4273 - t4253 * (t4225 * t4328 + t4228 * t4295)) * t4254 + t4251 * t4352 + (-(-t4225 * t4252 + t4228 * t4337 + t4225) * t4231 - t4225 * pkin(6)) * t4273 + (t4322 * t4225 - t4304 * t4228) * t4264 - t4352) * t4530 + (-t4228 * t4399 + (((-t4203 * t4504 + t4228 * t4346) * t4254 - t4228 * t4396) * t4253 + (t4225 * t4393 + t4331) * t4243 - t4331 + t4225 * t4136) * t4274 - t4207 * t4497 + (t4225 * t4325 - ((t4203 * t4225 * t4265 - t4387) * t4254 + t4197 * t4438) * t4253 + (t4340 + t4504) * t4243 - t4340) * t4277) * t4358 + (-(t4153 * t4227 + t4384) * t4543 + (t4290 * t4224 - t4293 * t4227) * t4270 - (t4224 * t4254 + t4227 * t4482) * t4539) * t4534 + ((t4284 * t4224 - t4563 * t4492) * t4243 + ((t4271 * t4442 + 0.2e1 * t4384) * t4543 + (-t4224 * t4311 - t4560 * t4227) * t4270 - t4253 * (t4224 * t4329 + t4227 * t4296)) * t4254 + t4249 * t4354 + (-(-t4224 * t4250 + t4227 * t4338 + t4224) * t4230 - t4224 * pkin(6)) * t4270 + (t4323 * t4224 - t4305 * t4227) * t4261 - t4354) * t4531 + (-t4227 * t4400 + (((-t4201 * t4506 + t4227 * t4347) * t4254 - t4227 * t4397) * t4253 + (t4224 * t4394 + t4333) * t4243 - t4333 + t4224 * t4135) * t4271 - t4206 * t4499 + (t4224 * t4326 - ((t4201 * t4224 * t4262 - t4389) * t4254 + t4196 * t4442) * t4253 + (t4342 + t4506) * t4243 - t4342) * t4277) * t4359 + (-(t4152 * t4226 + t4386) * t4544 + (t4291 * t4223 - t4294 * t4226) * t4267 - (t4223 * t4254 + t4226 * t4484) * t4541) * t4535 + ((t4285 * t4223 - t4564 * t4494) * t4243 + ((t4268 * t4446 + 0.2e1 * t4386) * t4544 + (-t4223 * t4312 - t4559 * t4226) * t4267 - t4253 * (t4223 * t4330 + t4226 * t4297)) * t4254 + t4247 * t4356 + (-(-t4223 * t4248 + t4226 * t4339 + t4223) * t4229 - t4223 * pkin(6)) * t4267 + (t4324 * t4223 - t4306 * t4226) * t4258 - t4356) * t4532 + (-t4226 * t4401 + (((-t4199 * t4508 + t4226 * t4348) * t4254 - t4226 * t4398) * t4253 + (t4223 * t4395 + t4335) * t4243 - t4335 + t4223 * t4134) * t4268 - t4205 * t4501 + (t4223 * t4327 - ((t4199 * t4223 * t4259 - t4391) * t4254 + t4195 * t4446) * t4253 + (t4344 + t4508) * t4243 - t4344) * t4277) * t4360 - m(4) * g(1); ((t4154 * t4225 - t4365) * t4542 + (t4225 * t4292 + t4228 * t4289) * t4273 + (t4225 * t4480 - t4228 * t4254) * t4537) * t4533 + ((t4283 * t4228 + t4562 * t4496) * t4243 + (-(t4225 * t4439 - 0.2e1 * t4365) * t4542 + (t4561 * t4225 - t4228 * t4310) * t4273 + t4253 * (t4225 * t4295 - t4402 * t4491)) * t4254 - t4251 * t4353 + ((t4225 * t4337 + t4228 * t4252 - t4228) * t4231 - t4228 * pkin(6)) * t4273 + (t4225 * t4304 + t4228 * t4322) * t4264 + t4353) * t4530 + (t4225 * t4399 + (((-t4203 * t4503 - t4225 * t4346) * t4254 + t4225 * t4396) * t4253 + (t4228 * t4393 - t4332) * t4243 + t4332 + t4228 * t4136) * t4274 - t4197 * t4503 + (t4228 * t4325 + ((-t4203 * t4491 - t4388) * t4254 + t4275 * t4497) * t4253 + (-t4341 + t4503) * t4243 + t4341) * t4277) * t4358 + ((t4153 * t4224 - t4366) * t4543 + (t4224 * t4293 + t4227 * t4290) * t4270 + (t4224 * t4482 - t4227 * t4254) * t4539) * t4534 + ((t4284 * t4227 + t4563 * t4498) * t4243 + (-(t4224 * t4443 - 0.2e1 * t4366) * t4543 + (t4560 * t4224 - t4227 * t4311) * t4270 + t4253 * (t4224 * t4296 - t4403 * t4493)) * t4254 - t4249 * t4355 + ((t4224 * t4338 + t4227 * t4250 - t4227) * t4230 - t4227 * pkin(6)) * t4270 + (t4224 * t4305 + t4227 * t4323) * t4261 + t4355) * t4531 + (t4224 * t4400 + (((-t4201 * t4505 - t4224 * t4347) * t4254 + t4224 * t4397) * t4253 + (t4227 * t4394 - t4334) * t4243 + t4334 + t4227 * t4135) * t4271 - t4196 * t4505 + (t4227 * t4326 + ((-t4201 * t4493 - t4390) * t4254 + t4272 * t4499) * t4253 + (-t4343 + t4505) * t4243 + t4343) * t4277) * t4359 + ((t4152 * t4223 - t4367) * t4544 + (t4223 * t4294 + t4226 * t4291) * t4267 + (t4223 * t4484 - t4226 * t4254) * t4541) * t4535 + ((t4285 * t4226 + t4564 * t4500) * t4243 + (-(t4223 * t4447 - 0.2e1 * t4367) * t4544 + (t4559 * t4223 - t4226 * t4312) * t4267 + t4253 * (t4223 * t4297 - t4404 * t4495)) * t4254 - t4247 * t4357 + ((t4223 * t4339 + t4226 * t4248 - t4226) * t4229 - t4226 * pkin(6)) * t4267 + (t4223 * t4306 + t4226 * t4324) * t4258 + t4357) * t4532 + (t4223 * t4401 + (((-t4199 * t4507 - t4223 * t4348) * t4254 + t4223 * t4398) * t4253 + (t4226 * t4395 - t4336) * t4243 + t4336 + t4226 * t4134) * t4268 - t4195 * t4507 + (t4226 * t4327 + ((-t4199 * t4495 - t4392) * t4254 + t4269 * t4501) * t4253 + (-t4345 + t4507) * t4243 + t4345) * t4277) * t4360 - m(4) * g(2); ((t4369 + t4450) * t4542 + (-pkin(3) * t4374 - t4163 * t4274 + (pkin(2) * t4468 + t4449) * t4265) * t4273 - pkin(2) * t4374) * t4533 + (((-t4349 * t4468 + t4409 * t4453) * t4273 + (-t4307 * t4468 + t4266 * (t4210 * t4265 + t4277 - t4571)) * t4264) * t4253 - 0.2e1 * (t4428 * t4437 + t4450 * t4556) * t4542 - (t4173 * t4266 - t4176 * t4468) * t4469 + t4275 * t4328 * t4502) * t4530 + ((t4371 + t4456) * t4543 + (-pkin(3) * t4375 - t4162 * t4271 + (pkin(2) * t4470 + t4455) * t4262) * t4270 - pkin(2) * t4375) * t4534 + (((-t4350 * t4470 + t4410 * t4459) * t4270 + (-t4308 * t4470 + t4263 * (t4209 * t4262 + t4277 - t4572)) * t4261) * t4253 - 0.2e1 * (t4428 * t4441 + t4456 * t4556) * t4543 - (t4172 * t4263 - t4175 * t4470) * t4471 + t4272 * t4329 * t4502) * t4531 + ((t4373 + t4462) * t4544 + (-pkin(3) * t4376 - t4161 * t4268 + (pkin(2) * t4472 + t4461) * t4259) * t4267 - pkin(2) * t4376) * t4535 + (((-t4351 * t4472 + t4411 * t4465) * t4267 + (-t4309 * t4472 + t4260 * (t4208 * t4259 + t4277 - t4573)) * t4258) * t4253 - 0.2e1 * (t4428 * t4445 + t4462 * t4556) * t4544 - (t4171 * t4260 - t4174 * t4472) * t4473 + t4269 * t4330 * t4502) * t4532 - m(4) * g(3) + (((t4275 * t4361 - 0.2e1 * t4163 * t4234 + (t4486 + t4489) * t4266) * t4252 + (t4369 * t4486 + ((t4265 * t4426 - t4509) * t4275 + t4452) * t4234 + (-pkin(2) * t4509 + t4265 * t4377) * t4275 + pkin(2) * t4452) * t4274 + t4277 * (t4266 * t4197 + (pkin(3) * t4469 - t4203 * t4481 + t4220) * t4275)) / ((t4512 + (t4204 * t4545 + t4264 * t4282 - t4420) * t4253) * t4274 - t4415 * t4542 + (t4164 * t4265 + t4253 * t4423) * t4273 + t4192 * t4253 * t4537) * t4115 + ((t4272 * t4361 - 0.2e1 * t4162 * t4233 + (t4487 + t4489) * t4263) * t4250 + (t4371 * t4487 + ((t4262 * t4426 - t4510) * t4272 + t4458) * t4233 + (-pkin(2) * t4510 + t4262 * t4377) * t4272 + pkin(2) * t4458) * t4271 + t4277 * (t4263 * t4196 + (pkin(3) * t4471 - t4201 * t4483 + t4220) * t4272)) / ((t4513 + (t4202 * t4546 + t4261 * t4282 - t4421) * t4253) * t4271 - t4416 * t4543 + (t4164 * t4262 + t4253 * t4424) * t4270 + t4191 * t4253 * t4539) * t4114 + ((t4269 * t4361 - 0.2e1 * t4161 * t4232 + (t4488 + t4489) * t4260) * t4248 + (t4373 * t4488 + ((t4259 * t4426 - t4511) * t4269 + t4464) * t4232 + (-pkin(2) * t4511 + t4259 * t4377) * t4269 + pkin(2) * t4464) * t4268 + t4277 * (t4260 * t4195 + (pkin(3) * t4473 - t4199 * t4485 + t4220) * t4269)) / ((t4514 + (t4200 * t4547 + t4258 * t4282 - t4422) * t4253) * t4268 - t4417 * t4544 + (t4164 * t4259 + t4253 * t4425) * t4267 + t4190 * t4253 * t4541) * t4113) * t4253 * t4555;];
taugX  = t1;
