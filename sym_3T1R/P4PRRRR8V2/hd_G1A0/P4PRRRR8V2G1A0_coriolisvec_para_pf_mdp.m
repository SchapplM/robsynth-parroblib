% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P4PRRRR8V2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% xDP [4x1]
%   Generalized platform velocities
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4PRRRR8V2G1A0_convert_par2_MPV_fixb.m

% Output:
% taucX [4x1]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:13
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(8,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [4 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: xDP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: pkin has to be [8x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P4PRRRR8V2G1A0_coriolisvec_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 11:12:36
% EndTime: 2020-08-07 11:13:05
% DurationCPUTime: 30.77s
% Computational Cost: add. (172836->854), mult. (352110->1660), div. (11204->22), fcn. (325172->30), ass. (0->631)
t4364 = cos(qJ(3,4));
t4318 = pkin(3) * t4364 + pkin(2);
t4365 = cos(qJ(2,4));
t4381 = pkin(7) + pkin(6);
t4319 = t4365 * t4381;
t4363 = sin(qJ(2,4));
t4259 = t4318 * t4363 - t4319;
t4357 = cos(pkin(4));
t4362 = sin(qJ(3,4));
t4635 = t4357 * t4362;
t4282 = t4318 * t4635;
t4355 = sin(pkin(4));
t4650 = t4355 * t4364;
t4213 = t4259 * t4650 + t4282;
t4210 = 0.1e1 / t4213;
t4372 = cos(qJ(3,3));
t4320 = pkin(3) * t4372 + pkin(2);
t4373 = cos(qJ(2,3));
t4323 = t4373 * t4381;
t4367 = sin(qJ(2,3));
t4263 = t4320 * t4367 - t4323;
t4366 = sin(qJ(3,3));
t4633 = t4357 * t4366;
t4283 = t4320 * t4633;
t4642 = t4355 * t4372;
t4223 = t4263 * t4642 + t4283;
t4214 = 0.1e1 / t4223;
t4374 = cos(qJ(3,2));
t4321 = pkin(3) * t4374 + pkin(2);
t4375 = cos(qJ(2,2));
t4324 = t4375 * t4381;
t4369 = sin(qJ(2,2));
t4264 = t4321 * t4369 - t4324;
t4368 = sin(qJ(3,2));
t4631 = t4357 * t4368;
t4284 = t4321 * t4631;
t4640 = t4355 * t4374;
t4224 = t4264 * t4640 + t4284;
t4216 = 0.1e1 / t4224;
t4376 = cos(qJ(3,1));
t4322 = pkin(3) * t4376 + pkin(2);
t4377 = cos(qJ(2,1));
t4325 = t4377 * t4381;
t4371 = sin(qJ(2,1));
t4265 = t4322 * t4371 - t4325;
t4370 = sin(qJ(3,1));
t4629 = t4357 * t4370;
t4285 = t4322 * t4629;
t4638 = t4355 * t4376;
t4225 = t4265 * t4638 + t4285;
t4218 = 0.1e1 / t4225;
t4294 = pkin(2) * t4363 - t4319;
t4651 = t4355 * t4363;
t4346 = t4364 ^ 2;
t4739 = pkin(3) * t4346;
t4186 = 0.1e1 / ((pkin(3) * t4635 + t4294 * t4355) * t4364 + pkin(2) * t4635 + t4651 * t4739);
t4296 = pkin(2) * t4367 - t4323;
t4647 = t4355 * t4367;
t4350 = t4372 ^ 2;
t4738 = pkin(3) * t4350;
t4187 = 0.1e1 / ((pkin(3) * t4633 + t4296 * t4355) * t4372 + pkin(2) * t4633 + t4647 * t4738);
t4297 = pkin(2) * t4369 - t4324;
t4645 = t4355 * t4369;
t4351 = t4374 ^ 2;
t4737 = pkin(3) * t4351;
t4188 = 0.1e1 / ((pkin(3) * t4631 + t4297 * t4355) * t4374 + pkin(2) * t4631 + t4645 * t4737);
t4298 = pkin(2) * t4371 - t4325;
t4643 = t4355 * t4371;
t4352 = t4376 ^ 2;
t4736 = pkin(3) * t4352;
t4189 = 0.1e1 / ((pkin(3) * t4629 + t4298 * t4355) * t4376 + pkin(2) * t4629 + t4643 * t4736);
t4382 = xP(4);
t4343 = sin(t4382);
t4344 = cos(t4382);
t4383 = koppelP(4,2);
t4387 = koppelP(4,1);
t4278 = -t4343 * t4383 + t4344 * t4387;
t4378 = xDP(4);
t4379 = xDP(2);
t4230 = -t4278 * t4378 - t4379;
t4274 = t4343 * t4387 + t4344 * t4383;
t4380 = xDP(1);
t4233 = t4274 * t4378 - t4380;
t4358 = legFrame(4,3);
t4331 = sin(t4358);
t4335 = cos(t4358);
t4354 = sin(pkin(8));
t4356 = cos(pkin(8));
t4242 = -t4331 * t4354 + t4335 * t4356;
t4246 = t4331 * t4356 + t4335 * t4354;
t4634 = t4357 * t4363;
t4404 = t4242 * t4634 + t4246 * t4365;
t4405 = -t4242 * t4365 + t4246 * t4634;
t4131 = ((t4230 * t4246 + t4233 * t4242) * t4650 + (t4230 * t4405 + t4233 * t4404) * t4362) * t4210;
t4775 = 0.2e1 * t4131;
t4385 = koppelP(2,2);
t4389 = koppelP(2,1);
t4280 = -t4343 * t4385 + t4344 * t4389;
t4231 = -t4280 * t4378 - t4379;
t4276 = t4343 * t4389 + t4344 * t4385;
t4234 = t4276 * t4378 - t4380;
t4360 = legFrame(2,3);
t4333 = sin(t4360);
t4337 = cos(t4360);
t4244 = -t4333 * t4354 + t4337 * t4356;
t4248 = t4333 * t4356 + t4337 * t4354;
t4630 = t4357 * t4369;
t4400 = t4244 * t4630 + t4248 * t4375;
t4401 = -t4244 * t4375 + t4248 * t4630;
t4136 = ((t4231 * t4248 + t4234 * t4244) * t4640 + (t4231 * t4401 + t4234 * t4400) * t4368) * t4216;
t4774 = 0.2e1 * t4136;
t4386 = koppelP(1,2);
t4390 = koppelP(1,1);
t4281 = -t4343 * t4386 + t4344 * t4390;
t4232 = -t4281 * t4378 - t4379;
t4277 = t4343 * t4390 + t4344 * t4386;
t4235 = t4277 * t4378 - t4380;
t4361 = legFrame(1,3);
t4334 = sin(t4361);
t4338 = cos(t4361);
t4245 = -t4334 * t4354 + t4338 * t4356;
t4249 = t4334 * t4356 + t4338 * t4354;
t4628 = t4357 * t4371;
t4398 = t4245 * t4628 + t4249 * t4377;
t4399 = -t4245 * t4377 + t4249 * t4628;
t4137 = ((t4232 * t4249 + t4235 * t4245) * t4638 + (t4232 * t4399 + t4235 * t4398) * t4370) * t4218;
t4773 = 0.2e1 * t4137;
t4384 = koppelP(3,2);
t4388 = koppelP(3,1);
t4279 = -t4343 * t4384 + t4344 * t4388;
t4236 = -t4279 * t4378 - t4379;
t4275 = t4343 * t4388 + t4344 * t4384;
t4237 = t4275 * t4378 - t4380;
t4359 = legFrame(3,3);
t4332 = sin(t4359);
t4336 = cos(t4359);
t4243 = -t4332 * t4354 + t4336 * t4356;
t4247 = t4332 * t4356 + t4336 * t4354;
t4632 = t4357 * t4367;
t4402 = t4243 * t4632 + t4247 * t4373;
t4403 = -t4243 * t4373 + t4247 * t4632;
t4138 = ((t4236 * t4247 + t4237 * t4243) * t4642 + (t4236 * t4403 + t4237 * t4402) * t4366) * t4214;
t4772 = 0.2e1 * t4138;
t4703 = t4131 * t4381;
t4552 = t4362 * t4703;
t4182 = t4230 * t4354 + t4233 * t4356;
t4425 = t4230 * t4356 - t4233 * t4354;
t4622 = t4363 * t4381;
t4675 = (t4318 * t4365 + t4622) * t4357;
t4143 = (t4182 * t4335 + t4331 * t4425) * t4675 - t4259 * (t4331 * t4182 - t4425 * t4335);
t4693 = t4143 * t4210;
t4108 = t4552 - t4693;
t4166 = -t4242 * t4650 - t4404 * t4362;
t4167 = -t4246 * t4650 - t4405 * t4362;
t4353 = t4378 ^ 2;
t4392 = 0.1e1 / pkin(3);
t4537 = t4392 * t4693;
t4458 = t4365 * t4537;
t4623 = t4363 * t4364;
t4493 = t4355 * t4623;
t4683 = t4210 * t4362;
t4538 = t4143 * t4683;
t4649 = t4355 * t4365;
t4744 = pkin(2) * t4365;
t4066 = (-((t4131 * t4649 + t4357 * t4537) * t4739 + ((-t4538 + t4703) * t4363 + t4131 * t4744) * t4650 + t4357 * t4108) * t4131 - (t4355 * t4458 + (t4346 * t4357 - t4362 * t4493 - t4357) * t4131) * t4693 + (-t4166 * t4278 - t4167 * t4274) * t4353) * t4186;
t4174 = -t4242 * t4675 + t4246 * t4259;
t4175 = -t4242 * t4259 - t4246 * t4675;
t4238 = pkin(3) * t4623 + t4294;
t4198 = 0.1e1 / (t4238 * t4650 + t4282);
t4326 = pkin(2) ^ 2 + t4381 ^ 2;
t4391 = pkin(3) ^ 2;
t4751 = 0.2e1 * pkin(2);
t4735 = pkin(3) * t4751;
t4557 = (-t4381 * t4538 + (t4346 * t4391 + t4364 * t4735 + t4326) * t4131) * t4131 * t4186;
t4626 = t4357 * t4392;
t4636 = t4355 * t4392;
t4653 = t4353 * t4392;
t4740 = pkin(2) * t4392;
t4070 = -t4557 * t4626 - (-t4357 * t4552 + (-t4238 * t4362 * t4636 + t4357 * (t4364 * t4740 + t4346)) * t4693) * t4198 * t4537 + (-t4174 * t4278 - t4175 * t4274) * t4210 * t4653;
t4052 = t4066 * t4649 + t4070 * t4357;
t4130 = t4131 ^ 2;
t4393 = 0.1e1 / pkin(3) ^ 2;
t4695 = t4143 ^ 2 / t4213 ^ 2;
t4539 = t4393 * t4695;
t4115 = t4130 + t4539;
t4417 = -0.2e1 * t4131 * t4458;
t4763 = t4052 * t4364 + (-t4115 * t4623 + t4362 * t4417) * t4355;
t4699 = t4138 * t4381;
t4545 = t4366 * t4699;
t4185 = t4236 * t4354 + t4237 * t4356;
t4422 = t4236 * t4356 - t4237 * t4354;
t4617 = t4367 * t4381;
t4674 = (t4320 * t4373 + t4617) * t4357;
t4147 = (t4185 * t4336 + t4332 * t4422) * t4674 - t4263 * (t4332 * t4185 - t4422 * t4336);
t4688 = t4147 * t4214;
t4116 = t4545 - t4688;
t4168 = -t4243 * t4642 - t4402 * t4366;
t4171 = -t4247 * t4642 - t4403 * t4366;
t4526 = t4392 * t4688;
t4448 = t4373 * t4526;
t4618 = t4367 * t4372;
t4492 = t4355 * t4618;
t4681 = t4214 * t4366;
t4527 = t4147 * t4681;
t4641 = t4355 * t4373;
t4743 = pkin(2) * t4373;
t4067 = (-((t4138 * t4641 + t4357 * t4526) * t4738 + ((-t4527 + t4699) * t4367 + t4138 * t4743) * t4642 + t4357 * t4116) * t4138 - (t4355 * t4448 + (t4350 * t4357 - t4366 * t4492 - t4357) * t4138) * t4688 + (-t4168 * t4279 - t4171 * t4275) * t4353) * t4187;
t4176 = -t4243 * t4674 + t4247 * t4263;
t4179 = -t4243 * t4263 - t4247 * t4674;
t4239 = pkin(3) * t4618 + t4296;
t4199 = 0.1e1 / (t4239 * t4642 + t4283);
t4556 = (-t4381 * t4527 + (t4350 * t4391 + t4372 * t4735 + t4326) * t4138) * t4138 * t4187;
t4071 = -t4556 * t4626 - (-t4357 * t4545 + (-t4239 * t4366 * t4636 + t4357 * (t4372 * t4740 + t4350)) * t4688) * t4199 * t4526 + (-t4176 * t4279 - t4179 * t4275) * t4214 * t4653;
t4059 = t4067 * t4641 + t4071 * t4357;
t4134 = t4138 ^ 2;
t4692 = t4147 ^ 2 / t4223 ^ 2;
t4534 = t4393 * t4692;
t4121 = t4134 + t4534;
t4414 = -0.2e1 * t4138 * t4448;
t4762 = t4059 * t4372 + (-t4121 * t4618 + t4366 * t4414) * t4355;
t4701 = t4136 * t4381;
t4547 = t4368 * t4701;
t4183 = t4231 * t4354 + t4234 * t4356;
t4424 = t4231 * t4356 - t4234 * t4354;
t4613 = t4369 * t4381;
t4673 = (t4321 * t4375 + t4613) * t4357;
t4148 = (t4183 * t4337 + t4333 * t4424) * t4673 - t4264 * (t4333 * t4183 - t4424 * t4337);
t4686 = t4148 * t4216;
t4117 = t4547 - t4686;
t4169 = -t4244 * t4640 - t4400 * t4368;
t4172 = -t4248 * t4640 - t4401 * t4368;
t4524 = t4392 * t4686;
t4447 = t4375 * t4524;
t4614 = t4369 * t4374;
t4491 = t4355 * t4614;
t4679 = t4216 * t4368;
t4525 = t4148 * t4679;
t4639 = t4355 * t4375;
t4742 = pkin(2) * t4375;
t4068 = (-((t4136 * t4639 + t4357 * t4524) * t4737 + ((-t4525 + t4701) * t4369 + t4136 * t4742) * t4640 + t4357 * t4117) * t4136 - (t4355 * t4447 + (t4351 * t4357 - t4368 * t4491 - t4357) * t4136) * t4686 + (-t4169 * t4280 - t4172 * t4276) * t4353) * t4188;
t4177 = -t4244 * t4673 + t4248 * t4264;
t4180 = -t4244 * t4264 - t4248 * t4673;
t4240 = pkin(3) * t4614 + t4297;
t4200 = 0.1e1 / (t4240 * t4640 + t4284);
t4555 = (-t4381 * t4525 + (t4351 * t4391 + t4374 * t4735 + t4326) * t4136) * t4136 * t4188;
t4072 = -t4555 * t4626 - (-t4357 * t4547 + (-t4240 * t4368 * t4636 + t4357 * (t4374 * t4740 + t4351)) * t4686) * t4200 * t4524 + (-t4177 * t4280 - t4180 * t4276) * t4216 * t4653;
t4060 = t4068 * t4639 + t4072 * t4357;
t4132 = t4136 ^ 2;
t4691 = t4148 ^ 2 / t4224 ^ 2;
t4531 = t4393 * t4691;
t4119 = t4132 + t4531;
t4416 = -0.2e1 * t4136 * t4447;
t4761 = t4060 * t4374 + (-t4119 * t4614 + t4368 * t4416) * t4355;
t4700 = t4137 * t4381;
t4546 = t4370 * t4700;
t4184 = t4232 * t4354 + t4235 * t4356;
t4423 = t4232 * t4356 - t4235 * t4354;
t4609 = t4371 * t4381;
t4672 = (t4322 * t4377 + t4609) * t4357;
t4149 = (t4184 * t4338 + t4334 * t4423) * t4672 - t4265 * (t4334 * t4184 - t4423 * t4338);
t4684 = t4149 * t4218;
t4118 = t4546 - t4684;
t4170 = -t4245 * t4638 - t4398 * t4370;
t4173 = -t4249 * t4638 - t4399 * t4370;
t4522 = t4392 * t4684;
t4446 = t4377 * t4522;
t4610 = t4371 * t4376;
t4490 = t4355 * t4610;
t4677 = t4218 * t4370;
t4523 = t4149 * t4677;
t4637 = t4355 * t4377;
t4741 = pkin(2) * t4377;
t4069 = (-((t4137 * t4637 + t4357 * t4522) * t4736 + ((-t4523 + t4700) * t4371 + t4137 * t4741) * t4638 + t4357 * t4118) * t4137 - (t4355 * t4446 + (t4352 * t4357 - t4370 * t4490 - t4357) * t4137) * t4684 + (-t4170 * t4281 - t4173 * t4277) * t4353) * t4189;
t4178 = -t4245 * t4672 + t4249 * t4265;
t4181 = -t4245 * t4265 - t4249 * t4672;
t4241 = pkin(3) * t4610 + t4298;
t4201 = 0.1e1 / (t4241 * t4638 + t4285);
t4554 = (-t4381 * t4523 + (t4352 * t4391 + t4376 * t4735 + t4326) * t4137) * t4137 * t4189;
t4073 = -t4554 * t4626 - (-t4357 * t4546 + (-t4241 * t4370 * t4636 + t4357 * (t4376 * t4740 + t4352)) * t4684) * t4201 * t4522 + (-t4178 * t4281 - t4181 * t4277) * t4218 * t4653;
t4061 = t4069 * t4637 + t4073 * t4357;
t4133 = t4137 ^ 2;
t4690 = t4149 ^ 2 / t4225 ^ 2;
t4528 = t4393 * t4690;
t4120 = t4133 + t4528;
t4415 = -0.2e1 * t4137 * t4446;
t4760 = t4061 * t4376 + (-t4120 * t4610 + t4370 * t4415) * t4355;
t4625 = t4362 * t4363;
t4759 = -t4052 * t4362 + (t4115 * t4625 + t4364 * t4417) * t4355;
t4616 = t4368 * t4369;
t4758 = -t4060 * t4368 + (t4119 * t4616 + t4374 * t4416) * t4355;
t4612 = t4370 * t4371;
t4757 = -t4061 * t4370 + (t4120 * t4612 + t4376 * t4415) * t4355;
t4620 = t4366 * t4367;
t4756 = -t4059 * t4366 + (t4121 * t4620 + t4372 * t4414) * t4355;
t4301 = t4609 + t4741;
t4644 = t4355 * t4370;
t4406 = pkin(3) * t4644 - t4298 * t4357;
t4755 = t4301 * t4356 + t4406 * t4354;
t4300 = t4613 + t4742;
t4646 = t4355 * t4368;
t4407 = pkin(3) * t4646 - t4297 * t4357;
t4754 = t4300 * t4356 + t4407 * t4354;
t4299 = t4617 + t4743;
t4648 = t4355 * t4366;
t4408 = pkin(3) * t4648 - t4296 * t4357;
t4753 = t4299 * t4356 + t4408 * t4354;
t4295 = t4622 + t4744;
t4652 = t4355 * t4362;
t4409 = pkin(3) * t4652 - t4294 * t4357;
t4752 = t4295 * t4356 + t4409 * t4354;
t4750 = -2 * MDP(10);
t4749 = -2 * MDP(11);
t4327 = 0.2e1 * t4346 - 0.1e1;
t4328 = 0.2e1 * t4350 - 0.1e1;
t4329 = 0.2e1 * t4351 - 0.1e1;
t4330 = 0.2e1 * t4352 - 0.1e1;
t4748 = pkin(2) * t4131;
t4747 = pkin(2) * t4136;
t4746 = pkin(2) * t4137;
t4745 = pkin(2) * t4138;
t4734 = MDP(3) * t4355;
t4733 = MDP(4) * t4355;
t4732 = MDP(9) * t4392;
t4731 = t4066 * t4186;
t4730 = t4067 * t4187;
t4729 = t4068 * t4188;
t4728 = t4069 * t4189;
t4727 = t4070 * t4210;
t4726 = t4070 * t4362;
t4725 = t4070 * t4364;
t4724 = t4071 * t4214;
t4723 = t4071 * t4366;
t4722 = t4071 * t4372;
t4721 = t4072 * t4216;
t4720 = t4072 * t4368;
t4719 = t4072 * t4374;
t4718 = t4073 * t4218;
t4717 = t4073 * t4370;
t4716 = t4073 * t4376;
t4190 = t4295 * t4354 - t4409 * t4356;
t4250 = t4354 * t4634 - t4356 * t4365;
t4251 = t4354 * t4365 + t4356 * t4634;
t4585 = pkin(2) * t4652;
t4154 = -(t4250 * t4335 + t4251 * t4331) * t4739 + (-t4331 * t4190 + t4752 * t4335) * t4364 + t4246 * t4585;
t4155 = (-t4250 * t4331 + t4251 * t4335) * t4739 + (t4190 * t4335 + t4752 * t4331) * t4364 - t4242 * t4585;
t4397 = (-t4154 * t4278 - t4155 * t4274) * t4186;
t4603 = t4364 * t4557 + (pkin(2) * t4537 - t4108 * t4364) * t4186 * t4693;
t4074 = t4353 * t4397 + t4603;
t4715 = t4074 * t4186;
t4714 = t4074 * t4363;
t4713 = t4074 * t4365;
t4191 = t4299 * t4354 - t4408 * t4356;
t4252 = t4354 * t4632 - t4356 * t4373;
t4255 = t4354 * t4373 + t4356 * t4632;
t4584 = pkin(2) * t4648;
t4156 = -(t4252 * t4336 + t4255 * t4332) * t4738 + (-t4332 * t4191 + t4753 * t4336) * t4372 + t4247 * t4584;
t4159 = (-t4252 * t4332 + t4255 * t4336) * t4738 + (t4191 * t4336 + t4753 * t4332) * t4372 - t4243 * t4584;
t4396 = (-t4156 * t4279 - t4159 * t4275) * t4187;
t4596 = t4372 * t4556 + (pkin(2) * t4526 - t4116 * t4372) * t4187 * t4688;
t4075 = t4353 * t4396 + t4596;
t4712 = t4075 * t4187;
t4711 = t4075 * t4367;
t4710 = t4075 * t4373;
t4192 = t4300 * t4354 - t4407 * t4356;
t4253 = t4354 * t4630 - t4356 * t4375;
t4256 = t4354 * t4375 + t4356 * t4630;
t4583 = pkin(2) * t4646;
t4157 = -(t4253 * t4337 + t4256 * t4333) * t4737 + (-t4333 * t4192 + t4754 * t4337) * t4374 + t4248 * t4583;
t4160 = (-t4253 * t4333 + t4256 * t4337) * t4737 + (t4192 * t4337 + t4754 * t4333) * t4374 - t4244 * t4583;
t4395 = (-t4157 * t4280 - t4160 * t4276) * t4188;
t4595 = t4374 * t4555 + (pkin(2) * t4524 - t4117 * t4374) * t4188 * t4686;
t4076 = t4353 * t4395 + t4595;
t4709 = t4076 * t4188;
t4708 = t4076 * t4369;
t4707 = t4076 * t4375;
t4193 = t4301 * t4354 - t4406 * t4356;
t4254 = t4354 * t4628 - t4356 * t4377;
t4257 = t4354 * t4377 + t4356 * t4628;
t4582 = pkin(2) * t4644;
t4158 = -(t4254 * t4338 + t4257 * t4334) * t4736 + (-t4334 * t4193 + t4755 * t4338) * t4376 + t4249 * t4582;
t4161 = (-t4254 * t4334 + t4257 * t4338) * t4736 + (t4193 * t4338 + t4755 * t4334) * t4376 - t4245 * t4582;
t4394 = (-t4158 * t4281 - t4161 * t4277) * t4189;
t4594 = t4376 * t4554 + (pkin(2) * t4522 - t4118 * t4376) * t4189 * t4684;
t4077 = t4353 * t4394 + t4594;
t4706 = t4077 * t4189;
t4705 = t4077 * t4371;
t4704 = t4077 * t4377;
t4286 = t4354 * t4383 + t4356 * t4387;
t4287 = -t4354 * t4387 + t4356 * t4383;
t4226 = t4286 * t4343 + t4287 * t4344;
t4421 = t4286 * t4344 - t4287 * t4343;
t4162 = t4226 * t4335 - t4421 * t4331;
t4266 = t4365 * t4383 - t4387 * t4634;
t4267 = t4365 * t4387 + t4383 * t4634;
t4202 = t4266 * t4356 - t4267 * t4354;
t4203 = t4266 * t4354 + t4267 * t4356;
t4702 = (((-t4202 * t4343 + t4203 * t4344) * t4335 + (t4202 * t4344 + t4203 * t4343) * t4331) * t4362 + t4162 * t4650) * t4198;
t4288 = t4354 * t4384 + t4356 * t4388;
t4289 = -t4354 * t4388 + t4356 * t4384;
t4227 = t4288 * t4343 + t4289 * t4344;
t4420 = t4288 * t4344 - t4289 * t4343;
t4163 = t4227 * t4336 - t4420 * t4332;
t4268 = t4373 * t4384 - t4388 * t4632;
t4271 = t4373 * t4388 + t4384 * t4632;
t4204 = t4268 * t4356 - t4271 * t4354;
t4207 = t4268 * t4354 + t4271 * t4356;
t4698 = (((-t4204 * t4343 + t4207 * t4344) * t4336 + (t4204 * t4344 + t4207 * t4343) * t4332) * t4366 + t4163 * t4642) * t4199;
t4290 = t4354 * t4385 + t4356 * t4389;
t4291 = -t4354 * t4389 + t4356 * t4385;
t4228 = t4290 * t4343 + t4291 * t4344;
t4419 = t4290 * t4344 - t4291 * t4343;
t4164 = t4228 * t4337 - t4419 * t4333;
t4269 = t4375 * t4385 - t4389 * t4630;
t4272 = t4375 * t4389 + t4385 * t4630;
t4205 = t4269 * t4356 - t4272 * t4354;
t4208 = t4269 * t4354 + t4272 * t4356;
t4697 = (((-t4205 * t4343 + t4208 * t4344) * t4337 + t4333 * (t4205 * t4344 + t4208 * t4343)) * t4368 + t4164 * t4640) * t4200;
t4292 = t4354 * t4386 + t4356 * t4390;
t4293 = -t4354 * t4390 + t4356 * t4386;
t4229 = t4292 * t4343 + t4293 * t4344;
t4418 = t4292 * t4344 - t4293 * t4343;
t4165 = t4229 * t4338 - t4334 * t4418;
t4270 = t4377 * t4386 - t4390 * t4628;
t4273 = t4377 * t4390 + t4386 * t4628;
t4206 = t4270 * t4356 - t4273 * t4354;
t4209 = t4270 * t4354 + t4273 * t4356;
t4696 = (((-t4206 * t4343 + t4209 * t4344) * t4338 + (t4206 * t4344 + t4209 * t4343) * t4334) * t4370 + t4165 * t4638) * t4201;
t4694 = t4143 * t4186;
t4689 = t4147 * t4187;
t4687 = t4148 * t4188;
t4685 = t4149 * t4189;
t4682 = t4210 * t4364;
t4680 = t4214 * t4372;
t4678 = t4216 * t4374;
t4676 = t4218 * t4376;
t4663 = t4343 * t4353;
t4658 = t4344 * t4353;
t4627 = t4357 * t4381;
t4624 = t4362 * t4364;
t4621 = t4364 * t4365;
t4619 = t4366 * t4372;
t4615 = t4368 * t4374;
t4611 = t4370 * t4376;
t4608 = t4372 * t4373;
t4607 = t4374 * t4375;
t4606 = t4376 * t4377;
t4459 = t4357 * t4539;
t4488 = t4070 * t4623;
t4605 = -t4355 * t4488 - t4364 * t4459 + t4759;
t4489 = t4070 * t4625;
t4604 = -t4355 * t4489 - t4362 * t4459 + t4763;
t4452 = t4357 * t4531;
t4483 = t4072 * t4614;
t4602 = -t4355 * t4483 - t4374 * t4452 + t4758;
t4449 = t4357 * t4528;
t4482 = t4073 * t4610;
t4601 = -t4355 * t4482 - t4376 * t4449 + t4757;
t4455 = t4357 * t4534;
t4484 = t4071 * t4618;
t4600 = -t4355 * t4484 - t4372 * t4455 + t4756;
t4487 = t4071 * t4620;
t4599 = -t4355 * t4487 - t4366 * t4455 + t4762;
t4486 = t4072 * t4616;
t4598 = -t4355 * t4486 - t4368 * t4452 + t4761;
t4485 = t4073 * t4612;
t4597 = -t4355 * t4485 - t4370 * t4449 + t4760;
t4589 = -0.2e1 * t4694;
t4588 = -0.2e1 * t4689;
t4587 = -0.2e1 * t4687;
t4586 = -0.2e1 * t4685;
t4581 = MDP(5) * t4624;
t4580 = MDP(5) * t4619;
t4579 = MDP(5) * t4615;
t4578 = MDP(5) * t4611;
t4577 = t4066 * t4702;
t4345 = t4362 ^ 2;
t4576 = t4345 * t4731;
t4575 = t4067 * t4698;
t4347 = t4366 ^ 2;
t4574 = t4347 * t4730;
t4573 = t4068 * t4697;
t4348 = t4368 ^ 2;
t4572 = t4348 * t4729;
t4571 = t4069 * t4696;
t4349 = t4370 ^ 2;
t4570 = t4349 * t4728;
t4569 = t4186 * t4726;
t4568 = t4186 * t4725;
t4567 = t4187 * t4723;
t4566 = t4187 * t4722;
t4565 = t4188 * t4720;
t4564 = t4188 * t4719;
t4563 = t4189 * t4717;
t4562 = t4189 * t4716;
t4561 = t4074 * t4702;
t4560 = t4075 * t4698;
t4559 = t4076 * t4697;
t4558 = t4077 * t4696;
t4150 = t4162 * t4675 - t4259 * (t4226 * t4331 + t4421 * t4335);
t4553 = t4130 * t4150 * t4210;
t4152 = t4164 * t4673 - (t4228 * t4333 + t4419 * t4337) * t4264;
t4551 = t4132 * t4152 * t4216;
t4153 = t4165 * t4672 - (t4334 * t4229 + t4418 * t4338) * t4265;
t4550 = t4133 * t4153 * t4218;
t4151 = t4163 * t4674 - t4263 * (t4227 * t4332 + t4420 * t4336);
t4549 = t4134 * t4151 * t4214;
t4548 = t4364 * t4702;
t4544 = t4372 * t4698;
t4543 = t4374 * t4697;
t4542 = t4376 * t4696;
t4541 = t4362 * t4695;
t4540 = t4364 * t4695;
t4536 = t4366 * t4692;
t4535 = t4372 * t4692;
t4533 = t4368 * t4691;
t4532 = t4374 * t4691;
t4530 = t4370 * t4690;
t4529 = t4376 * t4690;
t4521 = t4174 * t4683;
t4520 = t4174 * t4682;
t4519 = t4175 * t4683;
t4518 = t4175 * t4682;
t4517 = t4176 * t4681;
t4516 = t4176 * t4680;
t4515 = t4177 * t4679;
t4514 = t4177 * t4678;
t4513 = t4178 * t4677;
t4512 = t4178 * t4676;
t4511 = t4179 * t4681;
t4510 = t4179 * t4680;
t4509 = t4180 * t4679;
t4508 = t4180 * t4678;
t4507 = t4181 * t4677;
t4506 = t4181 * t4676;
t4504 = t4210 * t4624;
t4502 = t4214 * t4619;
t4500 = t4216 * t4615;
t4498 = t4218 * t4611;
t4497 = t4318 * t4652;
t4496 = t4320 * t4648;
t4495 = t4321 * t4646;
t4494 = t4322 * t4644;
t4481 = t4694 * t4775;
t4480 = t4687 * t4774;
t4479 = t4685 * t4773;
t4478 = t4689 * t4772;
t4477 = t4166 * t4589;
t4476 = t4167 * t4589;
t4475 = t4168 * t4588;
t4474 = t4171 * t4588;
t4473 = t4169 * t4587;
t4472 = t4172 * t4587;
t4471 = t4170 * t4586;
t4470 = t4173 * t4586;
t4469 = pkin(6) * t4537;
t4468 = pkin(6) * t4526;
t4467 = pkin(6) * t4524;
t4466 = pkin(6) * t4522;
t4465 = t4695 * t4702;
t4464 = t4692 * t4698;
t4463 = t4691 * t4697;
t4462 = t4690 * t4696;
t4461 = t4186 * t4541;
t4460 = t4186 * t4540;
t4457 = t4187 * t4536;
t4456 = t4187 * t4535;
t4454 = t4188 * t4533;
t4453 = t4188 * t4532;
t4451 = t4189 * t4530;
t4450 = t4189 * t4529;
t4445 = t4318 * t4357;
t4444 = t4320 * t4357;
t4443 = t4321 * t4357;
t4442 = t4322 * t4357;
t4441 = 0.2e1 * t4624 * t4731;
t4440 = 0.2e1 * t4619 * t4730;
t4439 = 0.2e1 * t4615 * t4729;
t4438 = 0.2e1 * t4611 * t4728;
t4437 = t4327 * t4481;
t4436 = t4329 * t4480;
t4435 = t4330 * t4479;
t4434 = t4328 * t4478;
t4433 = t4066 * t4365 - t4130 * t4363;
t4432 = -t4066 * t4363 - t4130 * t4365;
t4431 = t4067 * t4373 - t4134 * t4367;
t4430 = -t4067 * t4367 - t4134 * t4373;
t4429 = t4068 * t4375 - t4132 * t4369;
t4428 = -t4068 * t4369 - t4132 * t4375;
t4427 = t4069 * t4377 - t4133 * t4371;
t4426 = -t4069 * t4371 - t4133 * t4377;
t4413 = t4066 * t4751 + t4074 * t4649;
t4412 = t4067 * t4751 + t4075 * t4641;
t4411 = t4068 * t4751 + t4076 * t4639;
t4410 = t4069 * t4751 + t4077 * t4637;
t4129 = (-t4277 * ((t4322 * t4245 + t4249 * t4627) * t4606 - (-t4245 * t4381 + t4249 * t4442) * t4610 + t4249 * t4494) + t4281 * ((-t4245 * t4627 + t4322 * t4249) * t4606 + (t4245 * t4442 + t4249 * t4381) * t4610 - t4245 * t4494)) / (-t4325 * t4638 + (t4490 + t4629) * t4322);
t4128 = (-t4276 * ((t4321 * t4244 + t4248 * t4627) * t4607 - (-t4244 * t4381 + t4248 * t4443) * t4614 + t4248 * t4495) + t4280 * ((-t4244 * t4627 + t4321 * t4248) * t4607 + (t4244 * t4443 + t4248 * t4381) * t4614 - t4244 * t4495)) / (-t4324 * t4640 + (t4491 + t4631) * t4321);
t4127 = (-t4275 * ((t4320 * t4243 + t4247 * t4627) * t4608 - (-t4243 * t4381 + t4247 * t4444) * t4618 + t4247 * t4496) + t4279 * ((-t4243 * t4627 + t4320 * t4247) * t4608 + (t4243 * t4444 + t4247 * t4381) * t4618 - t4243 * t4496)) / (-t4323 * t4642 + (t4492 + t4633) * t4320);
t4126 = (-t4274 * ((t4318 * t4242 + t4246 * t4627) * t4621 - (-t4242 * t4381 + t4246 * t4445) * t4623 + t4246 * t4497) + t4278 * ((-t4242 * t4627 + t4318 * t4246) * t4621 + (t4242 * t4445 + t4246 * t4381) * t4623 - t4242 * t4497)) / (-t4319 * t4650 + (t4493 + t4635) * t4318);
t4125 = t4330 * t4133;
t4124 = t4329 * t4132;
t4123 = t4328 * t4134;
t4122 = t4327 * t4130;
t4114 = t4376 * t4746 - t4370 * t4466 / 0.2e1;
t4113 = t4374 * t4747 - t4368 * t4467 / 0.2e1;
t4112 = t4372 * t4745 - t4366 * t4468 / 0.2e1;
t4111 = t4370 * t4746 + t4376 * t4466 / 0.2e1;
t4110 = t4368 * t4747 + t4374 * t4467 / 0.2e1;
t4109 = t4366 * t4745 + t4372 * t4468 / 0.2e1;
t4107 = t4364 * t4748 - t4362 * t4469 / 0.2e1;
t4106 = t4362 * t4748 + t4364 * t4469 / 0.2e1;
t4065 = pkin(6) * t4069 + t4077 * t4643;
t4064 = pkin(6) * t4068 + t4076 * t4645;
t4063 = pkin(6) * t4067 + t4075 * t4647;
t4062 = pkin(6) * t4066 + t4074 * t4651;
t4049 = -t4065 * t4376 - t4077 * t4629;
t4048 = -t4064 * t4374 - t4076 * t4631;
t4047 = -t4063 * t4372 - t4075 * t4633;
t4046 = t4077 * t4357 * t4376 - t4065 * t4370;
t4045 = t4076 * t4357 * t4374 - t4064 * t4368;
t4044 = t4075 * t4357 * t4372 - t4063 * t4366;
t4043 = -pkin(6) * t4716 - t4410 * t4370;
t4042 = -pkin(6) * t4717 + t4410 * t4376;
t4041 = -pkin(6) * t4722 - t4412 * t4366;
t4040 = -pkin(6) * t4723 + t4412 * t4372;
t4039 = -pkin(6) * t4719 - t4411 * t4368;
t4038 = -pkin(6) * t4720 + t4411 * t4374;
t4037 = -t4062 * t4364 - t4074 * t4635;
t4036 = t4074 * t4357 * t4364 - t4062 * t4362;
t4029 = -pkin(6) * t4725 - t4413 * t4362;
t4028 = -pkin(6) * t4726 + t4413 * t4364;
t1 = [(t4154 * t4715 + t4156 * t4712 + t4157 * t4709 + t4158 * t4706) * MDP(1) + (t4166 * t4731 + t4168 * t4730 + t4169 * t4729 + t4170 * t4728) * MDP(2) + ((t4427 * t4158 + t4170 * t4704) * t4189 + (t4429 * t4157 + t4169 * t4707) * t4188 + (t4431 * t4156 + t4168 * t4710) * t4187 + (t4433 * t4154 + t4166 * t4713) * t4186) * t4734 + ((t4426 * t4158 - t4170 * t4705) * t4189 + (t4428 * t4157 - t4169 * t4708) * t4188 + (t4430 * t4156 - t4168 * t4711) * t4187 + (t4432 * t4154 - t4166 * t4714) * t4186) * t4733 + (t4166 * t4576 + t4168 * t4574 + t4169 * t4572 + t4170 * t4570 + ((-t4133 * t4178 + t4170 * t4479) * t4498 + (-t4132 * t4177 + t4169 * t4480) * t4500 + (-t4134 * t4176 + t4168 * t4478) * t4502 + (-t4130 * t4174 + t4166 * t4481) * t4504) * t4392) * MDP(5) + (t4166 * t4441 + t4168 * t4440 + t4169 * t4439 + t4170 * t4438 + ((-t4125 * t4178 + t4170 * t4435) * t4218 + (-t4124 * t4177 + t4169 * t4436) * t4216 + (-t4123 * t4176 + t4168 * t4434) * t4214 + (-t4122 * t4174 + t4166 * t4437) * t4210) * t4392) * MDP(6) + (t4166 * t4569 + t4168 * t4567 + t4169 * t4565 + t4170 * t4563 + (t4166 * t4460 + t4168 * t4456 + t4169 * t4453 + t4170 * t4450) * t4393 + (t4066 * t4521 + t4067 * t4517 + t4068 * t4515 + t4069 * t4513) * t4392) * MDP(7) + (t4166 * t4568 + t4168 * t4566 + t4169 * t4564 + t4170 * t4562 + (-t4166 * t4461 - t4168 * t4457 - t4169 * t4454 - t4170 * t4451) * t4393 + (t4066 * t4520 + t4067 * t4516 + t4068 * t4514 + t4069 * t4512) * t4392) * MDP(8) + (t4174 * t4727 + t4176 * t4724 + t4177 * t4721 + t4178 * t4718) * t4732 + ((t4170 * t4042 + t4597 * t4158) * t4189 + (t4169 * t4038 + t4598 * t4157) * t4188 + (t4168 * t4040 + t4599 * t4156) * t4187 + (t4166 * t4028 + t4604 * t4154) * t4186 + ((t4178 * t4046 + t4111 * t4471) * t4218 + (t4177 * t4045 + t4110 * t4473) * t4216 + (t4176 * t4044 + t4109 * t4475) * t4214 + (t4174 * t4036 + t4106 * t4477) * t4210 + (t4130 * t4521 + t4132 * t4515 + t4133 * t4513 + t4134 * t4517) * pkin(2)) * t4392) * MDP(10) + ((t4170 * t4043 + t4601 * t4158) * t4189 + (t4169 * t4039 + t4602 * t4157) * t4188 + (t4168 * t4041 + t4600 * t4156) * t4187 + (t4166 * t4029 + t4605 * t4154) * t4186 + ((t4178 * t4049 + t4114 * t4471) * t4218 + (t4177 * t4048 + t4113 * t4473) * t4216 + (t4176 * t4047 + t4112 * t4475) * t4214 + (t4174 * t4037 + t4107 * t4477) * t4210 + (t4130 * t4520 + t4132 * t4514 + t4133 * t4512 + t4134 * t4516) * pkin(2)) * t4392) * MDP(11) - MDP(13) * t4658 + MDP(14) * t4663; (t4155 * t4715 + t4159 * t4712 + t4160 * t4709 + t4161 * t4706) * MDP(1) + (t4167 * t4731 + t4171 * t4730 + t4172 * t4729 + t4173 * t4728) * MDP(2) + ((t4427 * t4161 + t4173 * t4704) * t4189 + (t4429 * t4160 + t4172 * t4707) * t4188 + (t4431 * t4159 + t4171 * t4710) * t4187 + (t4433 * t4155 + t4167 * t4713) * t4186) * t4734 + ((t4426 * t4161 - t4173 * t4705) * t4189 + (t4428 * t4160 - t4172 * t4708) * t4188 + (t4430 * t4159 - t4171 * t4711) * t4187 + (t4432 * t4155 - t4167 * t4714) * t4186) * t4733 + (t4167 * t4576 + t4171 * t4574 + t4172 * t4572 + t4173 * t4570 + ((-t4133 * t4181 + t4173 * t4479) * t4498 + (-t4132 * t4180 + t4172 * t4480) * t4500 + (-t4134 * t4179 + t4171 * t4478) * t4502 + (-t4130 * t4175 + t4167 * t4481) * t4504) * t4392) * MDP(5) + (t4167 * t4441 + t4171 * t4440 + t4172 * t4439 + t4173 * t4438 + ((-t4181 * t4125 + t4173 * t4435) * t4218 + (-t4180 * t4124 + t4172 * t4436) * t4216 + (-t4179 * t4123 + t4171 * t4434) * t4214 + (-t4175 * t4122 + t4167 * t4437) * t4210) * t4392) * MDP(6) + (t4167 * t4569 + t4171 * t4567 + t4172 * t4565 + t4173 * t4563 + (t4167 * t4460 + t4171 * t4456 + t4172 * t4453 + t4173 * t4450) * t4393 + (t4066 * t4519 + t4067 * t4511 + t4068 * t4509 + t4069 * t4507) * t4392) * MDP(7) + (t4167 * t4568 + t4171 * t4566 + t4172 * t4564 + t4173 * t4562 + (-t4167 * t4461 - t4171 * t4457 - t4172 * t4454 - t4173 * t4451) * t4393 + (t4066 * t4518 + t4067 * t4510 + t4068 * t4508 + t4069 * t4506) * t4392) * MDP(8) + (t4175 * t4727 + t4179 * t4724 + t4180 * t4721 + t4181 * t4718) * t4732 + ((t4042 * t4173 + t4597 * t4161) * t4189 + (t4038 * t4172 + t4598 * t4160) * t4188 + (t4040 * t4171 + t4599 * t4159) * t4187 + (t4028 * t4167 + t4604 * t4155) * t4186 + ((t4181 * t4046 + t4111 * t4470) * t4218 + (t4180 * t4045 + t4110 * t4472) * t4216 + (t4179 * t4044 + t4109 * t4474) * t4214 + (t4175 * t4036 + t4106 * t4476) * t4210 + (t4130 * t4519 + t4132 * t4509 + t4133 * t4507 + t4134 * t4511) * pkin(2)) * t4392) * MDP(10) + ((t4043 * t4173 + t4601 * t4161) * t4189 + (t4039 * t4172 + t4602 * t4160) * t4188 + (t4041 * t4171 + t4600 * t4159) * t4187 + (t4029 * t4167 + t4605 * t4155) * t4186 + ((t4181 * t4049 + t4114 * t4470) * t4218 + (t4180 * t4048 + t4113 * t4472) * t4216 + (t4179 * t4047 + t4112 * t4474) * t4214 + (t4175 * t4037 + t4107 * t4476) * t4210 + (t4130 * t4518 + t4132 * t4508 + t4133 * t4506 + t4134 * t4510) * pkin(2)) * t4392) * MDP(11) - MDP(13) * t4663 - MDP(14) * t4658; (t4760 + t4761 + t4762 + t4763) * MDP(10) + (t4756 + t4757 + t4758 + t4759) * MDP(11) + ((-t4530 - t4533 - t4536 - t4541) * MDP(10) + (-t4529 - t4532 - t4535 - t4540) * MDP(11)) * t4393 * t4357 + ((t4427 + t4429 + t4431 + t4433) * MDP(3) + (t4426 + t4428 + t4430 + t4432) * MDP(4) + (-t4485 - t4486 - t4487 - t4489) * MDP(10) + (-t4482 - t4483 - t4484 - t4488) * MDP(11)) * t4355 + (t4594 + t4595 + t4596 + t4603 + (t4394 + t4395 + t4396 + t4397) * t4353) * MDP(1); (t4571 + t4573 + t4575 + t4577) * MDP(2) + (t4345 * t4577 + t4347 * t4575 + t4348 * t4573 + t4349 * t4571) * MDP(5) + 0.2e1 * (t4066 * t4362 * t4548 + t4067 * t4366 * t4544 + t4068 * t4368 * t4543 + t4069 * t4370 * t4542) * MDP(6) + (t4696 * t4717 + t4697 * t4720 + t4698 * t4723 + t4702 * t4726) * MDP(7) + (t4070 * t4548 + t4071 * t4544 + t4072 * t4543 + t4073 * t4542) * MDP(8) + (t4028 * t4702 + t4038 * t4697 + t4040 * t4698 + t4042 * t4696) * MDP(10) + (t4029 * t4702 + t4039 * t4697 + t4041 * t4698 + t4043 * t4696) * MDP(11) + (t4077 * MDP(1) + MDP(10) * t4597 + MDP(11) * t4601) * t4129 + (t4076 * MDP(1) + MDP(10) * t4598 + MDP(11) * t4602) * t4128 + (t4075 * MDP(1) + MDP(10) * t4599 + MDP(11) * t4600) * t4127 + (t4074 * MDP(1) + MDP(10) * t4604 + MDP(11) * t4605) * t4126 + ((t4364 * t4465 + t4372 * t4464 + t4374 * t4463 + t4376 * t4462) * MDP(7) + (-t4362 * t4465 - t4366 * t4464 - t4368 * t4463 - t4370 * t4462) * MDP(8)) * t4393 + ((t4365 * t4561 + t4373 * t4560 + t4375 * t4559 + t4377 * t4558) * MDP(3) + (-t4363 * t4561 - t4367 * t4560 - t4369 * t4559 - t4371 * t4558) * MDP(4) + (MDP(3) * t4427 + MDP(4) * t4426) * t4129 + (MDP(3) * t4429 + MDP(4) * t4428) * t4128 + (MDP(3) * t4431 + MDP(4) * t4430) * t4127 + (MDP(3) * t4433 + MDP(4) * t4432) * t4126) * t4355 + (((t4111 * t4750 + t4114 * t4749 + (MDP(6) * t4330 + t4578) * t4773) * t4149 * t4696 + (-t4133 * t4578 + t4046 * MDP(10) + t4049 * MDP(11) - t4125 * MDP(6) + t4073 * MDP(9) + (MDP(7) * t4370 + MDP(8) * t4376) * t4069) * t4153) * t4218 + ((t4110 * t4750 + t4113 * t4749 + (MDP(6) * t4329 + t4579) * t4774) * t4148 * t4697 + (-t4132 * t4579 + t4045 * MDP(10) + t4048 * MDP(11) - t4124 * MDP(6) + t4072 * MDP(9) + (MDP(7) * t4368 + MDP(8) * t4374) * t4068) * t4152) * t4216 + ((t4109 * t4750 + t4112 * t4749 + (MDP(6) * t4328 + t4580) * t4772) * t4147 * t4698 + (-t4134 * t4580 + t4044 * MDP(10) + t4047 * MDP(11) - t4123 * MDP(6) + t4071 * MDP(9) + (MDP(7) * t4366 + MDP(8) * t4372) * t4067) * t4151) * t4214 + ((t4106 * t4750 + t4107 * t4749 + (MDP(6) * t4327 + t4581) * t4775) * t4143 * t4702 + (-t4130 * t4581 + t4036 * MDP(10) + t4037 * MDP(11) - t4122 * MDP(6) + t4070 * MDP(9) + (MDP(7) * t4362 + MDP(8) * t4364) * t4066) * t4150) * t4210 + ((t4362 * t4553 + t4366 * t4549 + t4368 * t4551 + t4370 * t4550) * MDP(10) + (t4364 * t4553 + t4372 * t4549 + t4374 * t4551 + t4376 * t4550) * MDP(11)) * pkin(2)) * t4392;];
taucX  = t1;
