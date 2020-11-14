% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [15x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRPRR12V2G2A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:23
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G2A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V2G2A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:23:53
% EndTime: 2020-08-06 19:23:56
% DurationCPUTime: 2.97s
% Computational Cost: add. (1503->215), mult. (2184->379), div. (180->6), fcn. (1953->18), ass. (0->157)
t4424 = legFrame(3,2);
t4400 = sin(t4424);
t4403 = cos(t4424);
t4382 = g(1) * t4403 - g(2) * t4400;
t4428 = sin(qJ(1,3));
t4434 = cos(qJ(1,3));
t4358 = -g(3) * t4428 + t4382 * t4434;
t4427 = sin(qJ(2,3));
t4433 = cos(qJ(2,3));
t4406 = t4433 * pkin(2);
t4397 = t4427 * qJ(3,3);
t4479 = t4406 + t4397;
t4511 = MDP(12) - MDP(3);
t4520 = t4397 + pkin(1);
t4546 = g(3) * t4434 + t4382 * t4428;
t4551 = -t4511 * t4546 + MDP(14) * (g(3) * (-t4434 * pkin(5) + (pkin(1) + t4479) * t4428) - t4382 * (t4428 * pkin(5) + (t4406 + t4520) * t4434));
t4552 = MDP(10) - MDP(13);
t4555 = (-(-t4552 * t4427 + MDP(2)) * t4358 + t4551) * t4434;
t4425 = legFrame(2,2);
t4401 = sin(t4425);
t4404 = cos(t4425);
t4383 = g(1) * t4404 - g(2) * t4401;
t4430 = sin(qJ(1,2));
t4436 = cos(qJ(1,2));
t4359 = -g(3) * t4430 + t4383 * t4436;
t4429 = sin(qJ(2,2));
t4435 = cos(qJ(2,2));
t4407 = t4435 * pkin(2);
t4398 = t4429 * qJ(3,2);
t4478 = t4407 + t4398;
t4519 = t4398 + pkin(1);
t4547 = g(3) * t4436 + t4383 * t4430;
t4549 = -t4511 * t4547 + MDP(14) * ((-t4436 * pkin(5) + (pkin(1) + t4478) * t4430) * g(3) - t4383 * (t4430 * pkin(5) + (t4407 + t4519) * t4436));
t4554 = (-(-t4552 * t4429 + MDP(2)) * t4359 + t4549) * t4436;
t4426 = legFrame(1,2);
t4402 = sin(t4426);
t4405 = cos(t4426);
t4384 = g(1) * t4405 - g(2) * t4402;
t4432 = sin(qJ(1,1));
t4438 = cos(qJ(1,1));
t4360 = -g(3) * t4432 + t4384 * t4438;
t4431 = sin(qJ(2,1));
t4437 = cos(qJ(2,1));
t4408 = t4437 * pkin(2);
t4399 = t4431 * qJ(3,1);
t4477 = t4408 + t4399;
t4518 = t4399 + pkin(1);
t4548 = g(3) * t4438 + t4384 * t4432;
t4550 = -t4511 * t4548 + MDP(14) * ((-t4438 * pkin(5) + (pkin(1) + t4477) * t4432) * g(3) - (t4432 * pkin(5) + (t4408 + t4518) * t4438) * t4384);
t4553 = (-(-t4552 * t4431 + MDP(2)) * t4360 + t4550) * t4438;
t4379 = g(1) * t4400 + g(2) * t4403;
t4337 = -t4379 * t4433 + t4427 * t4546;
t4441 = 0.1e1 / qJ(3,3);
t4506 = t4337 * t4441;
t4380 = g(1) * t4401 + g(2) * t4404;
t4341 = -t4380 * t4435 + t4429 * t4547;
t4442 = 0.1e1 / qJ(3,2);
t4505 = t4341 * t4442;
t4381 = g(1) * t4402 + g(2) * t4405;
t4345 = -t4381 * t4437 + t4431 * t4548;
t4443 = 0.1e1 / qJ(3,1);
t4504 = t4345 * t4443;
t4440 = pkin(2) + pkin(3);
t4530 = -0.2e1 * t4440;
t4517 = qJ(3,1) * t4402;
t4516 = qJ(3,1) * t4405;
t4515 = qJ(3,2) * t4401;
t4514 = qJ(3,2) * t4404;
t4513 = qJ(3,3) * t4400;
t4512 = qJ(3,3) * t4403;
t4510 = MDP(9) + MDP(11);
t4509 = MDP(14) * t4337;
t4508 = MDP(14) * t4341;
t4507 = MDP(14) * t4345;
t4439 = pkin(5) - pkin(6);
t4394 = t4428 * t4439;
t4485 = t4433 * t4440;
t4452 = t4485 + t4520;
t4503 = (t4452 * t4434 + t4394) * t4441;
t4395 = t4430 * t4439;
t4484 = t4435 * t4440;
t4451 = t4484 + t4519;
t4502 = (t4451 * t4436 + t4395) * t4442;
t4396 = t4432 * t4439;
t4483 = t4437 * t4440;
t4450 = t4483 + t4518;
t4501 = (t4450 * t4438 + t4396) * t4443;
t4455 = pkin(1) * t4428 - t4434 * t4439;
t4494 = t4427 * t4428;
t4471 = qJ(3,3) * t4494;
t4500 = (t4455 + 0.2e1 * t4471) * t4440;
t4454 = pkin(1) * t4430 - t4436 * t4439;
t4491 = t4429 * t4430;
t4472 = qJ(3,2) * t4491;
t4499 = (t4454 + 0.2e1 * t4472) * t4440;
t4453 = pkin(1) * t4432 - t4438 * t4439;
t4488 = t4431 * t4432;
t4473 = qJ(3,1) * t4488;
t4498 = (t4453 + 0.2e1 * t4473) * t4440;
t4497 = (qJ(3,3) + t4440) * (-qJ(3,3) + t4440);
t4496 = (qJ(3,2) + t4440) * (-qJ(3,2) + t4440);
t4495 = (qJ(3,1) + t4440) * (-qJ(3,1) + t4440);
t4493 = t4427 * t4440;
t4492 = t4428 * t4440;
t4490 = t4429 * t4440;
t4489 = t4430 * t4440;
t4487 = t4431 * t4440;
t4486 = t4432 * t4440;
t4391 = pkin(1) * t4427 + qJ(3,3);
t4482 = t4440 * t4391;
t4392 = pkin(1) * t4429 + qJ(3,2);
t4481 = t4440 * t4392;
t4393 = pkin(1) * t4431 + qJ(3,1);
t4480 = t4440 * t4393;
t4476 = qJ(3,1) * t4530;
t4475 = qJ(3,2) * t4530;
t4474 = qJ(3,3) * t4530;
t4470 = t4433 * t4503;
t4469 = t4435 * t4502;
t4468 = t4437 * t4501;
t4467 = t4358 * t4433 * t4434;
t4466 = t4359 * t4435 * t4436;
t4465 = t4360 * t4437 * t4438;
t4464 = t4428 * t4497;
t4463 = t4430 * t4496;
t4462 = t4432 * t4495;
t4338 = t4379 * t4427 + t4433 * t4546;
t4342 = t4380 * t4429 + t4435 * t4547;
t4346 = t4381 * t4431 + t4437 * t4548;
t4331 = -t4379 * t4479 + t4546 * (pkin(2) * t4427 - qJ(3,3) * t4433);
t4446 = t4331 * MDP(14) + t4552 * t4338;
t4332 = -t4380 * t4478 + t4547 * (pkin(2) * t4429 - qJ(3,2) * t4435);
t4445 = t4332 * MDP(14) + t4552 * t4342;
t4333 = -t4381 * t4477 + t4548 * (pkin(2) * t4431 - qJ(3,1) * t4437);
t4444 = t4333 * MDP(14) + t4552 * t4346;
t4423 = t4437 ^ 2;
t4422 = t4435 ^ 2;
t4421 = t4433 ^ 2;
t4378 = 0.1e1 / t4450;
t4377 = 0.1e1 / t4451;
t4376 = 0.1e1 / t4452;
t4375 = pkin(1) * qJ(3,1) - t4431 * t4495;
t4374 = pkin(1) * qJ(3,2) - t4429 * t4496;
t4373 = pkin(1) * qJ(3,3) - t4427 * t4497;
t4372 = t4453 + t4473;
t4370 = t4454 + t4472;
t4368 = t4455 + t4471;
t4366 = t4432 * qJ(3,1) + t4453 * t4431;
t4365 = t4430 * qJ(3,2) + t4454 * t4429;
t4364 = t4428 * qJ(3,3) + t4455 * t4427;
t4330 = (t4405 * t4486 - t4517) * t4423 + (t4372 * t4405 + t4402 * t4487) * t4437 + t4402 * t4393;
t4329 = (-t4402 * t4486 - t4516) * t4423 + (-t4372 * t4402 + t4405 * t4487) * t4437 + t4405 * t4393;
t4328 = (t4404 * t4489 - t4515) * t4422 + (t4370 * t4404 + t4401 * t4490) * t4435 + t4401 * t4392;
t4327 = (-t4401 * t4489 - t4514) * t4422 + (-t4370 * t4401 + t4404 * t4490) * t4435 + t4404 * t4392;
t4326 = (t4403 * t4492 - t4513) * t4421 + (t4368 * t4403 + t4400 * t4493) * t4433 + t4400 * t4391;
t4325 = (-t4400 * t4492 - t4512) * t4421 + (-t4368 * t4400 + t4403 * t4493) * t4433 + t4403 * t4391;
t1 = [-g(1) * MDP(15) + t4510 * ((t4330 * t4504 - t4405 * t4465) * t4378 + (t4328 * t4505 - t4404 * t4466) * t4377 + (t4326 * t4506 - t4403 * t4467) * t4376) + ((-((t4402 * t4476 + t4405 * t4462) * t4423 + (-t4375 * t4402 + t4405 * t4498) * t4437 + t4366 * t4516 + t4402 * t4480) * t4507 + t4444 * t4330) * t4443 + t4405 * t4553) * t4378 + ((-((t4401 * t4475 + t4404 * t4463) * t4422 + (-t4374 * t4401 + t4404 * t4499) * t4435 + t4365 * t4514 + t4401 * t4481) * t4508 + t4445 * t4328) * t4442 + t4404 * t4554) * t4377 + ((-((t4400 * t4474 + t4403 * t4464) * t4421 + (-t4373 * t4400 + t4403 * t4500) * t4433 + t4364 * t4512 + t4400 * t4482) * t4509 + t4446 * t4326) * t4441 + t4403 * t4555) * t4376; -g(2) * MDP(15) + t4510 * ((t4329 * t4504 + t4402 * t4465) * t4378 + (t4327 * t4505 + t4401 * t4466) * t4377 + (t4325 * t4506 + t4400 * t4467) * t4376) + ((-((-t4402 * t4462 + t4405 * t4476) * t4423 + (-t4375 * t4405 - t4402 * t4498) * t4437 - t4366 * t4517 + t4405 * t4480) * t4507 + t4444 * t4329) * t4443 - t4402 * t4553) * t4378 + ((-((-t4401 * t4463 + t4404 * t4475) * t4422 + (-t4374 * t4404 - t4401 * t4499) * t4435 - t4365 * t4515 + t4404 * t4481) * t4508 + t4445 * t4327) * t4442 - t4401 * t4554) * t4377 + ((-((-t4400 * t4464 + t4403 * t4474) * t4421 + (-t4373 * t4403 - t4400 * t4500) * t4433 - t4364 * t4513 + t4403 * t4482) * t4509 + t4446 * t4325) * t4441 - t4400 * t4555) * t4376; -g(3) * MDP(15) + t4510 * ((t4345 * t4501 + t4360 * t4432) * t4437 * t4378 + (t4341 * t4502 + t4359 * t4430) * t4435 * t4377 + (t4337 * t4503 + t4358 * t4428) * t4433 * t4376) + ((t4333 * t4468 - (t4438 * t4423 * t4495 + ((0.2e1 * t4399 + pkin(1)) * t4438 + t4396) * t4483 + qJ(3,1) * (t4393 * t4438 + t4431 * t4396)) * t4504) * MDP(14) + (t4360 * MDP(2) - t4550) * t4432 + t4552 * (t4346 * t4468 - t4360 * t4488)) * t4378 + ((t4332 * t4469 - (t4436 * t4422 * t4496 + ((0.2e1 * t4398 + pkin(1)) * t4436 + t4395) * t4484 + qJ(3,2) * (t4392 * t4436 + t4429 * t4395)) * t4505) * MDP(14) + (t4359 * MDP(2) - t4549) * t4430 + t4552 * (t4342 * t4469 - t4359 * t4491)) * t4377 + ((t4331 * t4470 - (t4434 * t4421 * t4497 + ((0.2e1 * t4397 + pkin(1)) * t4434 + t4394) * t4485 + qJ(3,3) * (t4391 * t4434 + t4427 * t4394)) * t4506) * MDP(14) + (t4358 * MDP(2) - t4551) * t4428 + t4552 * (t4338 * t4470 - t4358 * t4494)) * t4376;];
taugX  = t1;
