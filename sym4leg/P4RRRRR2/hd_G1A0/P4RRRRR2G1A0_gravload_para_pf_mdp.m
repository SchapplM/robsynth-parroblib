% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P4RRRRR2G1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [4x1]
%   Generalized platform coordinates
% qJ [3x4]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [4x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% MDP [17x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P4RRRRR2G1A0_convert_par2_MPV_fixb.m

% Output:
% taugX [4x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 17:26
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P4RRRRR2G1A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(3,1),zeros(4,3),zeros(4,3),zeros(2,1),zeros(17,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: qJ has to be [3x4] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: Koppelpunkt has to be [4x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [17 1]), ...
  'P4RRRRR2G1A0_gravload_para_pf_mdp: MDP has to be [17x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 17:26:53
% EndTime: 2020-08-07 17:26:56
% DurationCPUTime: 2.68s
% Computational Cost: add. (2057->285), mult. (1785->449), div. (484->18), fcn. (1804->90), ass. (0->222)
t4453 = legFrame(4,3);
t4427 = sin(t4453);
t4560 = cos(t4453);
t4357 = g(1) * t4427 - t4560 * g(2);
t4361 = t4560 * g(1) + t4427 * g(2);
t4433 = qJ(1,4) + qJ(2,4);
t4419 = sin(t4433);
t4420 = cos(t4433);
t4329 = t4357 * t4420 + t4361 * t4419;
t4457 = sin(qJ(3,4));
t4565 = t4329 * t4457 ^ 2;
t4454 = legFrame(3,3);
t4428 = sin(t4454);
t4559 = cos(t4454);
t4358 = g(1) * t4428 - t4559 * g(2);
t4362 = t4559 * g(1) + t4428 * g(2);
t4450 = qJ(1,3) + qJ(2,3);
t4421 = sin(t4450);
t4424 = cos(t4450);
t4331 = t4358 * t4424 + t4362 * t4421;
t4461 = sin(qJ(3,3));
t4564 = t4331 * t4461 ^ 2;
t4455 = legFrame(2,3);
t4429 = sin(t4455);
t4558 = cos(t4455);
t4359 = g(1) * t4429 - t4558 * g(2);
t4363 = t4558 * g(1) + t4429 * g(2);
t4451 = qJ(1,2) + qJ(2,2);
t4422 = sin(t4451);
t4425 = cos(t4451);
t4332 = t4359 * t4425 + t4363 * t4422;
t4463 = sin(qJ(3,2));
t4563 = t4332 * t4463 ^ 2;
t4456 = legFrame(1,3);
t4430 = sin(t4456);
t4557 = cos(t4456);
t4360 = g(1) * t4430 - t4557 * g(2);
t4364 = t4557 * g(1) + t4430 * g(2);
t4452 = qJ(1,1) + qJ(2,1);
t4423 = sin(t4452);
t4426 = cos(t4452);
t4333 = t4360 * t4426 + t4364 * t4423;
t4465 = sin(qJ(3,1));
t4562 = t4333 * t4465 ^ 2;
t4561 = 2 * pkin(1);
t4516 = qJ(1,4) + t4453;
t4556 = pkin(1) * sin(t4516);
t4550 = qJ(1,3) + t4454;
t4555 = pkin(1) * sin(t4550);
t4551 = qJ(1,2) + t4455;
t4554 = pkin(1) * sin(t4551);
t4552 = qJ(1,1) + t4456;
t4553 = pkin(1) * sin(t4552);
t4549 = t4329 * t4457;
t4459 = cos(qJ(3,4));
t4548 = t4329 * t4459;
t4547 = t4331 * t4461;
t4467 = cos(qJ(3,3));
t4546 = t4331 * t4467;
t4545 = t4332 * t4463;
t4469 = cos(qJ(3,2));
t4544 = t4332 * t4469;
t4543 = t4333 * t4465;
t4471 = cos(qJ(3,1));
t4542 = t4333 * t4471;
t4408 = qJ(2,4) + t4516;
t4393 = qJ(3,4) + t4408;
t4372 = sin(t4393);
t4394 = -qJ(3,4) + t4408;
t4373 = sin(t4394);
t4345 = -0.2e1 * t4556 + (-t4372 - t4373) * pkin(2);
t4367 = 0.1e1 / (sin(qJ(2,4) + qJ(3,4)) + sin(qJ(2,4) - qJ(3,4)));
t4541 = t4345 * t4367;
t4346 = cos(t4516) * t4561 + (cos(t4393) + cos(t4394)) * pkin(2);
t4540 = t4346 * t4367;
t4412 = qJ(2,3) + t4550;
t4401 = qJ(3,3) + t4412;
t4378 = sin(t4401);
t4402 = -qJ(3,3) + t4412;
t4379 = sin(t4402);
t4347 = -0.2e1 * t4555 + (-t4378 - t4379) * pkin(2);
t4368 = 0.1e1 / (sin(qJ(2,3) + qJ(3,3)) + sin(qJ(2,3) - qJ(3,3)));
t4539 = t4347 * t4368;
t4413 = qJ(2,2) + t4551;
t4403 = qJ(3,2) + t4413;
t4380 = sin(t4403);
t4404 = -qJ(3,2) + t4413;
t4381 = sin(t4404);
t4348 = -0.2e1 * t4554 + (-t4380 - t4381) * pkin(2);
t4369 = 0.1e1 / (sin(qJ(2,2) + qJ(3,2)) + sin(qJ(2,2) - qJ(3,2)));
t4538 = t4348 * t4369;
t4414 = qJ(2,1) + t4552;
t4405 = qJ(3,1) + t4414;
t4382 = sin(t4405);
t4406 = -qJ(3,1) + t4414;
t4383 = sin(t4406);
t4349 = -0.2e1 * t4553 + (-t4382 - t4383) * pkin(2);
t4370 = 0.1e1 / (sin(qJ(2,1) + qJ(3,1)) + sin(qJ(2,1) - qJ(3,1)));
t4537 = t4349 * t4370;
t4350 = cos(t4550) * t4561 + (cos(t4401) + cos(t4402)) * pkin(2);
t4536 = t4350 * t4368;
t4351 = cos(t4551) * t4561 + (cos(t4403) + cos(t4404)) * pkin(2);
t4535 = t4351 * t4369;
t4352 = cos(t4552) * t4561 + (cos(t4405) + cos(t4406)) * pkin(2);
t4534 = t4352 * t4370;
t4473 = xP(4);
t4377 = -t4473 + t4408;
t4474 = koppelP(4,2);
t4478 = koppelP(4,1);
t4353 = -t4474 * cos(t4377) + t4478 * sin(t4377);
t4435 = 0.1e1 / sin(qJ(2,4));
t4533 = t4353 * t4435;
t4385 = -t4473 + t4412;
t4475 = koppelP(3,2);
t4479 = koppelP(3,1);
t4354 = -t4475 * cos(t4385) + t4479 * sin(t4385);
t4439 = 0.1e1 / sin(qJ(2,3));
t4532 = t4354 * t4439;
t4386 = -t4473 + t4413;
t4476 = koppelP(2,2);
t4480 = koppelP(2,1);
t4355 = -t4476 * cos(t4386) + t4480 * sin(t4386);
t4441 = 0.1e1 / sin(qJ(2,2));
t4531 = t4355 * t4441;
t4387 = -t4473 + t4414;
t4477 = koppelP(1,2);
t4481 = koppelP(1,1);
t4356 = -t4477 * cos(t4387) + t4481 * sin(t4387);
t4443 = 0.1e1 / sin(qJ(2,1));
t4530 = t4356 * t4443;
t4391 = sin(t4408);
t4529 = t4391 * t4435;
t4392 = cos(t4408);
t4528 = t4392 * t4435;
t4395 = sin(t4412);
t4527 = t4395 * t4439;
t4396 = sin(t4413);
t4526 = t4396 * t4441;
t4397 = sin(t4414);
t4525 = t4397 * t4443;
t4398 = cos(t4412);
t4524 = t4398 * t4439;
t4399 = cos(t4413);
t4523 = t4399 * t4441;
t4400 = cos(t4414);
t4522 = t4400 * t4443;
t4521 = t4435 * t4457;
t4520 = t4439 * t4461;
t4519 = t4441 * t4463;
t4518 = t4443 * t4465;
t4482 = 0.1e1 / pkin(2);
t4483 = 1 / pkin(1);
t4517 = t4482 * t4483;
t4515 = t4367 * t4549;
t4514 = t4367 * t4548;
t4513 = t4435 * t4548;
t4512 = t4368 * t4547;
t4511 = t4368 * t4546;
t4510 = t4439 * t4546;
t4509 = t4369 * t4545;
t4508 = t4369 * t4544;
t4507 = t4441 * t4544;
t4506 = t4370 * t4543;
t4505 = t4370 * t4542;
t4504 = t4443 * t4542;
t4371 = cos(qJ(2,4)) * pkin(1) + t4459 * pkin(2);
t4503 = t4371 * t4435 / t4459 ^ 2;
t4374 = cos(qJ(2,3)) * pkin(1) + t4467 * pkin(2);
t4502 = t4374 * t4439 / t4467 ^ 2;
t4375 = cos(qJ(2,2)) * pkin(1) + t4469 * pkin(2);
t4501 = t4375 * t4441 / t4469 ^ 2;
t4376 = cos(qJ(2,1)) * pkin(1) + t4471 * pkin(2);
t4500 = t4376 * t4443 / t4471 ^ 2;
t4436 = 0.1e1 / t4459;
t4499 = t4436 * t4521;
t4444 = 0.1e1 / t4467;
t4498 = t4444 * t4520;
t4446 = 0.1e1 / t4469;
t4497 = t4446 * t4519;
t4448 = 0.1e1 / t4471;
t4496 = t4448 * t4518;
t4495 = t4329 * t4521;
t4494 = t4331 * t4520;
t4493 = t4332 * t4519;
t4492 = t4333 * t4518;
t4491 = t4436 * t4495;
t4490 = t4444 * t4494;
t4489 = t4446 * t4493;
t4488 = t4448 * t4492;
t4487 = t4457 * t4503;
t4486 = t4461 * t4502;
t4485 = t4463 * t4501;
t4484 = t4465 * t4500;
t4330 = -t4357 * t4419 + t4361 * t4420;
t4334 = -t4358 * t4421 + t4362 * t4424;
t4335 = -t4359 * t4422 + t4363 * t4425;
t4336 = -t4360 * t4423 + t4364 * t4426;
t4472 = cos(qJ(1,1));
t4470 = cos(qJ(1,2));
t4468 = cos(qJ(1,3));
t4466 = sin(qJ(1,1));
t4464 = sin(qJ(1,2));
t4462 = sin(qJ(1,3));
t4460 = cos(qJ(1,4));
t4458 = sin(qJ(1,4));
t4432 = cos(t4473);
t4431 = sin(t4473);
t4366 = g(1) * t4432 + g(2) * t4431;
t4365 = g(1) * t4431 - g(2) * t4432;
t4344 = -t4360 * t4466 + t4364 * t4472;
t4343 = -t4359 * t4464 + t4363 * t4470;
t4342 = -t4358 * t4462 + t4362 * t4468;
t4341 = t4360 * t4472 + t4364 * t4466;
t4340 = t4359 * t4470 + t4363 * t4464;
t4339 = t4358 * t4468 + t4362 * t4462;
t4338 = -t4357 * t4458 + t4361 * t4460;
t4337 = t4357 * t4460 + t4361 * t4458;
t4328 = ((t4431 * t4481 + t4432 * t4477) * t4352 - 0.2e1 * (-t4431 * t4477 + t4432 * t4481) * (t4553 + (t4383 / 0.2e1 + t4382 / 0.2e1) * pkin(2))) * t4370 * t4517;
t4327 = ((t4431 * t4480 + t4432 * t4476) * t4351 - 0.2e1 * (-t4431 * t4476 + t4432 * t4480) * (t4554 + (t4381 / 0.2e1 + t4380 / 0.2e1) * pkin(2))) * t4369 * t4517;
t4326 = ((t4431 * t4479 + t4432 * t4475) * t4350 - 0.2e1 * (-t4431 * t4475 + t4432 * t4479) * (t4555 + (t4379 / 0.2e1 + t4378 / 0.2e1) * pkin(2))) * t4368 * t4517;
t4325 = ((t4431 * t4478 + t4432 * t4474) * t4346 - 0.2e1 * (-t4431 * t4474 + t4432 * t4478) * (t4556 + (t4373 / 0.2e1 + t4372 / 0.2e1) * pkin(2))) * t4367 * t4517;
t1 = [(-t4365 * t4431 - t4366 * t4432) * MDP(17) + ((t4337 * t4528 + t4339 * t4524 + t4340 * t4523 + t4341 * t4522) * MDP(2) + (t4338 * t4528 + t4342 * t4524 + t4343 * t4523 + t4344 * t4522) * MDP(3) + (t4329 * t4528 + t4331 * t4524 + t4332 * t4523 + t4333 * t4522) * MDP(5) + (t4330 * t4528 + t4334 * t4524 + t4335 * t4523 + t4336 * t4522) * MDP(6) + (t4392 * t4513 + t4398 * t4510 + t4399 * t4507 + t4400 * t4504) * MDP(12) + (-t4392 * t4495 - t4398 * t4494 - t4399 * t4493 - t4400 * t4492) * MDP(13) + ((-t4329 * t4540 - t4331 * t4536 - t4332 * t4535 - t4333 * t4534) * MDP(5) + (-t4330 * t4540 - t4334 * t4536 - t4335 * t4535 - t4336 * t4534) * MDP(6) + (-t4346 * t4514 - t4350 * t4511 - t4351 * t4508 - t4352 * t4505) * MDP(12) + (t4346 * t4515 + t4350 * t4512 + t4351 * t4509 + t4352 * t4506) * MDP(13)) * t4482) * t4483; (t4365 * t4432 - t4366 * t4431) * MDP(17) + ((t4337 * t4529 + t4339 * t4527 + t4340 * t4526 + t4341 * t4525) * MDP(2) + (t4338 * t4529 + t4342 * t4527 + t4343 * t4526 + t4344 * t4525) * MDP(3) + (t4329 * t4529 + t4331 * t4527 + t4332 * t4526 + t4333 * t4525) * MDP(5) + (t4330 * t4529 + t4334 * t4527 + t4335 * t4526 + t4336 * t4525) * MDP(6) + (t4391 * t4513 + t4395 * t4510 + t4396 * t4507 + t4397 * t4504) * MDP(12) + (-t4391 * t4495 - t4395 * t4494 - t4396 * t4493 - t4397 * t4492) * MDP(13) + ((t4329 * t4541 + t4331 * t4539 + t4332 * t4538 + t4333 * t4537) * MDP(5) + (t4330 * t4541 + t4334 * t4539 + t4335 * t4538 + t4336 * t4537) * MDP(6) + (t4345 * t4514 + t4347 * t4511 + t4348 * t4508 + t4349 * t4505) * MDP(12) + (-t4345 * t4515 - t4347 * t4512 - t4348 * t4509 - t4349 * t4506) * MDP(13)) * t4482) * t4483; -g(3) * MDP(17) + ((t4448 * (-g(3) * t4471 + t4336 * t4465) + t4446 * (-g(3) * t4469 + t4335 * t4463) + t4444 * (-g(3) * t4467 + t4334 * t4461) + t4436 * (-g(3) * t4459 + t4330 * t4457)) * MDP(12) + (t4448 * (g(3) * t4465 + t4336 * t4471) + t4446 * (g(3) * t4463 + t4335 * t4469) + t4444 * (g(3) * t4461 + t4334 * t4467) + t4436 * (g(3) * t4457 + t4330 * t4459)) * MDP(13)) * t4482 + ((t4337 * t4499 + t4339 * t4498 + t4340 * t4497 + t4341 * t4496) * MDP(2) + (t4338 * t4499 + t4342 * t4498 + t4343 * t4497 + t4344 * t4496) * MDP(3) + (t4488 + t4489 + t4490 + t4491) * MDP(5) + (t4330 * t4499 + t4334 * t4498 + t4335 * t4497 + t4336 * t4496) * MDP(6) + (t4492 + t4493 + t4494 + t4495) * MDP(12) + (-t4436 * t4435 * t4565 - t4444 * t4439 * t4564 - t4446 * t4441 * t4563 - t4448 * t4443 * t4562) * MDP(13) + ((-t4329 * t4487 - t4331 * t4486 - t4332 * t4485 - t4333 * t4484) * MDP(5) + (-t4330 * t4487 - t4334 * t4486 - t4335 * t4485 - t4336 * t4484) * MDP(6) + (-t4371 * t4491 - t4374 * t4490 - t4375 * t4489 - t4376 * t4488) * MDP(12) + (t4500 * t4562 + t4501 * t4563 + t4502 * t4564 + t4503 * t4565) * MDP(13)) * t4482) * t4483; (t4325 * t4329 + t4326 * t4331 + t4327 * t4332 + t4328 * t4333) * MDP(5) + (t4325 * t4330 + t4326 * t4334 + t4327 * t4335 + t4328 * t4336) * MDP(6) + (t4325 * t4548 + t4326 * t4546 + t4327 * t4544 + t4328 * t4542) * MDP(12) + (-t4325 * t4549 - t4326 * t4547 - t4327 * t4545 - t4328 * t4543) * MDP(13) + t4365 * MDP(15) + t4366 * MDP(16) + ((t4337 * t4533 + t4339 * t4532 + t4340 * t4531 + t4341 * t4530) * MDP(2) + (t4338 * t4533 + t4342 * t4532 + t4343 * t4531 + t4344 * t4530) * MDP(3) + (t4329 * t4533 + t4331 * t4532 + t4332 * t4531 + t4333 * t4530) * MDP(5) + (t4330 * t4533 + t4334 * t4532 + t4335 * t4531 + t4336 * t4530) * MDP(6) + (t4353 * t4513 + t4354 * t4510 + t4355 * t4507 + t4356 * t4504) * MDP(12) + (-t4353 * t4495 - t4354 * t4494 - t4355 * t4493 - t4356 * t4492) * MDP(13)) * t4483;];
taugX  = t1;
