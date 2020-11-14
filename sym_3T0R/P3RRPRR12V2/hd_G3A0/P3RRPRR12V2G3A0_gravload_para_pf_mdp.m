% Calculate minimal parameter regressor of Gravitation load for parallel robot
% P3RRPRR12V2G3A0
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
%   see P3RRPRR12V2G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:30
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRPRR12V2G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(6,1),zeros(15,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [15 1]), ...
  'P3RRPRR12V2G3A0_gravload_para_pf_mdp: MDP has to be [15x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:30:35
% EndTime: 2020-08-06 19:30:38
% DurationCPUTime: 2.82s
% Computational Cost: add. (1503->206), mult. (2184->379), div. (180->6), fcn. (1953->18), ass. (0->154)
t4472 = sin(qJ(2,1));
t4473 = sin(qJ(1,1));
t4467 = legFrame(1,2);
t4440 = sin(t4467);
t4443 = cos(t4467);
t4419 = t4443 * g(1) - t4440 * g(2);
t4479 = cos(qJ(1,1));
t4494 = g(3) * t4479 + t4419 * t4473;
t4524 = MDP(10) - MDP(13);
t4437 = t4472 * qJ(3,1);
t4478 = cos(qJ(2,1));
t4521 = t4478 * pkin(2) + t4437;
t4422 = pkin(1) + t4521;
t4565 = MDP(12) - MDP(3);
t4591 = -g(3) * t4473 + t4419 * t4479;
t4596 = -MDP(14) * (t4419 * (-t4479 * pkin(5) + t4422 * t4473) + g(3) * (t4473 * pkin(5) + t4422 * t4479)) - t4494 * MDP(2) + t4565 * t4591;
t4599 = (-t4472 * t4494 * t4524 - t4596) * t4473;
t4470 = sin(qJ(2,2));
t4471 = sin(qJ(1,2));
t4466 = legFrame(2,2);
t4439 = sin(t4466);
t4442 = cos(t4466);
t4418 = t4442 * g(1) - t4439 * g(2);
t4477 = cos(qJ(1,2));
t4495 = g(3) * t4477 + t4418 * t4471;
t4436 = t4470 * qJ(3,2);
t4476 = cos(qJ(2,2));
t4522 = t4476 * pkin(2) + t4436;
t4421 = pkin(1) + t4522;
t4592 = -g(3) * t4471 + t4418 * t4477;
t4595 = -MDP(14) * (t4418 * (-t4477 * pkin(5) + t4421 * t4471) + g(3) * (t4471 * pkin(5) + t4421 * t4477)) - t4495 * MDP(2) + t4565 * t4592;
t4598 = (-t4470 * t4495 * t4524 - t4595) * t4471;
t4468 = sin(qJ(2,3));
t4469 = sin(qJ(1,3));
t4465 = legFrame(3,2);
t4438 = sin(t4465);
t4441 = cos(t4465);
t4417 = t4441 * g(1) - t4438 * g(2);
t4475 = cos(qJ(1,3));
t4496 = g(3) * t4475 + t4417 * t4469;
t4435 = t4468 * qJ(3,3);
t4474 = cos(qJ(2,3));
t4523 = t4474 * pkin(2) + t4435;
t4420 = pkin(1) + t4523;
t4593 = -g(3) * t4469 + t4417 * t4475;
t4594 = -MDP(14) * (t4417 * (-t4475 * pkin(5) + t4420 * t4469) + g(3) * (t4469 * pkin(5) + t4420 * t4475)) - t4496 * MDP(2) + t4565 * t4593;
t4597 = (-t4468 * t4496 * t4524 - t4594) * t4469;
t4414 = t4438 * g(1) + t4441 * g(2);
t4372 = -t4414 * t4474 + t4468 * t4593;
t4482 = 0.1e1 / qJ(3,3);
t4560 = t4372 * t4482;
t4415 = t4439 * g(1) + t4442 * g(2);
t4376 = -t4415 * t4476 + t4470 * t4592;
t4483 = 0.1e1 / qJ(3,2);
t4559 = t4376 * t4483;
t4416 = t4440 * g(1) + t4443 * g(2);
t4380 = -t4416 * t4478 + t4472 * t4591;
t4484 = 0.1e1 / qJ(3,1);
t4558 = t4380 * t4484;
t4481 = pkin(2) + pkin(3);
t4581 = -0.2e1 * t4481;
t4571 = t4438 * qJ(3,3);
t4570 = t4439 * qJ(3,2);
t4569 = t4440 * qJ(3,1);
t4568 = t4441 * qJ(3,3);
t4567 = t4442 * qJ(3,2);
t4566 = t4443 * qJ(3,1);
t4564 = MDP(9) + MDP(11);
t4563 = MDP(14) * t4372;
t4562 = MDP(14) * t4376;
t4561 = MDP(14) * t4380;
t4527 = t4481 * t4474;
t4493 = t4435 + pkin(1) + t4527;
t4480 = pkin(5) - pkin(6);
t4533 = t4480 * t4475;
t4557 = (t4493 * t4469 - t4533) * t4482;
t4526 = t4481 * t4476;
t4492 = t4436 + pkin(1) + t4526;
t4532 = t4480 * t4477;
t4556 = (t4492 * t4471 - t4532) * t4483;
t4525 = t4481 * t4478;
t4491 = t4437 + pkin(1) + t4525;
t4531 = t4480 * t4479;
t4555 = (t4491 * t4473 - t4531) * t4484;
t4542 = t4468 * t4475;
t4512 = qJ(3,3) * t4542;
t4520 = pkin(1) * t4475 + t4469 * t4480;
t4551 = (0.2e1 * t4512 + t4520) * t4481;
t4540 = t4470 * t4477;
t4513 = qJ(3,2) * t4540;
t4519 = pkin(1) * t4477 + t4471 * t4480;
t4550 = (0.2e1 * t4513 + t4519) * t4481;
t4538 = t4472 * t4479;
t4514 = qJ(3,1) * t4538;
t4518 = pkin(1) * t4479 + t4473 * t4480;
t4549 = (0.2e1 * t4514 + t4518) * t4481;
t4545 = (qJ(3,3) + t4481) * (-qJ(3,3) + t4481);
t4544 = (qJ(3,2) + t4481) * (-qJ(3,2) + t4481);
t4543 = (qJ(3,1) + t4481) * (-qJ(3,1) + t4481);
t4541 = t4468 * t4481;
t4539 = t4470 * t4481;
t4537 = t4472 * t4481;
t4536 = t4475 * t4481;
t4535 = t4477 * t4481;
t4534 = t4479 * t4481;
t4429 = pkin(1) * t4468 + qJ(3,3);
t4530 = t4481 * t4429;
t4430 = pkin(1) * t4470 + qJ(3,2);
t4529 = t4481 * t4430;
t4431 = pkin(1) * t4472 + qJ(3,1);
t4528 = t4481 * t4431;
t4517 = qJ(3,1) * t4581;
t4516 = qJ(3,2) * t4581;
t4515 = qJ(3,3) * t4581;
t4511 = t4474 * t4557;
t4510 = t4476 * t4556;
t4509 = t4478 * t4555;
t4507 = t4496 * t4469 * t4474;
t4505 = t4495 * t4471 * t4476;
t4503 = t4494 * t4473 * t4478;
t4502 = t4475 * t4545;
t4501 = t4477 * t4544;
t4500 = t4479 * t4543;
t4373 = t4414 * t4468 + t4474 * t4593;
t4377 = t4415 * t4470 + t4476 * t4592;
t4381 = t4416 * t4472 + t4478 * t4591;
t4366 = -t4414 * t4523 + t4593 * (t4468 * pkin(2) - t4474 * qJ(3,3));
t4487 = t4366 * MDP(14) + t4524 * t4373;
t4367 = -t4415 * t4522 + t4592 * (t4470 * pkin(2) - t4476 * qJ(3,2));
t4486 = t4367 * MDP(14) + t4524 * t4377;
t4368 = -t4416 * t4521 + t4591 * (t4472 * pkin(2) - t4478 * qJ(3,1));
t4485 = t4368 * MDP(14) + t4524 * t4381;
t4464 = t4478 ^ 2;
t4463 = t4476 ^ 2;
t4462 = t4474 ^ 2;
t4413 = 0.1e1 / t4491;
t4412 = 0.1e1 / t4492;
t4411 = 0.1e1 / t4493;
t4410 = pkin(1) * qJ(3,1) - t4472 * t4543;
t4409 = pkin(1) * qJ(3,2) - t4470 * t4544;
t4408 = pkin(1) * qJ(3,3) - t4468 * t4545;
t4407 = t4514 + t4518;
t4405 = t4513 + t4519;
t4403 = t4512 + t4520;
t4401 = qJ(3,1) * t4479 + t4518 * t4472;
t4400 = qJ(3,2) * t4477 + t4519 * t4470;
t4399 = qJ(3,3) * t4475 + t4520 * t4468;
t4365 = (t4443 * t4534 - t4569) * t4464 + (t4407 * t4443 + t4440 * t4537) * t4478 + t4440 * t4431;
t4364 = (-t4440 * t4534 - t4566) * t4464 + (-t4407 * t4440 + t4443 * t4537) * t4478 + t4443 * t4431;
t4363 = (t4442 * t4535 - t4570) * t4463 + (t4405 * t4442 + t4439 * t4539) * t4476 + t4439 * t4430;
t4362 = (-t4439 * t4535 - t4567) * t4463 + (-t4405 * t4439 + t4442 * t4539) * t4476 + t4442 * t4430;
t4361 = (t4441 * t4536 - t4571) * t4462 + (t4403 * t4441 + t4438 * t4541) * t4474 + t4438 * t4429;
t4360 = (-t4438 * t4536 - t4568) * t4462 + (-t4403 * t4438 + t4441 * t4541) * t4474 + t4441 * t4429;
t1 = [-g(1) * MDP(15) + t4564 * ((t4365 * t4558 - t4443 * t4503) * t4413 + (t4363 * t4559 - t4442 * t4505) * t4412 + (t4361 * t4560 - t4441 * t4507) * t4411) + ((-((t4440 * t4517 + t4443 * t4500) * t4464 + (-t4410 * t4440 + t4443 * t4549) * t4478 + t4401 * t4566 + t4440 * t4528) * t4561 + t4485 * t4365) * t4484 - t4443 * t4599) * t4413 + ((-((t4439 * t4516 + t4442 * t4501) * t4463 + (-t4409 * t4439 + t4442 * t4550) * t4476 + t4400 * t4567 + t4439 * t4529) * t4562 + t4486 * t4363) * t4483 - t4442 * t4598) * t4412 + ((-((t4438 * t4515 + t4441 * t4502) * t4462 + (-t4408 * t4438 + t4441 * t4551) * t4474 + t4399 * t4568 + t4438 * t4530) * t4563 + t4487 * t4361) * t4482 - t4441 * t4597) * t4411; -g(2) * MDP(15) + t4564 * ((t4364 * t4558 + t4440 * t4503) * t4413 + (t4362 * t4559 + t4439 * t4505) * t4412 + (t4360 * t4560 + t4438 * t4507) * t4411) + ((-((-t4440 * t4500 + t4443 * t4517) * t4464 + (-t4410 * t4443 - t4440 * t4549) * t4478 - t4401 * t4569 + t4443 * t4528) * t4561 + t4485 * t4364) * t4484 + t4440 * t4599) * t4413 + ((-((-t4439 * t4501 + t4442 * t4516) * t4463 + (-t4409 * t4442 - t4439 * t4550) * t4476 - t4400 * t4570 + t4442 * t4529) * t4562 + t4486 * t4362) * t4483 + t4439 * t4598) * t4412 + ((-((-t4438 * t4502 + t4441 * t4515) * t4462 + (-t4408 * t4441 - t4438 * t4551) * t4474 - t4399 * t4571 + t4441 * t4530) * t4563 + t4487 * t4360) * t4482 + t4438 * t4597) * t4411; -g(3) * MDP(15) + t4564 * ((-t4380 * t4555 - t4479 * t4494) * t4478 * t4413 + (-t4376 * t4556 - t4477 * t4495) * t4476 * t4412 + (-t4372 * t4557 - t4475 * t4496) * t4474 * t4411) + ((-t4368 * t4509 - (-t4473 * t4464 * t4543 - ((0.2e1 * t4437 + pkin(1)) * t4473 - t4531) * t4525 - qJ(3,1) * (t4431 * t4473 - t4472 * t4531)) * t4558) * MDP(14) + t4596 * t4479 - t4524 * (t4381 * t4509 - t4494 * t4538)) * t4413 + ((-t4367 * t4510 - (-t4471 * t4463 * t4544 - ((0.2e1 * t4436 + pkin(1)) * t4471 - t4532) * t4526 - qJ(3,2) * (t4430 * t4471 - t4470 * t4532)) * t4559) * MDP(14) + t4595 * t4477 - t4524 * (t4377 * t4510 - t4495 * t4540)) * t4412 + ((-t4366 * t4511 - (-t4469 * t4462 * t4545 - ((0.2e1 * t4435 + pkin(1)) * t4469 - t4533) * t4527 - qJ(3,3) * (t4429 * t4469 - t4468 * t4533)) * t4560) * MDP(14) + t4594 * t4475 - t4524 * (t4373 * t4511 - t4496 * t4542)) * t4411;];
taugX  = t1;
