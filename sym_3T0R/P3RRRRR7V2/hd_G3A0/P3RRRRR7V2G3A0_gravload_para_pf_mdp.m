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
% MDP [18x1]
%   Minimal dynamic parameter vector for parallel robot(fixed base model)
%   see P3RRRRR7V2G3A0_convert_par2_MPV_fixb.m

% Output:
% taugX [3x1]
%   minimal parameter regressor of Gravitation load for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 10:47
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taugX = P3RRRRR7V2G3A0_gravload_para_pf_mdp(xP, qJ, g, legFrame, ...
  koppelP, pkin, MDP)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,1),zeros(3,3),zeros(3,3),zeros(7,1),zeros(18,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: qJ has to be [3x3] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: Koppelpunkt has to be [3x3] (double)');
assert(isreal(MDP) && all(size(MDP) == [18 1]), ...
  'P3RRRRR7V2G3A0_gravload_para_pf_mdp: MDP has to be [18x1] (double)'); 

%% Symbolic Calculation
% From invdyn_para_plfcoord_taugreg_mdp_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:47:45
% EndTime: 2020-08-07 10:47:48
% DurationCPUTime: 3.70s
% Computational Cost: add. (2196->360), mult. (3723->541), div. (210->14), fcn. (2808->66), ass. (0->265)
t8504 = legFrame(1,2);
t8470 = sin(t8504);
t8473 = cos(t8504);
t8437 = t8473 * g(1) - t8470 * g(2);
t8513 = sin(qJ(1,1));
t8522 = cos(qJ(1,1));
t8554 = g(3) * t8522 + t8437 * t8513;
t8499 = qJ(2,1) + qJ(3,1);
t8461 = cos(t8499);
t8521 = cos(qJ(2,1));
t8431 = 0.1e1 / (t8521 * pkin(2) + pkin(3) * t8461 + pkin(1));
t8633 = t8431 * t8522;
t8572 = t8554 * t8633;
t8503 = legFrame(2,2);
t8469 = sin(t8503);
t8472 = cos(t8503);
t8436 = t8472 * g(1) - t8469 * g(2);
t8510 = sin(qJ(1,2));
t8519 = cos(qJ(1,2));
t8555 = g(3) * t8519 + t8436 * t8510;
t8498 = qJ(2,2) + qJ(3,2);
t8460 = cos(t8498);
t8518 = cos(qJ(2,2));
t8430 = 0.1e1 / (t8518 * pkin(2) + pkin(3) * t8460 + pkin(1));
t8635 = t8430 * t8519;
t8574 = t8555 * t8635;
t8502 = legFrame(3,2);
t8468 = sin(t8502);
t8471 = cos(t8502);
t8435 = t8471 * g(1) - t8468 * g(2);
t8507 = sin(qJ(1,3));
t8516 = cos(qJ(1,3));
t8556 = g(3) * t8516 + t8435 * t8507;
t8497 = qJ(2,3) + qJ(3,3);
t8459 = cos(t8497);
t8515 = cos(qJ(2,3));
t8429 = 0.1e1 / (t8515 * pkin(2) + pkin(3) * t8459 + pkin(1));
t8637 = t8429 * t8516;
t8576 = t8556 * t8637;
t8686 = -g(3) * t8513 + t8437 * t8522;
t8685 = -g(3) * t8510 + t8436 * t8519;
t8684 = -g(3) * t8507 + t8435 * t8516;
t8683 = 0.2e1 * pkin(3);
t8487 = (pkin(5) + pkin(6) + pkin(7));
t8682 = 2 * t8487;
t8681 = 0.2e1 * t8515 ^ 2;
t8680 = 0.2e1 * t8518 ^ 2;
t8679 = 0.2e1 * t8521 ^ 2;
t8514 = cos(qJ(3,3));
t8491 = t8514 ^ 2;
t8480 = pkin(3) * t8491;
t8517 = cos(qJ(3,2));
t8493 = t8517 ^ 2;
t8481 = pkin(3) * t8493;
t8520 = cos(qJ(3,1));
t8495 = t8520 ^ 2;
t8482 = pkin(3) * t8495;
t8505 = sin(qJ(3,3));
t8675 = t8505 * pkin(1);
t8674 = t8505 * pkin(3);
t8508 = sin(qJ(3,2));
t8673 = t8508 * pkin(1);
t8672 = t8508 * pkin(3);
t8511 = sin(qJ(3,1));
t8671 = t8511 * pkin(1);
t8670 = t8511 * pkin(3);
t8669 = t8514 * pkin(2);
t8474 = t8514 * pkin(3);
t8668 = t8517 * pkin(2);
t8475 = t8517 * pkin(3);
t8667 = t8520 * pkin(2);
t8476 = t8520 * pkin(3);
t8666 = -qJ(3,1) + qJ(1,1);
t8665 = qJ(3,1) + qJ(1,1);
t8664 = -qJ(3,2) + qJ(1,2);
t8663 = qJ(3,2) + qJ(1,2);
t8662 = -qJ(3,3) + qJ(1,3);
t8661 = qJ(3,3) + qJ(1,3);
t8660 = qJ(1,1) + 0.2e1 * qJ(2,1);
t8659 = qJ(1,1) - 0.2e1 * qJ(2,1);
t8658 = qJ(1,2) + 0.2e1 * qJ(2,2);
t8657 = qJ(1,2) - 0.2e1 * qJ(2,2);
t8656 = qJ(1,3) + 0.2e1 * qJ(2,3);
t8655 = -0.2e1 * qJ(2,3) + qJ(1,3);
t8654 = 0.2e1 * pkin(1);
t8456 = sin(t8497);
t8462 = qJ(2,3) + t8661;
t8463 = -qJ(2,3) + t8662;
t8506 = sin(qJ(2,3));
t8526 = 0.2e1 * qJ(3,3);
t8537 = pkin(2) ^ 2;
t8653 = ((-sin(t8463) - sin(t8462)) * t8654 + (cos(t8463) + cos(t8462)) * t8682 + (sin(0.2e1 * qJ(3,3) - t8655) - sin(t8526 + t8656) - 0.2e1 * t8507) * pkin(3) + (sin(qJ(3,3) - t8655) - sin(qJ(3,3) + t8656) - sin(t8662) - sin(t8661)) * pkin(2)) / (-t8537 * sin(qJ(2,3) - qJ(3,3)) + pkin(2) * (0.2e1 * t8675 + pkin(2) * t8456 + (sin(t8526 + qJ(2,3)) - t8506) * pkin(3)));
t8457 = sin(t8498);
t8464 = qJ(2,2) + t8663;
t8465 = -qJ(2,2) + t8664;
t8509 = sin(qJ(2,2));
t8529 = 0.2e1 * qJ(3,2);
t8652 = ((-sin(t8465) - sin(t8464)) * t8654 + (cos(t8465) + cos(t8464)) * t8682 + (sin(0.2e1 * qJ(3,2) - t8657) - sin(t8529 + t8658) - 0.2e1 * t8510) * pkin(3) + (sin(qJ(3,2) - t8657) - sin(qJ(3,2) + t8658) - sin(t8664) - sin(t8663)) * pkin(2)) / (-t8537 * sin(qJ(2,2) - qJ(3,2)) + pkin(2) * (0.2e1 * t8673 + pkin(2) * t8457 + (sin(t8529 + qJ(2,2)) - t8509) * pkin(3)));
t8458 = sin(t8499);
t8466 = qJ(2,1) + t8665;
t8467 = -qJ(2,1) + t8666;
t8512 = sin(qJ(2,1));
t8532 = 0.2e1 * qJ(3,1);
t8651 = ((-sin(t8467) - sin(t8466)) * t8654 + (cos(t8467) + cos(t8466)) * t8682 + (sin(0.2e1 * qJ(3,1) - t8659) - sin(t8532 + t8660) - 0.2e1 * t8513) * pkin(3) + (sin(qJ(3,1) - t8659) - sin(qJ(3,1) + t8660) - sin(t8666) - sin(t8665)) * pkin(2)) / (-t8537 * sin(qJ(2,1) - qJ(3,1)) + pkin(2) * (0.2e1 * t8671 + pkin(2) * t8458 + (sin(t8532 + qJ(2,1)) - t8512) * pkin(3)));
t8602 = pkin(1) * t8516 + t8507 * t8487;
t8608 = t8505 * t8506;
t8399 = t8602 * t8608 + (t8491 - 0.1e1) * t8516 * pkin(3);
t8650 = t8399 * t8468;
t8649 = t8399 * t8471;
t8601 = pkin(1) * t8519 + t8510 * t8487;
t8607 = t8508 * t8509;
t8400 = t8601 * t8607 + (t8493 - 0.1e1) * t8519 * pkin(3);
t8648 = t8400 * t8469;
t8647 = t8400 * t8472;
t8600 = pkin(1) * t8522 + t8513 * t8487;
t8606 = t8511 * t8512;
t8401 = t8600 * t8606 + (t8495 - 0.1e1) * t8522 * pkin(3);
t8646 = t8401 * t8470;
t8645 = t8401 * t8473;
t8644 = t8556 * t8429;
t8643 = t8555 * t8430;
t8642 = t8554 * t8431;
t8592 = pkin(3) * t8608;
t8453 = t8474 + pkin(2);
t8624 = t8453 * t8515;
t8641 = 0.1e1 / (pkin(1) - t8592 + t8624) / t8505;
t8591 = pkin(3) * t8607;
t8454 = t8475 + pkin(2);
t8621 = t8454 * t8518;
t8640 = 0.1e1 / (pkin(1) - t8591 + t8621) / t8508;
t8590 = pkin(3) * t8606;
t8455 = t8476 + pkin(2);
t8618 = t8455 * t8521;
t8639 = 0.1e1 / (pkin(1) - t8590 + t8618) / t8511;
t8638 = t8429 * t8507;
t8636 = t8430 * t8510;
t8634 = t8431 * t8513;
t8525 = pkin(2) / 0.2e1;
t8629 = (t8474 + t8525) * t8505;
t8628 = (t8475 + t8525) * t8508;
t8627 = (t8476 + t8525) * t8511;
t8626 = t8453 * t8468;
t8625 = t8453 * t8471;
t8623 = t8454 * t8469;
t8622 = t8454 * t8472;
t8620 = t8455 * t8470;
t8619 = t8455 * t8473;
t8617 = t8468 * t8516;
t8616 = t8469 * t8519;
t8615 = t8470 * t8522;
t8614 = t8471 * t8516;
t8613 = t8472 * t8519;
t8612 = t8473 * t8522;
t8535 = pkin(3) ^ 2;
t8611 = t8491 * t8535;
t8610 = t8493 * t8535;
t8609 = t8495 * t8535;
t8444 = pkin(1) * t8506 - t8674;
t8605 = t8514 * t8444;
t8445 = pkin(1) * t8509 - t8672;
t8604 = t8517 * t8445;
t8446 = pkin(1) * t8512 - t8670;
t8603 = t8520 * t8446;
t8599 = -t8535 / 0.2e1 + t8537 / 0.2e1;
t8598 = pkin(2) * t8474;
t8597 = pkin(2) * t8475;
t8596 = pkin(2) * t8476;
t8595 = t8453 * t8674;
t8594 = t8454 * t8672;
t8593 = t8455 * t8670;
t8432 = t8468 * g(1) + t8471 * g(2);
t8384 = -t8432 * t8459 + t8456 * t8684;
t8589 = t8384 * t8641;
t8385 = t8432 * t8456 + t8459 * t8684;
t8588 = t8385 * t8641;
t8433 = t8469 * g(1) + t8472 * g(2);
t8386 = -t8433 * t8460 + t8457 * t8685;
t8587 = t8386 * t8640;
t8387 = t8433 * t8457 + t8460 * t8685;
t8586 = t8387 * t8640;
t8434 = t8470 * g(1) + t8473 * g(2);
t8388 = -t8434 * t8461 + t8458 * t8686;
t8585 = t8388 * t8639;
t8389 = t8434 * t8458 + t8461 * t8686;
t8584 = t8389 * t8639;
t8390 = -t8432 * t8515 + t8506 * t8684;
t8583 = t8390 * t8641;
t8391 = t8432 * t8506 + t8515 * t8684;
t8582 = t8391 * t8641;
t8392 = -t8433 * t8518 + t8509 * t8685;
t8581 = t8392 * t8640;
t8393 = t8433 * t8509 + t8518 * t8685;
t8580 = t8393 * t8640;
t8394 = -t8434 * t8521 + t8512 * t8686;
t8579 = t8394 * t8639;
t8395 = t8434 * t8512 + t8521 * t8686;
t8578 = t8395 * t8639;
t8577 = t8556 * t8638;
t8575 = t8555 * t8636;
t8573 = t8554 * t8634;
t8571 = t8456 * t8644;
t8570 = t8459 * t8644;
t8569 = t8468 * t8638;
t8568 = t8471 * t8638;
t8567 = t8457 * t8643;
t8566 = t8460 * t8643;
t8565 = t8469 * t8636;
t8564 = t8472 * t8636;
t8563 = t8458 * t8642;
t8562 = t8461 * t8642;
t8561 = t8470 * t8634;
t8560 = t8473 * t8634;
t8559 = t8516 * t8608;
t8558 = t8519 * t8607;
t8557 = t8522 * t8606;
t8553 = t8506 * t8577;
t8552 = t8515 * t8577;
t8551 = t8509 * t8575;
t8550 = t8518 * t8575;
t8549 = t8512 * t8573;
t8548 = t8521 * t8573;
t8547 = t8507 * t8571;
t8546 = t8507 * t8570;
t8545 = t8510 * t8567;
t8544 = t8510 * t8566;
t8543 = t8513 * t8563;
t8542 = t8513 * t8562;
t8414 = t8559 * t8683 - t8602;
t8541 = pkin(2) * t8559 + t8414 * t8514;
t8415 = t8558 * t8683 - t8601;
t8540 = pkin(2) * t8558 + t8415 * t8517;
t8416 = t8557 * t8683 - t8600;
t8539 = pkin(2) * t8557 + t8416 * t8520;
t8538 = 0.1e1 / pkin(2);
t8536 = 0.1e1 / pkin(3);
t8524 = -pkin(3) / 0.2e1;
t8486 = -t8535 + t8537;
t8440 = t8482 + t8667 / 0.2e1 + t8524;
t8439 = t8481 + t8668 / 0.2e1 + t8524;
t8438 = t8480 + t8669 / 0.2e1 + t8524;
t8428 = t8596 + t8599 + t8609;
t8427 = t8597 + t8599 + t8610;
t8426 = t8598 + t8599 + t8611;
t8407 = t8671 + (-pkin(3) + t8667 + 0.2e1 * t8482) * t8512;
t8406 = t8673 + (-pkin(3) + t8668 + 0.2e1 * t8481) * t8509;
t8405 = t8675 + (-pkin(3) + t8669 + 0.2e1 * t8480) * t8506;
t8404 = pkin(1) * t8670 + (t8486 + 0.2e1 * t8596 + 0.2e1 * t8609) * t8512;
t8403 = pkin(1) * t8672 + (t8486 + 0.2e1 * t8597 + 0.2e1 * t8610) * t8509;
t8402 = pkin(1) * t8674 + (t8486 + 0.2e1 * t8598 + 0.2e1 * t8611) * t8506;
t8383 = t8428 * t8513 * t8679 + ((pkin(1) - 0.2e1 * t8590) * t8513 - t8522 * t8487) * t8618 - pkin(3) * ((pkin(1) * t8606 - pkin(3) + t8482) * t8513 - t8487 * t8557);
t8382 = t8427 * t8510 * t8680 + ((pkin(1) - 0.2e1 * t8591) * t8510 - t8519 * t8487) * t8621 - pkin(3) * ((pkin(1) * t8607 - pkin(3) + t8481) * t8510 - t8487 * t8558);
t8381 = t8426 * t8507 * t8681 + ((pkin(1) - 0.2e1 * t8592) * t8507 - t8516 * t8487) * t8624 - pkin(3) * ((pkin(1) * t8608 - pkin(3) + t8480) * t8507 - t8487 * t8559);
t8377 = (t8440 * t8612 + t8470 * t8627) * t8679 + (t8470 * t8407 - t8473 * t8539) * t8521 - t8645 + t8470 * t8603;
t8376 = (-t8440 * t8615 + t8473 * t8627) * t8679 + (t8473 * t8407 + t8470 * t8539) * t8521 + t8646 + t8473 * t8603;
t8375 = (t8439 * t8613 + t8469 * t8628) * t8680 + (t8469 * t8406 - t8472 * t8540) * t8518 - t8647 + t8469 * t8604;
t8374 = (-t8439 * t8616 + t8472 * t8628) * t8680 + (t8472 * t8406 + t8469 * t8540) * t8518 + t8648 + t8472 * t8604;
t8373 = (t8438 * t8614 + t8468 * t8629) * t8681 + (t8468 * t8405 - t8471 * t8541) * t8515 - t8649 + t8468 * t8605;
t8372 = (-t8438 * t8617 + t8471 * t8629) * t8681 + (t8471 * t8405 + t8468 * t8541) * t8515 + t8650 + t8471 * t8605;
t8371 = (-t8428 * t8612 - t8470 * t8593) * t8679 + (-t8470 * t8404 + t8416 * t8619) * t8521 + pkin(3) * t8645 - t8446 * t8620;
t8370 = (t8428 * t8615 - t8473 * t8593) * t8679 + (-t8473 * t8404 - t8416 * t8620) * t8521 - pkin(3) * t8646 - t8446 * t8619;
t8369 = (-t8427 * t8613 - t8469 * t8594) * t8680 + (-t8469 * t8403 + t8415 * t8622) * t8518 + pkin(3) * t8647 - t8445 * t8623;
t8368 = (t8427 * t8616 - t8472 * t8594) * t8680 + (-t8472 * t8403 - t8415 * t8623) * t8518 - pkin(3) * t8648 - t8445 * t8622;
t8367 = (-t8426 * t8614 - t8468 * t8595) * t8681 + (-t8468 * t8402 + t8414 * t8625) * t8515 + pkin(3) * t8649 - t8444 * t8626;
t8366 = (t8426 * t8617 - t8471 * t8595) * t8681 + (-t8471 * t8402 - t8414 * t8626) * t8515 - pkin(3) * t8650 - t8444 * t8625;
t1 = [(-t8554 * t8560 - t8555 * t8564 - t8556 * t8568) * MDP(2) + (-t8560 * t8686 - t8564 * t8685 - t8568 * t8684) * MDP(3) + (-t8471 * t8552 - t8472 * t8550 - t8473 * t8548) * MDP(9) + (t8471 * t8553 + t8472 * t8551 + t8473 * t8549) * MDP(10) + (-t8471 * t8546 - t8472 * t8544 - t8473 * t8542) * MDP(16) + (t8471 * t8547 + t8472 * t8545 + t8473 * t8543) * MDP(17) - g(1) * MDP(18) + ((t8373 * t8583 + t8375 * t8581 + t8377 * t8579) * MDP(9) + (t8373 * t8582 + t8375 * t8580 + t8377 * t8578) * MDP(10) + (t8373 * t8589 + t8375 * t8587 + t8377 * t8585) * MDP(16) + (t8373 * t8588 + t8375 * t8586 + t8377 * t8584) * MDP(17) + ((t8367 * t8589 + t8369 * t8587 + t8371 * t8585) * MDP(16) + (t8367 * t8588 + t8369 * t8586 + t8371 * t8584) * MDP(17)) * t8536) * t8538; (t8554 * t8561 + t8555 * t8565 + t8556 * t8569) * MDP(2) + (t8561 * t8686 + t8565 * t8685 + t8569 * t8684) * MDP(3) + (t8468 * t8552 + t8469 * t8550 + t8470 * t8548) * MDP(9) + (-t8468 * t8553 - t8469 * t8551 - t8470 * t8549) * MDP(10) + (t8468 * t8546 + t8469 * t8544 + t8470 * t8542) * MDP(16) + (-t8468 * t8547 - t8469 * t8545 - t8470 * t8543) * MDP(17) - g(2) * MDP(18) + ((t8372 * t8583 + t8374 * t8581 + t8376 * t8579) * MDP(9) + (t8372 * t8582 + t8374 * t8580 + t8376 * t8578) * MDP(10) + (t8372 * t8589 + t8374 * t8587 + t8376 * t8585) * MDP(16) + (t8372 * t8588 + t8374 * t8586 + t8376 * t8584) * MDP(17) + ((t8366 * t8589 + t8368 * t8587 + t8370 * t8585) * MDP(16) + (t8366 * t8588 + t8368 * t8586 + t8370 * t8584) * MDP(17)) * t8536) * t8538; (-t8574 - t8576 - t8572) * MDP(2) + (-t8633 * t8686 - t8635 * t8685 - t8637 * t8684) * MDP(3) + (-t8515 * t8576 - t8518 * t8574 - t8521 * t8572) * MDP(9) + (t8506 * t8576 + t8509 * t8574 + t8512 * t8572) * MDP(10) + (-t8516 * t8570 - t8519 * t8566 - t8522 * t8562) * MDP(16) + (t8516 * t8571 + t8519 * t8567 + t8522 * t8563) * MDP(17) - g(3) * MDP(18) + ((t8381 * t8589 + t8382 * t8587 + t8383 * t8585) * MDP(16) + (t8381 * t8588 + t8382 * t8586 + t8383 * t8584) * MDP(17)) * t8538 * t8536 + (t8390 * t8653 + t8392 * t8652 + t8394 * t8651) * MDP(9) / 0.2e1 + (t8391 * t8653 + t8393 * t8652 + t8395 * t8651) * MDP(10) / 0.2e1 + (t8384 * t8653 + t8386 * t8652 + t8388 * t8651) * MDP(16) / 0.2e1 + (t8385 * t8653 + t8387 * t8652 + t8389 * t8651) * MDP(17) / 0.2e1;];
taugX  = t1;
