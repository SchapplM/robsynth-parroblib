% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR6V1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x12]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:32
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR6V1G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR6V1G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:32:18
% EndTime: 2020-08-06 18:32:20
% DurationCPUTime: 1.89s
% Computational Cost: add. (7092->250), mult. (5442->455), div. (516->14), fcn. (3096->65), ass. (0->248)
t550 = sin(qJ(3,1));
t559 = 0.2e1 * qJ(3,1);
t702 = 2 * pkin(2);
t481 = pkin(3) * sin(t559) + t550 * t702 + (sin(pkin(7) + qJ(3,1)) + sin(-pkin(7) + qJ(3,1))) * pkin(1);
t478 = 0.1e1 / t481 ^ 2;
t553 = cos(qJ(3,1));
t708 = t478 * t553 / 0.2e1;
t707 = -t478 * t550 / 0.2e1;
t549 = sin(qJ(3,2));
t558 = 0.2e1 * qJ(3,2);
t480 = pkin(3) * sin(t558) + t549 * t702 + (sin(pkin(7) + qJ(3,2)) + sin(-pkin(7) + qJ(3,2))) * pkin(1);
t476 = 0.1e1 / t480 ^ 2;
t552 = cos(qJ(3,2));
t706 = t476 * t552 / 0.2e1;
t705 = -t476 * t549 / 0.2e1;
t548 = sin(qJ(3,3));
t557 = 0.2e1 * qJ(3,3);
t479 = pkin(3) * sin(t557) + t548 * t702 + (sin(pkin(7) + qJ(3,3)) + sin(-pkin(7) + qJ(3,3))) * pkin(1);
t474 = 0.1e1 / t479 ^ 2;
t551 = cos(qJ(3,3));
t704 = t474 * t551 / 0.2e1;
t703 = -t474 * t548 / 0.2e1;
t473 = 0.1e1 / t479;
t659 = t473 / 0.2e1;
t701 = t474 / 0.4e1;
t475 = 0.1e1 / t480;
t654 = t475 / 0.2e1;
t700 = t476 / 0.4e1;
t477 = 0.1e1 / t481;
t649 = t477 / 0.2e1;
t699 = t478 / 0.4e1;
t540 = qJ(1,1) + legFrame(1,3);
t529 = pkin(7) + t540;
t511 = sin(t529);
t524 = qJ(3,1) + t529;
t525 = -qJ(3,1) + t529;
t622 = cos(t524) + cos(t525);
t556 = (pkin(6) + pkin(5));
t669 = -2 * t556;
t514 = cos(t529);
t671 = -0.2e1 * t514;
t683 = -0.2e1 * pkin(1);
t472 = t511 * t669 + cos(t540) * t683 + pkin(2) * t671 - t622 * pkin(3);
t698 = t472 * t708;
t697 = t472 * t707;
t539 = qJ(1,2) + legFrame(2,3);
t528 = pkin(7) + t539;
t510 = sin(t528);
t520 = qJ(3,2) + t528;
t521 = -qJ(3,2) + t528;
t623 = cos(t520) + cos(t521);
t513 = cos(t528);
t673 = -0.2e1 * t513;
t471 = t510 * t669 + cos(t539) * t683 + pkin(2) * t673 - t623 * pkin(3);
t696 = t471 * t706;
t695 = t471 * t705;
t538 = qJ(1,3) + legFrame(3,3);
t527 = pkin(7) + t538;
t509 = sin(t527);
t516 = qJ(3,3) + t527;
t517 = -qJ(3,3) + t527;
t624 = cos(t516) + cos(t517);
t512 = cos(t527);
t675 = -0.2e1 * t512;
t470 = t509 * t669 + cos(t538) * t683 + pkin(2) * t675 - t624 * pkin(3);
t694 = t470 * t704;
t693 = t470 * t703;
t625 = sin(t524) + sin(t525);
t668 = 2 * t556;
t677 = -0.2e1 * t511;
t469 = t514 * t668 + sin(t540) * t683 + pkin(2) * t677 - t625 * pkin(3);
t692 = t469 * t708;
t691 = t469 * t707;
t626 = sin(t520) + sin(t521);
t679 = -0.2e1 * t510;
t468 = t513 * t668 + sin(t539) * t683 + pkin(2) * t679 - t626 * pkin(3);
t690 = t468 * t706;
t689 = t468 * t705;
t627 = sin(t516) + sin(t517);
t681 = -0.2e1 * t509;
t467 = t512 * t668 + sin(t538) * t683 + pkin(2) * t681 - t627 * pkin(3);
t688 = t467 * t704;
t687 = t467 * t703;
t686 = t548 * t551;
t685 = t549 * t552;
t684 = t550 * t553;
t680 = 0.2e1 * t509;
t678 = 0.2e1 * t510;
t676 = 0.2e1 * t511;
t674 = 0.2e1 * t512;
t672 = 0.2e1 * t513;
t670 = 0.2e1 * t514;
t537 = cos(pkin(7)) * pkin(1) + pkin(2);
t667 = t467 * t473;
t666 = t468 * t475;
t665 = t469 * t477;
t664 = t470 * t473;
t663 = t471 * t475;
t662 = t472 * t477;
t661 = t473 * t548;
t660 = t473 * t551;
t656 = t475 * t549;
t655 = t475 * t552;
t651 = t477 * t550;
t650 = t477 * t553;
t488 = t551 * pkin(3) + t537;
t482 = 0.1e1 / t488;
t646 = t482 * t509;
t645 = t482 * t512;
t644 = t482 * t548;
t643 = t482 * t551;
t489 = t552 * pkin(3) + t537;
t484 = 0.1e1 / t489;
t642 = t484 * t510;
t641 = t484 * t513;
t640 = t484 * t549;
t639 = t484 * t552;
t490 = t553 * pkin(3) + t537;
t486 = 0.1e1 / t490;
t638 = t486 * t511;
t637 = t486 * t514;
t636 = t486 * t550;
t635 = t486 * t553;
t483 = 0.1e1 / t488 ^ 2;
t503 = t509 ^ 2;
t634 = t503 * t483;
t485 = 0.1e1 / t489 ^ 2;
t504 = t510 ^ 2;
t633 = t504 * t485;
t487 = 0.1e1 / t490 ^ 2;
t505 = t511 ^ 2;
t632 = t505 * t487;
t506 = t512 ^ 2;
t631 = t506 * t483;
t507 = t513 ^ 2;
t630 = t507 * t485;
t508 = t514 ^ 2;
t629 = t508 * t487;
t536 = sin(pkin(7)) * pkin(1) + pkin(5);
t560 = 0.1e1 / pkin(3);
t628 = t536 * t560;
t621 = 0.2e1 * pkin(1);
t619 = 0.2e1 * t560;
t606 = t473 * t628;
t605 = t548 * t659;
t604 = t551 * t659;
t603 = t475 * t628;
t602 = t549 * t654;
t601 = t552 * t654;
t600 = t477 * t628;
t599 = t550 * t649;
t598 = t553 * t649;
t597 = t536 * t644;
t596 = t536 * t643;
t595 = t537 * t644;
t594 = t537 * t643;
t593 = t483 * t686;
t592 = t536 * t640;
t591 = t536 * t639;
t590 = t537 * t640;
t589 = t537 * t639;
t588 = t485 * t685;
t587 = t536 * t636;
t586 = t536 * t635;
t585 = t537 * t636;
t584 = t537 * t635;
t583 = t487 * t684;
t582 = t509 * t483 * t512;
t581 = t510 * t485 * t513;
t580 = t511 * t487 * t514;
t579 = t645 * t667;
t578 = t641 * t666;
t577 = t637 * t665;
t576 = t646 * t664;
t575 = t642 * t663;
t574 = t638 * t662;
t573 = t548 * t606;
t572 = t551 * t606;
t571 = t549 * t603;
t570 = t552 * t603;
t569 = t550 * t600;
t568 = t553 * t600;
t567 = t473 * t482 * (-t467 * t509 + t470 * t512);
t566 = t475 * t484 * (-t468 * t510 + t471 * t513);
t565 = t477 * t486 * (-t469 * t511 + t472 * t514);
t564 = t632 + t633 + t634;
t563 = t629 + t630 + t631;
t466 = -t580 - t581 - t582;
t562 = pkin(1) ^ 2;
t561 = 0.1e1 / pkin(3) ^ 2;
t544 = t550 ^ 2;
t543 = t549 ^ 2;
t542 = t548 ^ 2;
t535 = -qJ(3,1) + t540;
t534 = qJ(3,1) + t540;
t533 = -qJ(3,2) + t539;
t532 = qJ(3,2) + t539;
t531 = -qJ(3,3) + t538;
t530 = qJ(3,3) + t538;
t526 = -0.2e1 * qJ(3,1) + t529;
t523 = t559 + t529;
t522 = -0.2e1 * qJ(3,2) + t528;
t519 = t558 + t528;
t518 = -0.2e1 * qJ(3,3) + t527;
t515 = t557 + t527;
t465 = -t542 * t582 - t543 * t581 - t544 * t580;
t464 = -0.2e1 * t580 * t684 - 0.2e1 * t581 * t685 - 0.2e1 * t582 * t686;
t463 = t625 * t702 + (sin(t535) + sin(t534)) * t621 + t622 * t669 + (sin(t526) + sin(t523) + t676) * pkin(3);
t462 = t626 * t702 + (sin(t533) + sin(t532)) * t621 + t623 * t669 + (sin(t522) + sin(t519) + t678) * pkin(3);
t461 = t627 * t702 + (sin(t531) + sin(t530)) * t621 + t624 * t669 + (sin(t518) + sin(t515) + t680) * pkin(3);
t460 = t622 * t702 + (cos(t535) + cos(t534)) * t621 + t625 * t668 + (cos(t526) + cos(t523) + t670) * pkin(3);
t459 = t623 * t702 + (cos(t533) + cos(t532)) * t621 + t626 * t668 + (cos(t522) + cos(t519) + t672) * pkin(3);
t458 = t624 * t702 + (cos(t531) + cos(t530)) * t621 + t627 * t668 + (cos(t518) + cos(t515) + t674) * pkin(3);
t457 = -t469 * t569 + t584 * t670;
t456 = -t472 * t569 + t584 * t677;
t455 = -t467 * t573 + t594 * t674;
t454 = -t470 * t573 + t594 * t681;
t453 = -t468 * t571 + t589 * t672;
t452 = -t471 * t571 + t589 * t679;
t451 = -t472 * t568 + t585 * t676;
t450 = -t469 * t568 + t585 * t671;
t449 = -t471 * t570 + t590 * t678;
t448 = -t468 * t570 + t590 * t673;
t447 = -t470 * t572 + t595 * t680;
t446 = -t467 * t572 + t595 * t675;
t445 = -t463 * t599 - t514 * t586;
t444 = -t460 * t599 + t511 * t586;
t443 = -t462 * t602 - t513 * t591;
t442 = t462 * t601 - t513 * t592;
t441 = -t461 * t605 - t512 * t596;
t440 = t461 * t604 - t512 * t597;
t439 = t463 * t598 - t514 * t587;
t438 = t460 * t598 + t511 * t587;
t437 = -t459 * t602 + t510 * t591;
t436 = t459 * t601 + t510 * t592;
t435 = -t458 * t605 + t509 * t596;
t434 = t458 * t604 + t509 * t597;
t433 = (t470 * t660 + t471 * t655 + t472 * t650) * t560;
t432 = (-t470 * t661 - t471 * t656 - t472 * t651) * t560;
t431 = (t467 * t660 + t468 * t655 + t469 * t650) * t560;
t430 = (-t467 * t661 - t468 * t656 - t469 * t651) * t560;
t429 = t461 * t659 + t462 * t654 + t463 * t649;
t428 = t458 * t659 + t459 * t654 + t460 * t649;
t427 = (t467 * t470 * t474 + t468 * t471 * t476 + t469 * t472 * t478) * t561;
t426 = (t551 * t567 + t552 * t566 + t553 * t565) * t560;
t425 = (t548 * t567 + t549 * t566 + t550 * t565) * t560;
t424 = t458 * t461 * t701 + t459 * t462 * t700 + t460 * t463 * t699 + t466 * t562;
t1 = [t564, 0, 0, t458 ^ 2 * t701 + t459 ^ 2 * t700 + t460 ^ 2 * t699 + t564 * t562, t542 * t634 + t543 * t633 + t544 * t632, 0.2e1 * t503 * t593 + 0.2e1 * t504 * t588 + 0.2e1 * t505 * t583, (-t548 * t576 - t549 * t575 - t550 * t574) * t619, (-t551 * t576 - t552 * t575 - t553 * t574) * t619, (t470 ^ 2 * t474 + t471 ^ 2 * t476 + t472 ^ 2 * t478) * t561, -t452 * t642 - t454 * t646 - t456 * t638 + (t434 * t664 + t436 * t663 + t438 * t662 + t458 * t694 + t459 * t696 + t460 * t698) * t560, -t447 * t646 - t449 * t642 - t451 * t638 + (t435 * t664 + t437 * t663 + t444 * t662 + t458 * t693 + t459 * t695 + t460 * t697) * t560, 1; t466, 0, 0, t424, t465, t464, t425, t426, t427, -t453 * t642 - t455 * t646 - t457 * t638 + (t439 * t662 + t440 * t664 + t442 * t663 + t458 * t688 + t459 * t690 + t460 * t692) * t560, -t446 * t646 - t448 * t642 - t450 * t638 + (t441 * t664 + t443 * t663 + t445 * t662 + t458 * t687 + t459 * t689 + t460 * t691) * t560, 0; 0, 0, 0, t428, 0, 0, 0, 0, 0, t433, t432, 0; t466, 0, 0, t424, t465, t464, t425, t426, t427, t452 * t641 + t454 * t645 + t456 * t637 + (t434 * t667 + t436 * t666 + t438 * t665 + t461 * t694 + t462 * t696 + t463 * t698) * t560, t447 * t645 + t449 * t641 + t451 * t637 + (t435 * t667 + t437 * t666 + t444 * t665 + t461 * t693 + t462 * t695 + t463 * t697) * t560, 0; t563, 0, 0, t461 ^ 2 * t701 + t462 ^ 2 * t700 + t463 ^ 2 * t699 + t563 * t562, t542 * t631 + t543 * t630 + t544 * t629, 0.2e1 * t506 * t593 + 0.2e1 * t507 * t588 + 0.2e1 * t508 * t583, (t548 * t579 + t549 * t578 + t550 * t577) * t619, (t551 * t579 + t552 * t578 + t553 * t577) * t619, (t467 ^ 2 * t474 + t468 ^ 2 * t476 + t469 ^ 2 * t478) * t561, t453 * t641 + t455 * t645 + t457 * t637 + (t439 * t665 + t440 * t667 + t442 * t666 + t461 * t688 + t462 * t690 + t463 * t692) * t560, t446 * t645 + t448 * t641 + t450 * t637 + (t441 * t667 + t443 * t666 + t445 * t665 + t461 * t687 + t462 * t689 + t463 * t691) * t560, 1; 0, 0, 0, t429, 0, 0, 0, 0, 0, t431, t430, 0; 0, 0, 0, t428, 0, 0, 0, 0, 0, t433, t432, 0; 0, 0, 0, t429, 0, 0, 0, 0, 0, t431, t430, 0; 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0, 1;];
tau_reg  = t1;
