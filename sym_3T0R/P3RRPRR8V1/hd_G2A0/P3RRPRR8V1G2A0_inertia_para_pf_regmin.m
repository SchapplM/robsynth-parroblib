% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x13]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 19:59
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RRPRR8V1G2A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: pkin has to be [5x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:59:31
% EndTime: 2020-08-06 19:59:34
% DurationCPUTime: 2.29s
% Computational Cost: add. (3313->238), mult. (6336->584), div. (807->18), fcn. (5793->23), ass. (0->268)
t702 = 2 * qJ(3,1);
t701 = 2 * qJ(3,2);
t700 = 2 * qJ(3,3);
t517 = cos(pkin(5)) * pkin(2) + pkin(1);
t546 = legFrame(3,2);
t521 = sin(t546);
t524 = cos(t546);
t549 = sin(qJ(2,3));
t555 = cos(qJ(2,3));
t550 = sin(qJ(1,3));
t690 = pkin(2) * sin(pkin(5));
t653 = t550 * t690;
t668 = t517 * t550;
t478 = (-t521 * t668 + t524 * t690) * t555 + t549 * (t517 * t524 + t521 * t653);
t481 = (t521 * t690 + t524 * t668) * t555 + (t517 * t521 - t524 * t653) * t549;
t699 = t478 * t481;
t547 = legFrame(2,2);
t522 = sin(t547);
t525 = cos(t547);
t551 = sin(qJ(2,2));
t557 = cos(qJ(2,2));
t552 = sin(qJ(1,2));
t652 = t552 * t690;
t667 = t517 * t552;
t479 = (-t522 * t667 + t525 * t690) * t557 + t551 * (t517 * t525 + t522 * t652);
t482 = (t522 * t690 + t525 * t667) * t557 + (t517 * t522 - t525 * t652) * t551;
t698 = t479 * t482;
t548 = legFrame(1,2);
t523 = sin(t548);
t526 = cos(t548);
t553 = sin(qJ(2,1));
t559 = cos(qJ(2,1));
t554 = sin(qJ(1,1));
t651 = t554 * t690;
t666 = t517 * t554;
t480 = (-t523 * t666 + t526 * t690) * t559 + t553 * (t517 * t526 + t523 * t651);
t483 = (t523 * t690 + t526 * t666) * t559 + (t517 * t523 - t526 * t651) * t553;
t697 = t480 * t483;
t505 = t517 * t555 - t549 * t690;
t496 = 0.1e1 / t505;
t556 = cos(qJ(1,3));
t696 = t496 * t556;
t506 = t517 * t557 - t551 * t690;
t498 = 0.1e1 / t506;
t558 = cos(qJ(1,2));
t695 = t498 * t558;
t507 = t517 * t559 - t553 * t690;
t500 = 0.1e1 / t507;
t560 = cos(qJ(1,1));
t694 = t500 * t560;
t693 = pkin(1) * t555;
t692 = pkin(1) * t557;
t691 = pkin(1) * t559;
t689 = t478 * t496;
t688 = t479 * t498;
t687 = t480 * t500;
t686 = t481 * t496;
t685 = t482 * t498;
t684 = t483 * t500;
t514 = pkin(2) * cos(qJ(2,3) + pkin(5)) + t693;
t543 = pkin(4) + qJ(3,3);
t493 = t514 * t556 + t543 * t550;
t527 = 1 / t543;
t564 = pkin(1) ^ 2;
t623 = t555 ^ 2 * t564 + (qJ(3,3) ^ 2);
t484 = (-t493 * t693 + t623 * t556) * t527;
t683 = t484 * t496;
t515 = pkin(2) * cos(qJ(2,2) + pkin(5)) + t692;
t544 = pkin(4) + qJ(3,2);
t494 = t515 * t558 + t544 * t552;
t529 = 1 / t544;
t622 = t557 ^ 2 * t564 + (qJ(3,2) ^ 2);
t485 = (-t494 * t692 + t622 * t558) * t529;
t682 = t485 * t498;
t516 = pkin(2) * cos(qJ(2,1) + pkin(5)) + t691;
t545 = pkin(4) + qJ(3,1);
t495 = t516 * t560 + t545 * t554;
t531 = 1 / t545;
t621 = t559 ^ 2 * t564 + (qJ(3,1) ^ 2);
t486 = (-t495 * t691 + t621 * t560) * t531;
t681 = t486 * t500;
t680 = t496 * t527;
t497 = 0.1e1 / t505 ^ 2;
t528 = 1 / t543 ^ 2;
t679 = t497 * t528;
t678 = t498 * t529;
t499 = 0.1e1 / t506 ^ 2;
t530 = 1 / t544 ^ 2;
t677 = t499 * t530;
t676 = t500 * t531;
t501 = 0.1e1 / t507 ^ 2;
t532 = 1 / t545 ^ 2;
t675 = t501 * t532;
t508 = 0.1e1 / t514;
t674 = t508 * t521;
t673 = t508 * t524;
t510 = 0.1e1 / t515;
t672 = t510 * t522;
t671 = t510 * t525;
t512 = 0.1e1 / t516;
t670 = t512 * t523;
t669 = t512 * t526;
t665 = t527 * t556;
t537 = t556 ^ 2;
t664 = t528 * t537;
t663 = t529 * t558;
t539 = t558 ^ 2;
t662 = t530 * t539;
t661 = t531 * t560;
t541 = t560 ^ 2;
t660 = t532 * t541;
t659 = t496 * t693;
t658 = t498 * t692;
t657 = t500 * t691;
t656 = pkin(1) * t508 * t549;
t655 = pkin(1) * t510 * t551;
t654 = pkin(1) * t512 * t553;
t650 = t478 * t680;
t649 = t479 * t678;
t648 = t480 * t676;
t647 = t481 * t680;
t646 = t482 * t678;
t645 = t483 * t676;
t644 = t508 * t680;
t643 = t549 * t680;
t642 = t528 * t696;
t533 = t549 ^ 2;
t641 = t533 * t679;
t640 = t510 * t678;
t639 = t551 * t678;
t638 = t530 * t695;
t534 = t551 ^ 2;
t637 = t534 * t677;
t636 = t512 * t676;
t635 = t553 * t676;
t634 = t532 * t694;
t535 = t553 ^ 2;
t633 = t535 * t675;
t632 = t508 * t665;
t631 = t510 * t663;
t630 = t512 * t661;
t629 = t528 * t549 * t555;
t628 = t530 * t551 * t557;
t627 = t532 * t553 * t559;
t626 = t676 * t702;
t625 = t678 * t701;
t624 = t680 * t700;
t620 = t521 * t656;
t619 = t524 * t656;
t618 = t522 * t655;
t617 = t525 * t655;
t616 = t523 * t654;
t615 = t526 * t654;
t614 = qJ(3,1) * t635;
t613 = qJ(3,2) * t639;
t612 = qJ(3,3) * t643;
t611 = t679 * t699;
t610 = t677 * t698;
t609 = t675 * t697;
t608 = t508 * t643;
t607 = t555 * t644;
t606 = t533 * t642;
t605 = t497 * t629;
t604 = t510 * t639;
t603 = t557 * t640;
t602 = t534 * t638;
t601 = t499 * t628;
t600 = t512 * t635;
t599 = t559 * t636;
t598 = t535 * t634;
t597 = t501 * t627;
t596 = t549 * t632;
t595 = t555 * t632;
t594 = t551 * t631;
t593 = t557 * t631;
t592 = t553 * t630;
t591 = t559 * t630;
t590 = t634 * t702;
t589 = t638 * t701;
t588 = t642 * t700;
t587 = t521 * t608;
t586 = t524 * t608;
t585 = t629 * t696;
t584 = t522 * t604;
t583 = t525 * t604;
t582 = t628 * t695;
t581 = t523 * t600;
t580 = t526 * t600;
t579 = t627 * t694;
t578 = t521 * t596;
t577 = t524 * t596;
t576 = t522 * t594;
t575 = t525 * t594;
t574 = t523 * t592;
t573 = t526 * t592;
t572 = t623 * t496;
t571 = t622 * t498;
t570 = t621 * t500;
t569 = (t478 * t521 + t481 * t524) * t644;
t568 = (t479 * t522 + t482 * t525) * t640;
t567 = (t480 * t523 + t483 * t526) * t636;
t461 = t574 + t576 + t578;
t463 = t573 + t575 + t577;
t566 = t478 * t586 + t479 * t583 + t480 * t580;
t565 = t481 * t587 + t482 * t584 + t483 * t581;
t513 = 0.1e1 / t516 ^ 2;
t511 = 0.1e1 / t515 ^ 2;
t509 = 0.1e1 / t514 ^ 2;
t504 = t553 * t517 + t559 * t690;
t503 = t551 * t517 + t557 * t690;
t502 = t549 * t517 + t555 * t690;
t492 = t507 * t554 - t545 * t560;
t491 = t506 * t552 - t544 * t558;
t490 = t505 * t550 - t543 * t556;
t489 = (-t560 * t691 + t495) * t531;
t488 = (-t558 * t692 + t494) * t529;
t487 = (-t556 * t693 + t493) * t527;
t477 = t483 ^ 2;
t476 = t482 ^ 2;
t475 = t481 ^ 2;
t474 = t480 ^ 2;
t473 = t479 ^ 2;
t472 = t478 ^ 2;
t471 = -t492 * t523 + t504 * t526;
t470 = t492 * t526 + t504 * t523;
t469 = -t491 * t522 + t503 * t525;
t468 = t491 * t525 + t503 * t522;
t467 = -t490 * t521 + t502 * t524;
t466 = t490 * t524 + t502 * t521;
t465 = t509 * t521 * t524 + t511 * t522 * t525 + t513 * t523 * t526;
t464 = t524 * t595 + t525 * t593 + t526 * t591;
t462 = t521 * t595 + t522 * t593 + t523 * t591;
t460 = t483 * t626 - t616;
t459 = t480 * t626 - t615;
t458 = t482 * t625 - t618;
t457 = t479 * t625 - t617;
t456 = t481 * t624 - t620;
t455 = t478 * t624 - t619;
t454 = pkin(1) * t670 - t483 * t614;
t453 = pkin(1) * t669 - t480 * t614;
t452 = pkin(1) * t672 - t482 * t613;
t451 = pkin(1) * t671 - t479 * t613;
t450 = pkin(1) * t674 - t481 * t612;
t449 = pkin(1) * t673 - t478 * t612;
t448 = (-t483 * t657 + t470) * t531;
t447 = (-t480 * t657 + t471) * t531;
t446 = (-t482 * t658 + t468) * t529;
t445 = (-t479 * t658 + t469) * t529;
t444 = (-t481 * t659 + t466) * t527;
t443 = (-t478 * t659 + t467) * t527;
t442 = t481 * t642 + t482 * t638 + t483 * t634;
t441 = t478 * t642 + t479 * t638 + t480 * t634;
t440 = -qJ(3,1) * t616 + (-t470 * t691 + t483 * t570) * t531;
t439 = -qJ(3,1) * t615 + (-t471 * t691 + t480 * t570) * t531;
t438 = -qJ(3,2) * t618 + (-t468 * t692 + t482 * t571) * t529;
t437 = -qJ(3,2) * t617 + (-t469 * t692 + t479 * t571) * t529;
t436 = -qJ(3,3) * t620 + (-t466 * t693 + t481 * t572) * t527;
t435 = -qJ(3,3) * t619 + (-t467 * t693 + t478 * t572) * t527;
t434 = t481 * t606 + t482 * t602 + t483 * t598;
t433 = t478 * t606 + t479 * t602 + t480 * t598;
t432 = 0.2e1 * t481 * t585 + 0.2e1 * t482 * t582 + 0.2e1 * t483 * t579;
t431 = 0.2e1 * t478 * t585 + 0.2e1 * t479 * t582 + 0.2e1 * t480 * t579;
t430 = t609 + t610 + t611;
t429 = t533 * t611 + t534 * t610 + t535 * t609;
t428 = 0.2e1 * t597 * t697 + 0.2e1 * t601 * t698 + 0.2e1 * t605 * t699;
t427 = t555 * t569 + t557 * t568 + t559 * t567;
t426 = t549 * t569 + t551 * t568 + t553 * t567;
t1 = [t475 * t679 + t476 * t677 + t477 * t675, 0, 0, t475 * t641 + t476 * t637 + t477 * t633, 0.2e1 * t475 * t605 + 0.2e1 * t476 * t601 + 0.2e1 * t477 * t597, 0.2e1 * t565, 0.2e1 * t481 * t521 * t607 + 0.2e1 * t482 * t522 * t603 + 0.2e1 * t483 * t523 * t599, t509 * t521 ^ 2 + t511 * t522 ^ 2 + t513 * t523 ^ 2, 0, 0, -pkin(1) * t565 + t456 * t647 + t458 * t646 + t460 * t645, (t440 * t684 + t448 * t470) * t531 + (t438 * t685 + t446 * t468) * t529 + (t436 * t686 + t444 * t466) * t527 + (t450 * t674 + t452 * t672 + t454 * t670) * pkin(1), 1; t430, 0, 0, t429, t428, t426, t427, t465, 0, 0, t455 * t647 + t457 * t646 + t459 * t645 + (-t478 * t587 - t479 * t584 - t480 * t581) * pkin(1), (t439 * t684 + t447 * t470) * t531 + (t437 * t685 + t445 * t468) * t529 + (t435 * t686 + t443 * t466) * t527 + (t449 * t674 + t451 * t672 + t453 * t670) * pkin(1), 0; t442, 0, 0, t434, t432, t461, t462, 0, 0, 0, -pkin(1) * t461 + t481 * t588 + t482 * t589 + t483 * t590, (t470 * t489 + t483 * t681) * t531 + (t468 * t488 + t482 * t682) * t529 + (t466 * t487 + t481 * t683) * t527 + (-qJ(3,1) * t574 - qJ(3,2) * t576 - qJ(3,3) * t578) * pkin(1), 0; t430, 0, 0, t429, t428, t426, t427, t465, 0, 0, t456 * t650 + t458 * t649 + t460 * t648 + (-t481 * t586 - t482 * t583 - t483 * t580) * pkin(1), (t440 * t687 + t448 * t471) * t531 + (t438 * t688 + t446 * t469) * t529 + (t436 * t689 + t444 * t467) * t527 + (t450 * t673 + t452 * t671 + t454 * t669) * pkin(1), 0; t472 * t679 + t473 * t677 + t474 * t675, 0, 0, t472 * t641 + t473 * t637 + t474 * t633, 0.2e1 * t472 * t605 + 0.2e1 * t473 * t601 + 0.2e1 * t474 * t597, 0.2e1 * t566, 0.2e1 * t478 * t524 * t607 + 0.2e1 * t479 * t525 * t603 + 0.2e1 * t480 * t526 * t599, t509 * t524 ^ 2 + t511 * t525 ^ 2 + t513 * t526 ^ 2, 0, 0, -pkin(1) * t566 + t455 * t650 + t457 * t649 + t459 * t648, (t439 * t687 + t447 * t471) * t531 + (t437 * t688 + t445 * t469) * t529 + (t435 * t689 + t443 * t467) * t527 + (t449 * t673 + t451 * t671 + t453 * t669) * pkin(1), 1; t441, 0, 0, t433, t431, t463, t464, 0, 0, 0, -pkin(1) * t463 + t478 * t588 + t479 * t589 + t480 * t590, (t471 * t489 + t480 * t681) * t531 + (t469 * t488 + t479 * t682) * t529 + (t467 * t487 + t478 * t683) * t527 + (-qJ(3,1) * t573 - qJ(3,2) * t575 - qJ(3,3) * t577) * pkin(1), 0; t442, 0, 0, t434, t432, t461, t462, 0, 0, 0, t456 * t665 + t458 * t663 + t460 * t661, (t440 * t560 + t448 * t495) * t531 + (t438 * t558 + t446 * t494) * t529 + (t436 * t556 + t444 * t493) * t527, 0; t441, 0, 0, t433, t431, t463, t464, 0, 0, 0, t455 * t665 + t457 * t663 + t459 * t661, (t439 * t560 + t447 * t495) * t531 + (t437 * t558 + t445 * t494) * t529 + (t435 * t556 + t443 * t493) * t527, 0; t660 + t662 + t664, 0, 0, t533 * t664 + t534 * t662 + t535 * t660, 0.2e1 * t537 * t629 + 0.2e1 * t539 * t628 + 0.2e1 * t541 * t627, 0, 0, 0, 0, 0, 0.2e1 * qJ(3,1) * t660 + 0.2e1 * qJ(3,2) * t662 + 0.2e1 * qJ(3,3) * t664, (t486 * t560 + t489 * t495) * t531 + (t485 * t558 + t488 * t494) * t529 + (t484 * t556 + t487 * t493) * t527, 1;];
tau_reg  = t1;
