% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P4PRRR1G1P1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4*4x11]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRR1G1P1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_regmin: pkin has to be [7x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:14:37
% EndTime: 2020-03-02 20:14:44
% DurationCPUTime: 7.31s
% Computational Cost: add. (11876->320), mult. (7348->580), div. (1668->11), fcn. (9228->34), ass. (0->282)
t531 = pkin(7) + qJ(2,1);
t588 = qJ(3,1) + t531;
t513 = sin(t588);
t709 = legFrame(1,3);
t525 = sin(t709);
t560 = cos(t588);
t592 = cos(t709);
t552 = -t513 * t525 + t592 * t560;
t518 = sin(t531);
t521 = cos(t531);
t484 = t513 * t521 - t518 * t560;
t710 = 0.1e1 / t484;
t657 = t710 * t552;
t530 = pkin(7) + qJ(2,2);
t587 = qJ(3,2) + t530;
t512 = sin(t587);
t708 = legFrame(2,3);
t524 = sin(t708);
t559 = cos(t587);
t591 = cos(t708);
t551 = -t512 * t524 + t591 * t559;
t517 = sin(t530);
t520 = cos(t530);
t483 = t512 * t520 - t517 * t559;
t711 = 0.1e1 / t483;
t661 = t711 * t551;
t529 = pkin(7) + qJ(2,3);
t586 = qJ(3,3) + t529;
t511 = sin(t586);
t707 = legFrame(3,3);
t523 = sin(t707);
t558 = cos(t586);
t590 = cos(t707);
t550 = -t511 * t523 + t590 * t558;
t516 = sin(t529);
t519 = cos(t529);
t482 = t511 * t519 - t516 * t558;
t712 = 0.1e1 / t482;
t665 = t712 * t550;
t528 = pkin(7) + qJ(2,4);
t585 = qJ(3,4) + t528;
t506 = sin(t585);
t706 = legFrame(4,3);
t522 = sin(t706);
t557 = cos(t585);
t589 = cos(t706);
t549 = -t506 * t522 + t589 * t557;
t514 = sin(t528);
t515 = cos(t528);
t472 = t506 * t515 - t514 * t557;
t713 = 0.1e1 / t472;
t675 = t713 * t549;
t540 = xP(4);
t526 = sin(t540);
t527 = cos(t540);
t541 = koppelP(4,1);
t714 = koppelP(4,2);
t493 = t526 * t541 + t527 * t714;
t497 = -t526 * t714 + t527 * t541;
t461 = t493 * t589 - t522 * t497;
t465 = t522 * t493 + t497 * t589;
t721 = t461 * t557 - t465 * t506;
t544 = koppelP(1,1);
t717 = koppelP(1,2);
t496 = t526 * t544 + t527 * t717;
t500 = -t526 * t717 + t527 * t544;
t464 = t496 * t592 - t525 * t500;
t468 = t525 * t496 + t500 * t592;
t720 = t464 * t560 - t468 * t513;
t543 = koppelP(2,1);
t716 = koppelP(2,2);
t495 = t526 * t543 + t527 * t716;
t499 = -t526 * t716 + t527 * t543;
t463 = t495 * t591 - t524 * t499;
t467 = t524 * t495 + t499 * t591;
t719 = t463 * t559 - t467 * t512;
t542 = koppelP(3,1);
t715 = koppelP(3,2);
t494 = t526 * t542 + t527 * t715;
t498 = -t526 * t715 + t527 * t542;
t462 = t494 * t590 - t523 * t498;
t466 = t523 * t494 + t498 * t590;
t718 = t462 * t558 - t466 * t511;
t401 = pkin(2) * (t461 * t515 - t465 * t514) + t721 * pkin(3);
t705 = t401 * t713;
t470 = 0.1e1 / t472 ^ 2;
t704 = t401 * t470;
t402 = pkin(2) * (t462 * t519 - t466 * t516) + t718 * pkin(3);
t703 = t402 * t712;
t474 = 0.1e1 / t482 ^ 2;
t702 = t402 * t474;
t403 = pkin(2) * (t463 * t520 - t467 * t517) + t719 * pkin(3);
t701 = t403 * t711;
t476 = 0.1e1 / t483 ^ 2;
t700 = t403 * t476;
t404 = pkin(2) * (t464 * t521 - t468 * t518) + t720 * pkin(3);
t699 = t404 * t710;
t478 = 0.1e1 / t484 ^ 2;
t698 = t404 * t478;
t697 = t721 * t713;
t696 = t721 * t470;
t695 = t718 * t712;
t694 = t718 * t474;
t693 = t719 * t711;
t692 = t719 * t476;
t691 = t720 * t710;
t690 = t720 * t478;
t486 = t589 * t506 + t522 * t557;
t453 = pkin(2) * (t589 * t514 + t522 * t515) + t486 * pkin(3);
t689 = t453 * t713;
t454 = -pkin(2) * (t514 * t522 - t589 * t515) + t549 * pkin(3);
t688 = t454 * t713;
t490 = t590 * t511 + t523 * t558;
t455 = pkin(2) * (t590 * t516 + t523 * t519) + t490 * pkin(3);
t687 = t455 * t712;
t491 = t591 * t512 + t524 * t559;
t456 = pkin(2) * (t591 * t517 + t524 * t520) + t491 * pkin(3);
t686 = t456 * t711;
t492 = t592 * t513 + t525 * t560;
t457 = pkin(2) * (t592 * t518 + t525 * t521) + t492 * pkin(3);
t685 = t457 * t710;
t458 = -pkin(2) * (t516 * t523 - t590 * t519) + t550 * pkin(3);
t684 = t458 * t712;
t459 = -pkin(2) * (t517 * t524 - t591 * t520) + t551 * pkin(3);
t683 = t459 * t711;
t460 = -pkin(2) * (t518 * t525 - t592 * t521) + t552 * pkin(3);
t682 = t460 * t710;
t677 = t713 ^ 2;
t546 = 0.1e1 / pkin(2);
t676 = t713 * t546;
t674 = t713 * t486;
t532 = sin(qJ(3,4));
t673 = t713 * t532;
t533 = cos(qJ(3,4));
t672 = t713 * t533;
t671 = t712 ^ 2;
t670 = t712 * t546;
t669 = t711 ^ 2;
t668 = t711 * t546;
t667 = t710 ^ 2;
t666 = t710 * t546;
t664 = t712 * t490;
t534 = sin(qJ(3,3));
t663 = t712 * t534;
t537 = cos(qJ(3,3));
t662 = t712 * t537;
t660 = t711 * t491;
t535 = sin(qJ(3,2));
t659 = t711 * t535;
t538 = cos(qJ(3,2));
t658 = t711 * t538;
t656 = t710 * t492;
t536 = sin(qJ(3,1));
t655 = t710 * t536;
t539 = cos(qJ(3,1));
t654 = t710 * t539;
t545 = 0.1e1 / pkin(3);
t653 = t545 * t546;
t616 = t713 * t653;
t583 = t454 * t616;
t584 = t453 * t616;
t652 = t493 * t583 - t497 * t584;
t611 = t712 * t653;
t579 = t458 * t611;
t582 = t455 * t611;
t651 = t494 * t579 - t498 * t582;
t606 = t711 * t653;
t578 = t459 * t606;
t581 = t456 * t606;
t650 = t495 * t578 - t499 * t581;
t601 = t710 * t653;
t577 = t460 * t601;
t580 = t457 * t601;
t649 = t496 * t577 - t500 * t580;
t600 = t549 * t676;
t445 = t493 * t600;
t599 = t486 * t676;
t446 = t497 * t599;
t405 = -t445 + t446;
t598 = t550 * t670;
t447 = t494 * t598;
t595 = t490 * t670;
t450 = t498 * t595;
t410 = -t447 + t450;
t597 = t551 * t668;
t448 = t495 * t597;
t594 = t491 * t668;
t451 = t499 * t594;
t411 = -t448 + t451;
t596 = t552 * t666;
t449 = t496 * t596;
t593 = t492 * t666;
t452 = t500 * t593;
t412 = -t449 + t452;
t648 = t405 * t705;
t647 = t532 * t704;
t646 = t533 * t704;
t645 = t410 * t703;
t644 = t534 * t702;
t643 = t537 * t702;
t642 = t411 * t701;
t641 = t535 * t700;
t640 = t538 * t700;
t639 = t412 * t699;
t638 = t536 * t698;
t637 = t539 * t698;
t636 = t405 * t673;
t635 = t405 * t672;
t634 = t410 * t663;
t633 = t410 * t662;
t632 = t411 * t659;
t631 = t411 * t658;
t630 = t412 * t655;
t629 = t412 * t654;
t628 = t532 * t697;
t627 = t533 * t697;
t626 = t534 * t695;
t625 = t537 * t695;
t624 = t535 * t693;
t623 = t538 * t693;
t622 = t536 * t691;
t621 = t539 * t691;
t620 = t549 * t673;
t619 = t549 * t672;
t618 = t486 * t673;
t617 = t486 * t672;
t615 = t550 * t663;
t614 = t550 * t662;
t613 = t490 * t663;
t612 = t490 * t662;
t610 = t551 * t659;
t609 = t551 * t658;
t608 = t491 * t659;
t607 = t491 * t658;
t605 = t552 * t655;
t604 = t552 * t654;
t603 = t492 * t655;
t602 = t492 * t654;
t576 = t713 * t620;
t575 = t713 * t619;
t574 = t713 * t618;
t573 = t713 * t617;
t572 = t712 * t615;
t571 = t712 * t614;
t570 = t712 * t613;
t569 = t712 * t612;
t568 = t711 * t610;
t567 = t711 * t609;
t566 = t711 * t608;
t565 = t711 * t607;
t564 = t710 * t605;
t563 = t710 * t604;
t562 = t710 * t603;
t561 = t710 * t602;
t547 = 0.1e1 / pkin(2) ^ 2;
t548 = (t656 * t657 + t660 * t661 + t664 * t665 + t674 * t675) * t547;
t505 = t526 ^ 2 + t527 ^ 2;
t424 = -t577 + t596;
t423 = -t580 + t593;
t422 = -t578 + t597;
t421 = -t581 + t594;
t420 = -t579 + t598;
t419 = -t582 + t595;
t418 = -t577 + 0.2e1 * t596;
t417 = -t578 + 0.2e1 * t597;
t416 = -t579 + 0.2e1 * t598;
t415 = -t580 + 0.2e1 * t593;
t414 = -t581 + 0.2e1 * t594;
t413 = -t582 + 0.2e1 * t595;
t409 = -t583 + t600;
t408 = -t584 + t599;
t407 = -t583 + 0.2e1 * t600;
t406 = -t584 + 0.2e1 * t599;
t400 = t412 + t649;
t399 = t411 + t650;
t398 = t410 + t651;
t397 = -0.2e1 * t449 + 0.2e1 * t452 + t649;
t396 = -0.2e1 * t448 + 0.2e1 * t451 + t650;
t395 = -0.2e1 * t447 + 0.2e1 * t450 + t651;
t394 = t405 + t652;
t393 = -0.2e1 * t445 + 0.2e1 * t446 + t652;
t1 = [0, (t549 ^ 2 * t677 + t550 ^ 2 * t671 + t551 ^ 2 * t669 + t552 ^ 2 * t667) * t547, 0, 0, (t409 * t675 + t420 * t665 + t422 * t661 + t424 * t657 + (-t409 * t688 - t420 * t684 - t422 * t683 - t424 * t682) * t545) * t546, t407 * t619 + t416 * t614 + t417 * t609 + t418 * t604 + (-t454 * t575 - t458 * t571 - t459 * t567 - t460 * t563) * t653, -t407 * t620 - t416 * t615 - t417 * t610 - t418 * t605 + (t454 * t576 + t458 * t572 + t459 * t568 + t460 * t564) * t653, 0, 0, 0, t505; 0, t548, 0, 0, (t408 * t675 + t419 * t665 + t421 * t661 + t423 * t657 + (-t408 * t688 - t419 * t684 - t421 * t683 - t423 * t682) * t545) * t546, t406 * t619 + t413 * t614 + t414 * t609 + t415 * t604 + (-t454 * t573 - t458 * t569 - t459 * t565 - t460 * t561) * t653, -t406 * t620 - t413 * t615 - t414 * t610 - t415 * t605 + (t454 * t574 + t458 * t570 + t459 * t566 + t460 * t562) * t653, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (t405 * t675 + t410 * t665 + t411 * t661 + t412 * t657) * t546, 0, 0, (t394 * t675 + t398 * t665 + t399 * t661 + t400 * t657 + (-t394 * t688 - t398 * t684 - t399 * t683 - t400 * t682) * t545) * t546, t393 * t619 + t395 * t614 + t396 * t609 + t397 * t604 + (-t454 * t635 - t458 * t633 - t459 * t631 - t460 * t629) * t545, -t393 * t620 - t395 * t615 - t396 * t610 - t397 * t605 + (t454 * t636 + t458 * t634 + t459 * t632 + t460 * t630) * t545, 0, -t526, -t527, 0; 0, t548, 0, 0, (t409 * t674 + t420 * t664 + t422 * t660 + t424 * t656 + (-t409 * t689 - t420 * t687 - t422 * t686 - t424 * t685) * t545) * t546, t407 * t617 + t416 * t612 + t417 * t607 + t418 * t602 + (-t453 * t575 - t455 * t571 - t456 * t567 - t457 * t563) * t653, -t407 * t618 - t416 * t613 - t417 * t608 - t418 * t603 + (t453 * t576 + t455 * t572 + t456 * t568 + t457 * t564) * t653, 0, 0, 0, 0; 0, (t486 ^ 2 * t677 + t490 ^ 2 * t671 + t491 ^ 2 * t669 + t492 ^ 2 * t667) * t547, 0, 0, (t408 * t674 + t419 * t664 + t421 * t660 + t423 * t656 + (-t408 * t689 - t419 * t687 - t421 * t686 - t423 * t685) * t545) * t546, t406 * t617 + t413 * t612 + t414 * t607 + t415 * t602 + (-t453 * t573 - t455 * t569 - t456 * t565 - t457 * t561) * t653, -t406 * t618 - t413 * t613 - t414 * t608 - t415 * t603 + (t453 * t574 + t455 * t570 + t456 * t566 + t457 * t562) * t653, 0, 0, 0, t505; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (t405 * t674 + t410 * t664 + t411 * t660 + t412 * t656) * t546, 0, 0, (t394 * t674 + t398 * t664 + t399 * t660 + t400 * t656 + (-t394 * t689 - t398 * t687 - t399 * t686 - t400 * t685) * t545) * t546, t393 * t617 + t395 * t612 + t396 * t607 + t397 * t602 + (-t453 * t635 - t455 * t633 - t456 * t631 - t457 * t629) * t545, -t393 * t618 - t395 * t613 - t396 * t608 - t397 * t603 + (t453 * t636 + t455 * t634 + t456 * t632 + t457 * t630) * t545, 0, t527, -t526, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (-t549 * t696 - t550 * t694 - t551 * t692 - t552 * t690) * t547, 0, 0, (-t409 * t697 - t420 * t695 - t422 * t693 - t424 * t691 + (t409 * t705 + t420 * t703 + t422 * t701 + t424 * t699) * t545) * t546, -t407 * t627 - t416 * t625 - t417 * t623 - t418 * t621 + (t549 * t646 + t550 * t643 + t551 * t640 + t552 * t637) * t653, t407 * t628 + t416 * t626 + t417 * t624 + t418 * t622 + (-t549 * t647 - t550 * t644 - t551 * t641 - t552 * t638) * t653, 0, -t526, -t527, 0; 0, (-t486 * t696 - t490 * t694 - t491 * t692 - t492 * t690) * t547, 0, 0, (-t408 * t697 - t419 * t695 - t421 * t693 - t423 * t691 + (t408 * t705 + t419 * t703 + t421 * t701 + t423 * t699) * t545) * t546, -t406 * t627 - t413 * t625 - t414 * t623 - t415 * t621 + (t486 * t646 + t490 * t643 + t491 * t640 + t492 * t637) * t653, t406 * t628 + t413 * t626 + t414 * t624 + t415 * t622 + (-t486 * t647 - t490 * t644 - t491 * t641 - t492 * t638) * t653, 0, t527, -t526, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, (-t405 * t697 - t410 * t695 - t411 * t693 - t412 * t691) * t546, 0, 0, (-t394 * t697 - t398 * t695 - t399 * t693 - t400 * t691 + (t394 * t705 + t398 * t703 + t399 * t701 + t400 * t699) * t545) * t546, -t393 * t627 - t395 * t625 - t396 * t623 - t397 * t621 + (t533 * t648 + t537 * t645 + t538 * t642 + t539 * t639) * t545, t393 * t628 + t395 * t626 + t396 * t624 + t397 * t622 + (-t532 * t648 - t534 * t645 - t535 * t642 - t536 * t639) * t545, 1, 0, 0, 0;];
tau_reg  = t1;
