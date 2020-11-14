% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRRR8V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 16:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR8V1G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 16:50:15
% EndTime: 2020-08-06 16:50:18
% DurationCPUTime: 3.55s
% Computational Cost: add. (4670->263), mult. (13044->651), div. (876->14), fcn. (15234->22), ass. (0->324)
t545 = cos(qJ(2,3));
t539 = sin(qJ(2,3));
t544 = cos(qJ(3,3));
t662 = t539 * t544;
t507 = pkin(2) * t662 - t545 * pkin(5);
t532 = sin(pkin(3));
t534 = cos(pkin(3));
t538 = sin(qJ(3,3));
t672 = t534 * t538;
t611 = pkin(2) * t672 + t507 * t532;
t481 = 0.1e1 / t611 ^ 2;
t526 = 0.1e1 / t544 ^ 2;
t745 = t481 * t526;
t547 = cos(qJ(2,2));
t541 = sin(qJ(2,2));
t546 = cos(qJ(3,2));
t658 = t541 * t546;
t508 = pkin(2) * t658 - t547 * pkin(5);
t540 = sin(qJ(3,2));
t670 = t534 * t540;
t610 = pkin(2) * t670 + t508 * t532;
t483 = 0.1e1 / t610 ^ 2;
t528 = 0.1e1 / t546 ^ 2;
t744 = t483 * t528;
t549 = cos(qJ(2,1));
t543 = sin(qJ(2,1));
t548 = cos(qJ(3,1));
t654 = t543 * t548;
t509 = pkin(2) * t654 - t549 * pkin(5);
t542 = sin(qJ(3,1));
t668 = t534 * t542;
t609 = pkin(2) * t668 + t509 * t532;
t485 = 0.1e1 / t609 ^ 2;
t530 = 0.1e1 / t548 ^ 2;
t743 = t485 * t530;
t742 = 0.1e1 / t609;
t741 = 0.1e1 / t610;
t740 = 0.1e1 / t611;
t535 = legFrame(3,3);
t516 = sin(t535);
t519 = cos(t535);
t531 = sin(pkin(6));
t533 = cos(pkin(6));
t489 = -t531 * t516 + t519 * t533;
t671 = t534 * t539;
t496 = t531 * t545 + t533 * t671;
t499 = -t531 * t671 + t533 * t545;
t679 = t532 * t544;
t459 = (-t496 * t519 - t516 * t499) * t538 - t489 * t679;
t492 = t533 * t516 + t519 * t531;
t462 = (-t516 * t496 + t499 * t519) * t538 - t492 * t679;
t739 = t459 * t462;
t536 = legFrame(2,3);
t517 = sin(t536);
t520 = cos(t536);
t490 = -t531 * t517 + t520 * t533;
t669 = t534 * t541;
t497 = t531 * t547 + t533 * t669;
t500 = -t531 * t669 + t533 * t547;
t677 = t532 * t546;
t460 = (-t497 * t520 - t517 * t500) * t540 - t490 * t677;
t493 = t533 * t517 + t520 * t531;
t463 = (-t517 * t497 + t500 * t520) * t540 - t493 * t677;
t738 = t460 * t463;
t537 = legFrame(1,3);
t518 = sin(t537);
t521 = cos(t537);
t491 = -t531 * t518 + t521 * t533;
t667 = t534 * t543;
t495 = t531 * t667 - t533 * t549;
t498 = t531 * t549 + t533 * t667;
t675 = t532 * t548;
t461 = (t518 * t495 - t498 * t521) * t542 - t491 * t675;
t494 = t533 * t518 + t521 * t531;
t464 = (-t495 * t521 - t518 * t498) * t542 - t494 * t675;
t737 = t461 * t464;
t727 = pkin(2) * t544;
t510 = pkin(5) * t539 + t545 * t727;
t728 = pkin(2) * t532;
t560 = -t507 * t534 + t538 * t728;
t465 = t510 * t533 + t560 * t531;
t468 = t531 * t510 - t560 * t533;
t444 = t465 * t516 + t468 * t519;
t718 = t444 * t459;
t441 = t465 * t519 - t516 * t468;
t724 = t441 * t462;
t736 = t481 * (t718 + t724);
t726 = pkin(2) * t546;
t511 = pkin(5) * t541 + t547 * t726;
t559 = -t508 * t534 + t540 * t728;
t466 = t511 * t533 + t559 * t531;
t469 = t531 * t511 - t559 * t533;
t445 = t466 * t517 + t469 * t520;
t716 = t445 * t460;
t442 = t466 * t520 - t517 * t469;
t722 = t442 * t463;
t735 = t483 * (t716 + t722);
t725 = pkin(2) * t548;
t512 = pkin(5) * t543 + t549 * t725;
t558 = -t509 * t534 + t542 * t728;
t467 = t512 * t533 + t558 * t531;
t470 = t531 * t512 - t558 * t533;
t446 = t467 * t518 + t470 * t521;
t714 = t446 * t461;
t443 = t467 * t521 - t518 * t470;
t720 = t443 * t464;
t734 = t485 * (t714 + t720);
t733 = t526 * t538;
t732 = t528 * t540;
t731 = t530 * t542;
t730 = 0.2e1 * t532;
t550 = 0.1e1 / pkin(2);
t729 = 0.2e1 * t550;
t723 = t441 * t740;
t721 = t442 * t741;
t719 = t443 * t742;
t717 = t444 * t740;
t715 = t445 * t741;
t713 = t446 * t742;
t712 = t740 ^ 2;
t525 = 0.1e1 / t544;
t711 = t740 * t525;
t710 = t740 * t539;
t708 = t741 ^ 2;
t527 = 0.1e1 / t546;
t707 = t741 * t527;
t706 = t741 * t541;
t704 = t742 ^ 2;
t529 = 0.1e1 / t548;
t703 = t742 * t529;
t702 = t742 * t543;
t700 = t740 * t545;
t699 = t481 * t525;
t697 = t481 * t545;
t696 = t741 * t547;
t695 = t483 * t527;
t693 = t483 * t547;
t692 = t742 * t549;
t691 = t485 * t529;
t689 = t485 * t549;
t688 = t525 * t538;
t687 = t525 * t539;
t686 = t525 * t545;
t685 = t527 * t540;
t684 = t527 * t541;
t683 = t527 * t547;
t682 = t529 * t542;
t681 = t529 * t543;
t680 = t529 * t549;
t678 = t532 * t545;
t676 = t532 * t547;
t674 = t532 * t549;
t673 = t532 * t550;
t666 = t534 * t550;
t665 = t538 * t539;
t664 = t539 * t489;
t663 = t539 * t492;
t661 = t540 * t541;
t660 = t541 * t490;
t659 = t541 * t493;
t657 = t542 * t543;
t656 = t543 * t491;
t655 = t543 * t494;
t653 = t545 * t489;
t652 = t545 * t492;
t651 = t547 * t490;
t650 = t547 * t493;
t649 = t549 * t491;
t648 = t549 * t494;
t647 = t740 * t711;
t501 = -t532 * t665 + t534 * t544;
t646 = t501 * t711;
t503 = t532 * t662 + t672;
t645 = t503 * t711;
t644 = t550 * t711;
t643 = t740 * t666;
t642 = t741 * t707;
t502 = -t532 * t661 + t534 * t546;
t641 = t502 * t707;
t504 = t532 * t658 + t670;
t640 = t504 * t707;
t639 = t550 * t707;
t638 = t741 * t666;
t637 = t742 * t703;
t505 = -t532 * t657 + t534 * t548;
t636 = t505 * t703;
t506 = t532 * t654 + t668;
t635 = t506 * t703;
t634 = t550 * t703;
t633 = t742 * t666;
t632 = t740 * t687;
t631 = t740 * t686;
t630 = t740 * t678;
t522 = t538 ^ 2;
t629 = t522 * t745;
t628 = t481 * t688;
t627 = t481 * t678;
t626 = t741 * t684;
t625 = t741 * t683;
t624 = t741 * t676;
t523 = t540 ^ 2;
t623 = t523 * t744;
t622 = t483 * t685;
t621 = t483 * t676;
t620 = t742 * t681;
t619 = t742 * t680;
t618 = t742 * t674;
t524 = t542 ^ 2;
t617 = t524 * t743;
t616 = t485 * t682;
t615 = t485 * t674;
t614 = t538 * t686;
t613 = t540 * t683;
t612 = t542 * t680;
t608 = t441 * t459 * t699;
t607 = t442 * t460 * t695;
t606 = t443 * t461 * t691;
t605 = t444 * t462 * t699;
t604 = t445 * t463 * t695;
t603 = t446 * t464 * t691;
t447 = -(t534 * t653 - t663) * t727 - pkin(5) * (t534 * t664 + t652);
t602 = t447 * t647;
t448 = -(t534 * t652 + t664) * t727 - (t534 * t663 - t653) * pkin(5);
t601 = t448 * t647;
t449 = -(t534 * t651 - t659) * t726 - pkin(5) * (t534 * t660 + t650);
t600 = t449 * t642;
t450 = -(t534 * t650 + t660) * t726 - (t534 * t659 - t651) * pkin(5);
t599 = t450 * t642;
t451 = -(t534 * t649 - t655) * t725 - pkin(5) * (t534 * t656 + t648);
t598 = t451 * t637;
t452 = -(t534 * t648 + t656) * t725 - (t534 * t655 - t649) * pkin(5);
t597 = t452 * t637;
t596 = t745 * t739;
t595 = t744 * t738;
t594 = t743 * t737;
t593 = t712 * t733;
t592 = t501 * t644;
t591 = t503 * t644;
t590 = t665 * t711;
t589 = t673 * t710;
t588 = t708 * t732;
t587 = t502 * t639;
t586 = t504 * t639;
t585 = t661 * t707;
t584 = t673 * t706;
t583 = t704 * t731;
t582 = t505 * t634;
t581 = t506 * t634;
t580 = t657 * t703;
t579 = t673 * t702;
t578 = t740 * t614;
t577 = t481 * t614;
t576 = t741 * t613;
t575 = t483 * t613;
t574 = t742 * t612;
t573 = t485 * t612;
t572 = t441 * t601;
t571 = t442 * t599;
t570 = t443 * t597;
t569 = t444 * t602;
t568 = t445 * t600;
t567 = t446 * t598;
t566 = t532 * t577;
t565 = t532 * t575;
t564 = t532 * t573;
t557 = t589 * t688;
t556 = t584 * t685;
t555 = t579 * t682;
t554 = (t447 * t462 + t448 * t459) * t712;
t553 = (t449 * t463 + t450 * t460) * t708;
t552 = (t451 * t464 + t452 * t461) * t704;
t551 = 0.1e1 / pkin(2) ^ 2;
t458 = t464 ^ 2;
t457 = t463 ^ 2;
t456 = t462 ^ 2;
t455 = t461 ^ 2;
t454 = t460 ^ 2;
t453 = t459 ^ 2;
t440 = (t452 * t633 + t464 * t618) * t529;
t439 = (t451 * t633 + t461 * t618) * t529;
t438 = (t450 * t638 + t463 * t624) * t527;
t437 = (t449 * t638 + t460 * t624) * t527;
t436 = (t448 * t643 + t462 * t630) * t525;
t435 = (t447 * t643 + t459 * t630) * t525;
t434 = t440 * t542;
t433 = t440 * t548;
t432 = t439 * t542;
t431 = t439 * t548;
t430 = t438 * t540;
t429 = t438 * t546;
t428 = t437 * t540;
t427 = t437 * t546;
t426 = t436 * t538;
t425 = t436 * t544;
t424 = t435 * t538;
t423 = t435 * t544;
t422 = (-t462 * t632 - t463 * t626 - t464 * t620) * t532;
t421 = (-t459 * t632 - t460 * t626 - t461 * t620) * t532;
t420 = (t462 * t631 + t463 * t625 + t464 * t619) * t532;
t419 = (t459 * t631 + t460 * t625 + t461 * t619) * t532;
t418 = -t452 * t579 - t434;
t417 = -t451 * t579 - t432;
t416 = -t450 * t584 - t430;
t415 = -t449 * t584 - t428;
t414 = -t448 * t589 - t426;
t413 = -t447 * t589 - t424;
t412 = -t452 * t555 + t433;
t411 = -t451 * t555 + t431;
t410 = -t450 * t556 + t429;
t409 = -t449 * t556 + t427;
t408 = -t448 * t557 + t425;
t407 = -t447 * t557 + t423;
t406 = t713 + t715 + t717;
t405 = t719 + t721 + t723;
t404 = t594 + t595 + t596;
t403 = t522 * t596 + t523 * t595 + t524 * t594;
t402 = 0.2e1 * t616 * t737 + 0.2e1 * t622 * t738 + 0.2e1 * t628 * t739;
t401 = (t447 * t448 * t745 + t449 * t450 * t744 + t451 * t452 * t743) * t551;
t400 = t441 * t481 * t444 + t442 * t483 * t445 + t443 * t485 * t446;
t399 = (t680 * t734 + t683 * t735 + t686 * t736) * t532;
t398 = (-t681 * t734 - t684 * t735 - t687 * t736) * t532;
t397 = (t525 * t554 + t527 * t553 + t529 * t552) * t550;
t396 = (t552 * t731 + t553 * t732 + t554 * t733) * t550;
t1 = [t441 ^ 2 * t481 + t442 ^ 2 * t483 + t443 ^ 2 * t485, t453 * t745 + t454 * t744 + t455 * t743, (t545 * t608 + t547 * t607 + t549 * t606) * t730, (-t539 * t608 - t541 * t607 - t543 * t606) * t730, t453 * t629 + t454 * t623 + t455 * t617, 0.2e1 * t453 * t628 + 0.2e1 * t454 * t622 + 0.2e1 * t455 * t616, (t447 * t459 * t593 + t449 * t460 * t588 + t451 * t461 * t583) * t729, (t459 * t602 + t460 * t600 + t461 * t598) * t729, (t447 ^ 2 * t745 + t449 ^ 2 * t744 + t451 ^ 2 * t743) * t551, (t461 * t615 + (t451 * t582 + t411) * t742) * t443 + (t460 * t621 + (t449 * t587 + t409) * t741) * t442 + (t459 * t627 + (t447 * t592 + t407) * t740) * t441, (-t461 * t564 + (-t451 * t581 + t417) * t742) * t443 + (-t460 * t565 + (-t449 * t586 + t415) * t741) * t442 + (-t459 * t566 + (-t447 * t591 + t413) * t740) * t441, 1; t400, t404, t399, t398, t403, t402, t396, t397, t401, t408 * t723 + t410 * t721 + t412 * t719 + (t501 * t569 + t502 * t568 + t505 * t567) * t550 + (t689 * t714 + t693 * t716 + t697 * t718) * t532, t414 * t723 + t416 * t721 + t418 * t719 + (-t503 * t569 - t504 * t568 - t506 * t567) * t550 + (-t573 * t714 - t575 * t716 - t577 * t718) * t532, 0; t405, 0, t419, t421, 0, 0, 0, 0, 0, (t447 * t646 + t449 * t641 + t451 * t636) * t550 + (t459 * t700 + t460 * t696 + t461 * t692) * t532, (-t447 * t645 - t449 * t640 - t451 * t635) * t550 + (-t459 * t578 - t460 * t576 - t461 * t574) * t532, 0; t400, t404, t399, t398, t403, t402, t396, t397, t401, t407 * t717 + t409 * t715 + t411 * t713 + (t501 * t572 + t502 * t571 + t505 * t570) * t550 + (t689 * t720 + t693 * t722 + t697 * t724) * t532, t413 * t717 + t415 * t715 + t417 * t713 + (-t503 * t572 - t504 * t571 - t506 * t570) * t550 + (-t573 * t720 - t575 * t722 - t577 * t724) * t532, 0; t444 ^ 2 * t481 + t445 ^ 2 * t483 + t446 ^ 2 * t485, t456 * t745 + t457 * t744 + t458 * t743, (t545 * t605 + t547 * t604 + t549 * t603) * t730, (-t539 * t605 - t541 * t604 - t543 * t603) * t730, t456 * t629 + t457 * t623 + t458 * t617, 0.2e1 * t456 * t628 + 0.2e1 * t457 * t622 + 0.2e1 * t458 * t616, (t448 * t462 * t593 + t450 * t463 * t588 + t452 * t464 * t583) * t729, (t462 * t601 + t463 * t599 + t464 * t597) * t729, (t448 ^ 2 * t745 + t450 ^ 2 * t744 + t452 ^ 2 * t743) * t551, (t464 * t615 + (t452 * t582 + t412) * t742) * t446 + (t463 * t621 + (t450 * t587 + t410) * t741) * t445 + (t462 * t627 + (t448 * t592 + t408) * t740) * t444, (-t464 * t564 + (-t452 * t581 + t418) * t742) * t446 + (-t463 * t565 + (-t450 * t586 + t416) * t741) * t445 + (-t462 * t566 + (-t448 * t591 + t414) * t740) * t444, 1; t406, 0, t420, t422, 0, 0, 0, 0, 0, (t448 * t646 + t450 * t641 + t452 * t636) * t550 + (t462 * t700 + t463 * t696 + t464 * t692) * t532, (-t448 * t645 - t450 * t640 - t452 * t635) * t550 + (-t462 * t578 - t463 * t576 - t464 * t574) * t532, 0; t405, 0, t419, t421, 0, 0, 0, 0, 0, t423 + t427 + t431 + (-t447 * t590 - t449 * t585 - t451 * t580) * t673, -t424 - t428 - t432 + (-t447 * t710 - t449 * t706 - t451 * t702) * t673, 0; t406, 0, t420, t422, 0, 0, 0, 0, 0, t425 + t429 + t433 + (-t448 * t590 - t450 * t585 - t452 * t580) * t673, -t426 - t430 - t434 + (-t448 * t710 - t450 * t706 - t452 * t702) * t673, 0; 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
tau_reg  = t1;
