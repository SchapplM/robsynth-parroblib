% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P4PRRRR1G3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [4*4x15]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P4PRRRR1G3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_regmin: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_inertia_para_pf_regmin: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_regmin: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_regmin: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:05:03
% EndTime: 2020-03-02 19:05:23
% DurationCPUTime: 21.52s
% Computational Cost: add. (2573->528), mult. (6740->1016), div. (3196->22), fcn. (8620->26), ass. (0->433)
t454 = legFrame(4,2);
t408 = cos(t454);
t455 = legFrame(3,2);
t409 = cos(t455);
t456 = legFrame(2,2);
t410 = cos(t456);
t457 = legFrame(1,2);
t411 = cos(t457);
t469 = cos(qJ(2,1));
t463 = sin(qJ(2,1));
t434 = 0.1e1 / t463 ^ 2;
t462 = sin(qJ(3,1));
t432 = t462 ^ 2;
t468 = cos(qJ(3,1));
t498 = t468 ^ 2;
t446 = 0.1e1 / t498;
t715 = t446 * t432;
t646 = t434 * t715;
t579 = t469 * t646;
t467 = cos(qJ(2,2));
t461 = sin(qJ(2,2));
t430 = 0.1e1 / t461 ^ 2;
t460 = sin(qJ(3,2));
t428 = t460 ^ 2;
t466 = cos(qJ(3,2));
t495 = t466 ^ 2;
t441 = 0.1e1 / t495;
t719 = t441 * t428;
t654 = t430 * t719;
t583 = t467 * t654;
t465 = cos(qJ(2,3));
t459 = sin(qJ(2,3));
t426 = 0.1e1 / t459 ^ 2;
t458 = sin(qJ(3,3));
t424 = t458 ^ 2;
t464 = cos(qJ(3,3));
t492 = t464 ^ 2;
t436 = 0.1e1 / t492;
t723 = t436 * t424;
t662 = t426 * t723;
t587 = t465 * t662;
t453 = cos(qJ(2,4));
t451 = sin(qJ(2,4));
t417 = 0.1e1 / t451 ^ 2;
t450 = sin(qJ(3,4));
t415 = t450 ^ 2;
t452 = cos(qJ(3,4));
t483 = t452 ^ 2;
t419 = 0.1e1 / t483;
t745 = t419 * t415;
t672 = t417 * t745;
t591 = t453 * t672;
t794 = t408 * t591 + t409 * t587 + t410 * t583 + t411 * t579;
t470 = xP(4);
t412 = sin(t470);
t413 = cos(t470);
t471 = koppelP(4,2);
t475 = koppelP(4,1);
t387 = t412 * t475 + t413 * t471;
t391 = -t412 * t471 + t413 * t475;
t404 = sin(t454);
t355 = t387 * t408 + t391 * t404;
t793 = t355 * t404;
t472 = koppelP(3,2);
t476 = koppelP(3,1);
t388 = t412 * t476 + t413 * t472;
t392 = -t412 * t472 + t413 * t476;
t405 = sin(t455);
t356 = t388 * t409 + t392 * t405;
t792 = t356 * t405;
t473 = koppelP(2,2);
t477 = koppelP(2,1);
t389 = t412 * t477 + t413 * t473;
t393 = -t412 * t473 + t413 * t477;
t406 = sin(t456);
t357 = t389 * t410 + t393 * t406;
t791 = t357 * t406;
t474 = koppelP(1,2);
t478 = koppelP(1,1);
t390 = t412 * t478 + t413 * t474;
t394 = -t412 * t474 + t413 * t478;
t407 = sin(t457);
t358 = t390 * t411 + t394 * t407;
t790 = t358 * t407;
t753 = t411 * t469;
t386 = t407 * t463 + t753;
t433 = 0.1e1 / t463;
t762 = t386 * t433;
t757 = t407 * t469;
t385 = t411 * t463 - t757;
t763 = t385 * t433;
t567 = -t390 * t762 + t394 * t763;
t754 = t410 * t467;
t384 = t406 * t461 + t754;
t429 = 0.1e1 / t461;
t765 = t384 * t429;
t758 = t406 * t467;
t383 = t410 * t461 - t758;
t766 = t383 * t429;
t568 = -t389 * t765 + t393 * t766;
t755 = t409 * t465;
t382 = t405 * t459 + t755;
t425 = 0.1e1 / t459;
t768 = t382 * t425;
t759 = t405 * t465;
t381 = t409 * t459 - t759;
t769 = t381 * t425;
t569 = -t388 * t768 + t392 * t769;
t756 = t408 * t453;
t380 = t404 * t451 + t756;
t416 = 0.1e1 / t451;
t771 = t380 * t416;
t760 = t404 * t453;
t379 = t408 * t451 - t760;
t772 = t379 * t416;
t570 = -t387 * t771 + t391 * t772;
t445 = 0.1e1 / t468;
t716 = t445 * t462;
t644 = t433 * t716;
t440 = 0.1e1 / t466;
t720 = t440 * t460;
t652 = t429 * t720;
t435 = 0.1e1 / t464;
t724 = t435 * t458;
t660 = t425 * t724;
t418 = 0.1e1 / t452;
t746 = t418 * t450;
t670 = t416 * t746;
t789 = t567 * t644 + t568 * t652 + t569 * t660 + t570 * t670;
t788 = t567 * t762 + t568 * t765 + t569 * t768 + t570 * t771;
t787 = t567 * t763 + t568 * t766 + t569 * t769 + t570 * t772;
t479 = 1 / pkin(2);
t786 = 2 * t479;
t480 = 1 / (pkin(2) ^ 2);
t785 = 2 * t480;
t784 = t570 * t418;
t783 = t570 * t416;
t782 = t570 * t453;
t781 = t569 * t435;
t780 = t568 * t440;
t779 = t567 * t445;
t778 = t569 * t425;
t777 = t569 * t465;
t776 = t568 * t429;
t775 = t568 * t467;
t774 = t567 * t433;
t773 = t567 * t469;
t770 = t380 * t417;
t767 = t382 * t426;
t764 = t384 * t430;
t761 = t386 * t434;
t752 = t416 * t418;
t422 = t453 ^ 2;
t751 = t416 * t422;
t750 = t416 * t453;
t749 = t417 * t419;
t748 = t417 * t422;
t747 = t417 * t453;
t744 = t419 * t450;
t420 = t418 * t419;
t743 = t420 * t453;
t742 = t425 * t435;
t439 = t465 ^ 2;
t741 = t425 * t439;
t740 = t425 * t465;
t739 = t426 * t436;
t738 = t426 * t439;
t737 = t426 * t465;
t736 = t429 * t440;
t444 = t467 ^ 2;
t735 = t429 * t444;
t734 = t429 * t467;
t733 = t430 * t441;
t732 = t430 * t444;
t731 = t430 * t467;
t730 = t433 * t445;
t449 = t469 ^ 2;
t729 = t433 * t449;
t728 = t433 * t469;
t727 = t434 * t446;
t726 = t434 * t449;
t725 = t434 * t469;
t722 = t436 * t458;
t437 = t435 * t436;
t721 = t437 * t465;
t718 = t441 * t460;
t442 = t440 * t441;
t717 = t442 * t467;
t714 = t446 * t462;
t447 = t445 * t446;
t713 = t447 * t469;
t712 = t450 * t453;
t711 = t458 * t465;
t710 = t460 * t467;
t709 = t462 * t469;
t708 = t570 * t750;
t688 = t355 * t752;
t349 = t479 * t688;
t707 = t349 * t782;
t706 = t569 * t740;
t705 = t568 * t734;
t704 = t567 * t728;
t687 = t356 * t742;
t350 = t479 * t687;
t703 = t350 * t777;
t686 = t357 * t736;
t351 = t479 * t686;
t702 = t351 * t775;
t685 = t358 * t730;
t352 = t479 * t685;
t701 = t352 * t773;
t700 = t349 * t416 * t450;
t699 = t349 * t750;
t698 = t349 * t712;
t697 = t350 * t425 * t458;
t696 = t350 * t740;
t695 = t350 * t711;
t694 = t351 * t429 * t460;
t693 = t351 * t734;
t692 = t351 * t710;
t691 = t352 * t433 * t462;
t690 = t352 * t728;
t689 = t352 * t709;
t684 = t404 * t752;
t683 = t405 * t742;
t682 = t406 * t736;
t681 = t407 * t730;
t680 = t408 * t752;
t679 = t408 * t749;
t678 = t409 * t742;
t677 = t409 * t739;
t676 = t410 * t736;
t675 = t410 * t733;
t674 = t411 * t730;
t673 = t411 * t727;
t671 = t422 * t745;
t669 = t416 * t744;
t668 = t417 * t746;
t667 = t420 * t748;
t666 = 0.1e1 / t483 ^ 2 * t748;
t665 = t417 * t712;
t664 = t418 * t712;
t663 = t419 * t712;
t661 = t439 * t723;
t659 = t425 * t722;
t658 = t426 * t724;
t657 = t437 * t738;
t656 = 0.1e1 / t492 ^ 2 * t738;
t655 = t426 * t711;
t653 = t444 * t719;
t651 = t429 * t718;
t650 = t430 * t720;
t649 = t442 * t732;
t648 = 0.1e1 / t495 ^ 2 * t732;
t647 = t430 * t710;
t645 = t449 * t715;
t643 = t433 * t714;
t642 = t434 * t716;
t641 = t447 * t726;
t640 = 0.1e1 / t498 ^ 2 * t726;
t639 = t434 * t709;
t638 = t435 * t711;
t637 = t436 * t711;
t636 = t440 * t710;
t635 = t441 * t710;
t634 = t445 * t709;
t633 = t446 * t709;
t632 = -0.1e1 - t748;
t631 = -0.1e1 - t738;
t630 = -0.1e1 - t732;
t629 = -0.1e1 - t726;
t628 = t418 * t708;
t627 = t435 * t706;
t626 = t440 * t705;
t625 = t445 * t704;
t624 = t349 * t415 * t752;
t623 = t350 * t424 * t742;
t622 = t351 * t428 * t736;
t621 = t352 * t432 * t730;
t620 = t749 * t793;
t619 = t739 * t792;
t618 = t733 * t791;
t617 = t727 * t790;
t616 = t379 * t404 * t747;
t615 = t379 * t699;
t614 = t380 * t408 * t747;
t613 = t380 * t699;
t612 = t381 * t405 * t737;
t611 = t381 * t696;
t610 = t382 * t409 * t737;
t609 = t382 * t696;
t608 = t383 * t406 * t731;
t607 = t383 * t693;
t606 = t384 * t410 * t731;
t605 = t384 * t693;
t604 = t385 * t407 * t725;
t603 = t385 * t690;
t602 = t386 * t411 * t725;
t601 = t386 * t690;
t600 = t408 * t672;
t599 = t408 * t668;
t598 = t409 * t662;
t597 = t409 * t658;
t596 = t410 * t654;
t595 = t410 * t650;
t594 = t411 * t646;
t593 = t411 * t642;
t414 = t450 * t415;
t592 = t414 * t417 * t743;
t590 = t417 * t664;
t589 = t420 * t665;
t423 = t458 * t424;
t588 = t423 * t426 * t721;
t586 = t426 * t638;
t585 = t437 * t655;
t427 = t460 * t428;
t584 = t427 * t430 * t717;
t582 = t430 * t636;
t581 = t442 * t647;
t431 = t462 * t432;
t580 = t431 * t434 * t713;
t578 = t434 * t634;
t577 = t447 * t639;
t576 = t416 * t698;
t575 = t425 * t695;
t574 = t429 * t692;
t573 = t433 * t689;
t559 = t407 * t578;
t560 = t406 * t582;
t561 = t405 * t586;
t562 = t404 * t590;
t572 = (t559 + t560 + t561 + t562) * t479;
t571 = t794 * t479;
t566 = t450 * t628;
t565 = t458 * t627;
t564 = t460 * t626;
t563 = t462 * t625;
t557 = t408 * t590;
t555 = t409 * t586;
t553 = t410 * t582;
t551 = t411 * t578;
t550 = t379 * t408 - t380 * t404;
t549 = t379 * t422 - t760;
t548 = t381 * t409 - t382 * t405;
t547 = t381 * t439 - t759;
t546 = t383 * t410 - t384 * t406;
t545 = t383 * t444 - t758;
t544 = t385 * t411 - t386 * t407;
t543 = t385 * t449 - t757;
t542 = t416 * t671 - t451;
t541 = t425 * t661 - t459;
t540 = t429 * t653 - t461;
t539 = t433 * t645 - t463;
t538 = (-t451 - t751) * t746;
t537 = (-t459 - t741) * t724;
t536 = (-t461 - t735) * t720;
t535 = (-t463 - t729) * t716;
t534 = t550 * t418;
t533 = (-t380 * t422 - t756) * t417;
t532 = t548 * t435;
t531 = (-t382 * t439 - t755) * t426;
t530 = t546 * t440;
t529 = (-t384 * t444 - t754) * t430;
t528 = t544 * t445;
t527 = (-t386 * t449 - t753) * t434;
t526 = t355 * t379 * t417 + t404 * t783;
t525 = -t355 * t770 + t408 * t783;
t524 = t356 * t381 * t426 + t405 * t778;
t523 = -t356 * t767 + t409 * t778;
t522 = t357 * t383 * t430 + t406 * t776;
t521 = -t357 * t764 + t410 * t776;
t520 = t358 * t385 * t434 + t407 * t774;
t519 = -t358 * t761 + t411 * t774;
t518 = t550 * t747;
t517 = t548 * t737;
t516 = t546 * t731;
t515 = t544 * t725;
t514 = t526 * t453;
t513 = t525 * t453;
t512 = t524 * t465;
t511 = t523 * t465;
t510 = t522 * t467;
t509 = t521 * t467;
t508 = t520 * t469;
t507 = t519 * t469;
t506 = t349 * t746 + t350 * t724 + t351 * t720 + t352 * t716;
t505 = t414 * t667 + t423 * t657 + t427 * t649 + t431 * t641;
t504 = t415 * t416 * t743 + t424 * t425 * t721 + t428 * t429 * t717 + t432 * t433 * t713;
t503 = -t453 * t624 - t465 * t623 - t467 * t622 - t469 * t621;
t502 = -t355 * t591 - t356 * t587 - t357 * t583 - t358 * t579;
t501 = -t404 * t591 - t405 * t587 - t406 * t583 - t407 * t579;
t403 = t411 ^ 2;
t402 = t410 ^ 2;
t401 = t409 ^ 2;
t400 = t408 ^ 2;
t399 = t407 ^ 2;
t398 = t406 ^ 2;
t397 = t405 ^ 2;
t396 = t404 ^ 2;
t395 = t412 ^ 2 + t413 ^ 2;
t370 = t539 * t479;
t369 = t540 * t479;
t368 = t541 * t479;
t367 = t542 * t479;
t366 = t479 * t535;
t365 = t479 * t536;
t364 = t479 * t537;
t363 = t479 * t538;
t354 = (-t674 - t676 - t678 - t680) * t480;
t353 = (t681 + t682 + t683 + t684) * t480;
t348 = (-t408 * t669 - t409 * t659 - t410 * t651 - t411 * t643) * t480;
t347 = (t404 * t669 + t405 * t659 + t406 * t651 + t407 * t643) * t480;
t346 = (-t404 * t679 - t405 * t677 - t406 * t675 - t407 * t673) * t480;
t339 = (t408 * t589 + t409 * t585 + t410 * t581 + t411 * t577) * t480;
t338 = (t408 * t592 + t409 * t588 + t410 * t584 + t411 * t580) * t480;
t337 = (-t404 * t589 - t405 * t585 - t406 * t581 - t407 * t577) * t480;
t336 = (-t404 * t592 - t405 * t588 - t406 * t584 - t407 * t580) * t480;
t333 = t794 * t785;
t332 = t501 * t785;
t331 = (-t404 * t600 - t405 * t598 - t406 * t596 - t407 * t594) * t480;
t330 = (-t404 * t599 - t405 * t597 - t406 * t595 - t407 * t593) * t785;
t329 = t380 * t668 + t382 * t658 + t384 * t650 + t386 * t642;
t328 = t379 * t668 + t381 * t658 + t383 * t650 + t385 * t642;
t327 = t379 * t770 + t381 * t767 + t383 * t764 + t385 * t761;
t326 = ((t386 * t469 + t411) * t643 + (t384 * t467 + t410) * t651 + (t382 * t465 + t409) * t659 + (t380 * t453 + t408) * t669) * t479;
t325 = ((t385 * t469 - t407) * t643 + (t383 * t467 - t406) * t651 + (t381 * t465 - t405) * t659 + (t379 * t453 - t404) * t669) * t479;
t324 = (t527 * t714 + t529 * t718 + t531 * t722 + t533 * t744) * t479;
t323 = (-t549 * t417 * t744 - t547 * t426 * t722 - t545 * t430 * t718 - t543 * t434 * t714) * t479;
t322 = (t416 * t534 + t425 * t532 + t429 * t530 + t433 * t528) * t479;
t321 = (-t515 - t516 - t517 - t518) * t479;
t320 = (-t418 * t518 - t435 * t517 - t440 * t516 - t445 * t515) * t479;
t319 = (t528 * t639 + t530 * t647 + t532 * t655 + t534 * t665) * t479;
t1 = [t380 ^ 2 * t417 + t382 ^ 2 * t426 + t384 ^ 2 * t430 + t386 ^ 2 * t434, (t400 * t749 + t401 * t739 + t402 * t733 + t403 * t727) * t480, (-t418 * t614 - t435 * t610 - t440 * t606 - t445 * t602) * t786, (t380 * t680 + t382 * t678 + t384 * t676 + t386 * t674) * t786, (t400 * t672 + t401 * t662 + t402 * t654 + t403 * t646) * t480, (t400 * t668 + t401 * t658 + t402 * t650 + t403 * t642) * t785, 0, 0, 0, (-t602 - t606 - t610 - t614) * t786, (t380 * t557 + t382 * t555 + t384 * t553 + t386 * t551) * t786, 0, 0, 0, t395; t327, t346, t320, t322, t331, t330, 0, 0, 0, t321, t319, 0, 0, 0, 0; t329, t339, t324, t326, t338, t333, t348, t354, 0, t363 * t771 + t364 * t768 + t365 * t765 + t366 * t762 + (-t551 - t553 - t555 - t557) * t479, t367 * t771 + t368 * t768 + t369 * t765 + t370 * t762 + t571, 0, 0, 0, 0; t788, (-t349 * t680 - t350 * t678 - t351 * t676 - t352 * t674) * t479, t613 + t609 + t605 + t601 + (-t408 * t628 - t409 * t627 - t410 * t626 - t411 * t625) * t479, -t380 * t349 - t382 * t350 - t384 * t351 - t386 * t352 + (t408 * t784 + t409 * t781 + t410 * t780 + t411 * t779) * t479, (-t408 * t624 - t409 * t623 - t410 * t622 - t411 * t621) * t479, (-t408 * t700 - t409 * t697 - t410 * t694 - t411 * t691) * t786, 0, 0, 0, t452 * t613 + t464 * t609 + t466 * t605 + t468 * t601 + (-t408 * t708 - t409 * t706 - t410 * t705 - t411 * t704) * t479, -t380 * t576 - t382 * t575 - t384 * t574 - t386 * t573 + (t408 * t566 + t409 * t565 + t410 * t564 + t411 * t563) * t479, 0, -t412, -t413, 0; t327, t346, t320, t322, t331, t330, 0, 0, 0, t321, t319, 0, 0, 0, 0; t379 ^ 2 * t417 + t381 ^ 2 * t426 + t383 ^ 2 * t430 + t385 ^ 2 * t434, (t396 * t749 + t397 * t739 + t398 * t733 + t399 * t727) * t480, (t418 * t616 + t435 * t612 + t440 * t608 + t445 * t604) * t786, (-t379 * t684 - t381 * t683 - t383 * t682 - t385 * t681) * t786, (t396 * t672 + t397 * t662 + t398 * t654 + t399 * t646) * t480, (t396 * t668 + t397 * t658 + t398 * t650 + t399 * t642) * t785, 0, 0, 0, (t604 + t608 + t612 + t616) * t786, (-t379 * t562 - t381 * t561 - t383 * t560 - t385 * t559) * t786, 0, 0, 0, t395; t328, t337, t323, t325, t336, t332, t347, t353, 0, t363 * t772 + t364 * t769 + t365 * t766 + t366 * t763 + t572, t367 * t772 + t368 * t769 + t369 * t766 + t370 * t763 + t479 * t501, 0, 0, 0, 0; t787, (t349 * t684 + t350 * t683 + t351 * t682 + t352 * t681) * t479, t615 + t611 + t607 + t603 + (t404 * t628 + t405 * t627 + t406 * t626 + t407 * t625) * t479, -t379 * t349 - t381 * t350 - t383 * t351 - t385 * t352 + (-t404 * t784 - t405 * t781 - t406 * t780 - t407 * t779) * t479, (t404 * t624 + t405 * t623 + t406 * t622 + t407 * t621) * t479, (t404 * t700 + t405 * t697 + t406 * t694 + t407 * t691) * t786, 0, 0, 0, t452 * t615 + t464 * t611 + t466 * t607 + t468 * t603 + (t404 * t708 + t405 * t706 + t406 * t705 + t407 * t704) * t479, -t379 * t576 - t381 * t575 - t383 * t574 - t385 * t573 + (-t404 * t566 - t405 * t565 - t406 * t564 - t407 * t563) * t479, 0, t413, -t412, 0; t329, t339, t324, t326, t338, t333, t348, t354, 0, ((-t386 + t527) * t716 + (-t384 + t529) * t720 + (-t382 + t531) * t724 + (-t380 + t533) * t746) * t479, (t645 * t761 + t653 * t764 + t661 * t767 + t671 * t770 - t380 - t382 - t384 - t386) * t479 + t571, 0, 0, 0, 0; t328, t337, t323, t325, t336, t332, t347, t353, 0, (t379 * t632 * t746 + t381 * t631 * t724 + t383 * t630 * t720 + t385 * t629 * t716) * t479 + t572, (t543 * t646 + t545 * t654 + t547 * t662 + t549 * t672 - t379 - t381 - t383 - t385) * t479, 0, 0, 0, 0; t646 + t654 + t662 + t672, (t415 * t666 + t424 * t656 + t428 * t648 + t432 * t640) * t480, (-t415 * t667 - t424 * t657 - t428 * t649 - t432 * t641) * t786, t504 * t786, (t415 ^ 2 * t666 + t424 ^ 2 * t656 + t428 ^ 2 * t648 + t432 ^ 2 * t640) * t480, t505 * t785, -0.2e1 * t504 * t480, (-t416 * t663 - t425 * t637 - t429 * t635 - t433 * t633) * t785, (t419 + t436 + t441 + t446) * t480, t363 * t670 + t364 * t660 + t365 * t652 + t366 * t644 + (t629 * t715 + t630 * t719 + t631 * t723 + t632 * t745) * t479, t367 * t670 + t368 * t660 + t369 * t652 + t370 * t644 + (t505 - t716 - t720 - t724 - t746) * t479, 0, 0, 0, 1; t789, (-t419 * t576 - t436 * t575 - t441 * t574 - t446 * t573) * t479, t418 * t576 + t435 * t575 + t440 * t574 + t445 * t573 + (-t422 * t570 * t669 - t439 * t569 * t659 - t444 * t568 * t651 - t449 * t567 * t643) * t479, (t567 * t633 + t568 * t635 + t569 * t637 + t570 * t663) * t479 - t506, (-t414 * t419 * t699 - t423 * t436 * t696 - t427 * t441 * t693 - t431 * t446 * t690) * t479, t503 * t786, t506 * t479, (t349 + t350 + t351 + t352) * t479, 0, t576 + t575 + t574 + t573 + (t535 * t567 + t536 * t568 + t537 * t569 + t538 * t570) * t479, (t539 * t567 + t540 * t568 + t541 * t569 + t542 * t570) * t479 + t503, 0, 0, 0, 0; t788, (-t355 * t679 - t356 * t677 - t357 * t675 - t358 * t673) * t480, (-t418 * t513 - t435 * t511 - t440 * t509 - t445 * t507) * t479, ((-t358 * t762 + t411 * t567) * t445 + (-t357 * t765 + t410 * t568) * t440 + (-t356 * t768 + t409 * t569) * t435 + (-t355 * t771 + t408 * t570) * t418) * t479, (-t355 * t600 - t356 * t598 - t357 * t596 - t358 * t594) * t480, (-t355 * t599 - t356 * t597 - t357 * t595 - t358 * t593) * t785, 0, 0, 0, (-t507 - t509 - t511 - t513) * t479, (t519 * t634 + t521 * t636 + t523 * t638 + t525 * t664) * t479, 0, -t412, -t413, 0; t787, (t617 + t618 + t619 + t620) * t480, (t418 * t514 + t435 * t512 + t440 * t510 + t445 * t508) * t479, ((-t358 * t763 - t407 * t567) * t445 + (-t357 * t766 - t406 * t568) * t440 + (-t356 * t769 - t405 * t569) * t435 + (-t355 * t772 - t404 * t570) * t418) * t479, (t415 * t620 + t424 * t619 + t428 * t618 + t432 * t617) * t480, (t642 * t790 + t650 * t791 + t658 * t792 + t668 * t793) * t785, 0, 0, 0, (t508 + t510 + t512 + t514) * t479, (-t520 * t634 - t522 * t636 - t524 * t638 - t526 * t664) * t479, 0, t413, -t412, 0; t789, (-t355 * t589 - t356 * t585 - t357 * t581 - t358 * t577) * t480, ((t358 * t725 - t567 * t729) * t714 + (t357 * t731 - t568 * t735) * t718 + (t356 * t737 - t569 * t741) * t722 + (t355 * t747 - t570 * t751) * t744) * t479, ((-t358 * t433 + t773) * t714 + (-t357 * t429 + t775) * t718 + (-t356 * t425 + t777) * t722 + (-t355 * t416 + t782) * t744) * t479, (-t355 * t592 - t356 * t588 - t357 * t584 - t358 * t580) * t480, t502 * t785, (t355 * t669 + t356 * t659 + t357 * t651 + t358 * t643) * t480, (t685 + t686 + t687 + t688) * t480, 0, t570 * t363 + t569 * t364 + t568 * t365 + t567 * t366 + (t355 * t590 + t356 * t586 + t357 * t582 + t358 * t578) * t479, t367 * t570 + t368 * t569 + t369 * t568 + t370 * t567 + t479 * t502, 0, 0, 0, 0; t567 ^ 2 + t568 ^ 2 + t569 ^ 2 + t570 ^ 2, (t349 * t688 + t350 * t687 + t351 * t686 + t352 * t685) * t479, t707 + t703 + t702 + t701 + (t355 * t628 + t356 * t627 + t357 * t626 + t358 * t625) * t479, -t570 * t349 * t451 - t569 * t350 * t459 - t568 * t351 * t461 - t567 * t352 * t463 + (-t355 * t784 - t356 * t781 - t357 * t780 - t358 * t779) * t479, (t355 * t624 + t356 * t623 + t357 * t622 + t358 * t621) * t479, (t355 * t700 + t356 * t697 + t357 * t694 + t358 * t691) * t786, 0, 0, 0, t452 * t707 + t464 * t703 + t466 * t702 + t468 * t701 + (t355 * t708 + t356 * t706 + t357 * t705 + t358 * t704) * t479, -t570 * t698 - t569 * t695 - t568 * t692 - t567 * t689 + (-t355 * t566 - t356 * t565 - t357 * t564 - t358 * t563) * t479, 1, 0, 0, 0;];
tau_reg  = t1;
