% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRRR1G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [3x3]
%   Generalized joint coordinates (joint angles)
%   rows: links of the robot
%   columns: number of leg
% legFrame [3x3]
%   base frame orientation for each leg
%   row: number of leg
%   column: Euler angles for the orientation.
%   Euler angle convention from robot definition ("leg_frame")
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3x12]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRRR1G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:27
% EndTime: 2020-03-09 21:02:28
% DurationCPUTime: 1.64s
% Computational Cost: add. (1228->202), mult. (4788->496), div. (1707->17), fcn. (4095->18), ass. (0->219)
t532 = sin(qJ(2,3));
t496 = 0.1e1 / t532;
t537 = cos(qJ(3,3));
t560 = t537 ^ 2;
t510 = 0.1e1 / t560;
t531 = sin(qJ(3,3));
t495 = t531 ^ 2;
t512 = 0.1e1 / t560 ^ 2;
t680 = t495 * t512;
t701 = t496 * (t510 + t680);
t534 = sin(qJ(2,2));
t501 = 0.1e1 / t534;
t539 = cos(qJ(3,2));
t564 = t539 ^ 2;
t516 = 0.1e1 / t564;
t533 = sin(qJ(3,2));
t500 = t533 ^ 2;
t518 = 0.1e1 / t564 ^ 2;
t674 = t500 * t518;
t700 = t501 * (t516 + t674);
t536 = sin(qJ(2,1));
t506 = 0.1e1 / t536;
t541 = cos(qJ(3,1));
t568 = t541 ^ 2;
t522 = 0.1e1 / t568;
t535 = sin(qJ(3,1));
t505 = t535 ^ 2;
t524 = 0.1e1 / t568 ^ 2;
t668 = t505 * t524;
t699 = t506 * (t522 + t668);
t546 = 0.1e1 / pkin(2);
t509 = 0.1e1 / t537;
t515 = 0.1e1 / t539;
t521 = 0.1e1 / t541;
t547 = 0.1e1 / pkin(2) ^ 2;
t497 = 0.1e1 / t532 ^ 2;
t502 = 0.1e1 / t534 ^ 2;
t507 = 0.1e1 / t536 ^ 2;
t528 = legFrame(3,2);
t488 = sin(t528);
t491 = cos(t528);
t544 = xDP(2);
t545 = xDP(1);
t538 = cos(qJ(2,3));
t543 = xDP(3);
t648 = t538 * t543;
t476 = -t531 * t648 + (t488 * t544 - t491 * t545) * t537;
t511 = t509 * t510;
t675 = t497 * t547;
t607 = t512 * t675;
t677 = t496 * t538;
t614 = t476 * t677;
t645 = t543 * t546 ^ 2;
t651 = t531 * t532;
t679 = t496 * t509;
t455 = -(-t543 * t651 + t614) * t476 * t607 - (-t476 * t531 + t648) * t511 * t645 * t679;
t698 = t455 * t538;
t529 = legFrame(2,2);
t489 = sin(t529);
t492 = cos(t529);
t540 = cos(qJ(2,2));
t647 = t540 * t543;
t477 = -t533 * t647 + (t489 * t544 - t492 * t545) * t539;
t517 = t515 * t516;
t669 = t502 * t547;
t605 = t518 * t669;
t671 = t501 * t540;
t612 = t477 * t671;
t650 = t533 * t534;
t673 = t501 * t515;
t456 = -(-t543 * t650 + t612) * t477 * t605 - (-t477 * t533 + t647) * t517 * t645 * t673;
t697 = t456 * t540;
t530 = legFrame(1,2);
t490 = sin(t530);
t493 = cos(t530);
t542 = cos(qJ(2,1));
t646 = t542 * t543;
t478 = -t535 * t646 + (t490 * t544 - t493 * t545) * t541;
t523 = t521 * t522;
t663 = t507 * t547;
t603 = t524 * t663;
t665 = t506 * t542;
t610 = t478 * t665;
t649 = t535 * t536;
t667 = t506 * t521;
t457 = -(-t543 * t649 + t610) * t478 * t603 - (-t478 * t535 + t646) * t523 * t645 * t667;
t696 = t457 * t542;
t527 = t543 ^ 2;
t473 = t476 ^ 2;
t498 = t496 * t497;
t686 = t473 * t498;
t464 = (t496 * t527 + t686) * t546 * t511;
t695 = t464 * t496;
t694 = t464 * t509;
t693 = t464 * t546;
t474 = t477 ^ 2;
t503 = t501 * t502;
t685 = t474 * t503;
t465 = (t501 * t527 + t685) * t546 * t517;
t692 = t465 * t501;
t691 = t465 * t515;
t690 = t465 * t546;
t475 = t478 ^ 2;
t508 = t506 * t507;
t684 = t475 * t508;
t466 = (t506 * t527 + t684) * t546 * t523;
t689 = t466 * t506;
t688 = t466 * t521;
t687 = t466 * t546;
t683 = t476 * t497;
t682 = t477 * t502;
t681 = t478 * t507;
t514 = t538 ^ 2;
t678 = t496 * t514;
t676 = t497 * t512;
t520 = t540 ^ 2;
t672 = t501 * t520;
t670 = t502 * t518;
t526 = t542 ^ 2;
t666 = t506 * t526;
t664 = t507 * t524;
t662 = t509 * t531;
t661 = t509 * t538;
t660 = t511 * t531;
t659 = t515 * t533;
t658 = t515 * t540;
t657 = t517 * t533;
t656 = t521 * t535;
t655 = t521 * t542;
t654 = t523 * t535;
t653 = t527 * t547;
t548 = t546 * t547;
t652 = t527 * t548;
t644 = t543 * t548;
t470 = t473 * t607;
t467 = t510 * t653 + t470;
t636 = -0.2e1 * t543 * t547;
t574 = t614 * t636;
t602 = t532 * t653;
t643 = -t495 * t511 * t602 + t574 * t660 + (-t467 * t532 + t698) * t537;
t471 = t474 * t605;
t468 = t516 * t653 + t471;
t573 = t612 * t636;
t601 = t534 * t653;
t642 = -t500 * t517 * t601 + t573 * t657 + (-t468 * t534 + t697) * t539;
t472 = t475 * t603;
t469 = t522 * t653 + t472;
t572 = t610 * t636;
t600 = t536 * t653;
t641 = -t505 * t523 * t600 + t572 * t654 + (-t469 * t536 + t696) * t541;
t640 = (-t510 * t602 - t698) * t531 + t467 * t651 + t510 * t574;
t639 = (-t516 * t601 - t697) * t533 + t468 * t650 + t516 * t573;
t638 = (-t522 * t600 - t696) * t535 + t469 * t649 + t522 * t572;
t637 = 0.2e1 * t543;
t635 = t455 * t679;
t634 = t455 * t496 * t531;
t633 = t455 * t677;
t632 = t456 * t673;
t631 = t456 * t501 * t533;
t630 = t456 * t671;
t629 = t457 * t667;
t628 = t457 * t506 * t535;
t627 = t457 * t665;
t626 = t510 * t693;
t625 = t538 * t693;
t624 = t516 * t690;
t623 = t540 * t690;
t622 = t522 * t687;
t621 = t542 * t687;
t620 = t473 * t676;
t513 = t509 * t512;
t619 = t473 * t513 * t547;
t618 = t474 * t670;
t519 = t515 * t518;
t617 = t474 * t519 * t547;
t616 = t475 * t664;
t525 = t521 * t524;
t615 = t475 * t525 * t547;
t613 = t538 * t683;
t611 = t540 * t682;
t609 = t542 * t681;
t608 = t496 * t662;
t606 = t501 * t659;
t604 = t506 * t656;
t599 = t548 * t637;
t595 = t495 * t635;
t594 = t510 * t633;
t593 = t500 * t632;
t592 = t516 * t630;
t591 = t505 * t629;
t590 = t522 * t627;
t589 = t661 * t695;
t588 = t658 * t692;
t587 = t655 * t689;
t586 = t512 * t538 * t686;
t585 = t518 * t540 * t685;
t584 = t524 * t542 * t684;
t485 = -0.1e1 + 0.2e1 * t560;
t583 = t476 * t485 * t676;
t582 = t660 * t683;
t486 = -0.1e1 + 0.2e1 * t564;
t581 = t477 * t486 * t670;
t580 = t657 * t682;
t487 = -0.1e1 + 0.2e1 * t568;
t579 = t478 * t487 * t664;
t578 = t654 * t681;
t577 = t625 * t662;
t576 = t623 * t659;
t575 = t621 * t656;
t504 = t535 * t505;
t499 = t533 * t500;
t494 = t531 * t495;
t484 = t490 * t536 + t493 * t542;
t483 = -t490 * t542 + t493 * t536;
t482 = t489 * t534 + t492 * t540;
t481 = -t489 * t540 + t492 * t534;
t480 = t488 * t532 + t491 * t538;
t479 = -t488 * t538 + t491 * t532;
t1 = [t480 * t695 + t482 * t692 + t484 * t689, (-t491 * t635 - t492 * t632 - t493 * t629) * t546, t480 * t633 + t482 * t630 + t484 * t627 + (-t480 * t620 - t482 * t618 - t484 * t616) * t547 + (-t491 * t589 - t492 * t588 - t493 * t587) * t546, -t480 * t455 - t482 * t456 - t484 * t457 + (-t480 * t586 - t482 * t585 - t484 * t584) * t547 + (t491 * t694 + t492 * t691 + t493 * t688) * t546, (-t491 * t595 - t492 * t593 - t493 * t591) * t546 + (-t491 * t582 - t492 * t580 - t493 * t578) * t599, 0.2e1 * (-t491 * t634 - t492 * t631 - t493 * t628) * t546 + 0.2e1 * (-t491 * t583 - t492 * t581 - t493 * t579) * t644, (-t491 * t701 - t492 * t700 - t493 * t699) * t652, 0, 0, (t641 * t484 - t493 * t621) * t506 + (t642 * t482 - t492 * t623) * t501 + (t643 * t480 - t491 * t625) * t496, (t638 * t484 + t493 * t575) * t506 + (t639 * t482 + t492 * t576) * t501 + (t640 * t480 + t491 * t577) * t496, 0; t479 * t695 + t481 * t692 + t483 * t689, (t488 * t635 + t489 * t632 + t490 * t629) * t546, t479 * t633 + t481 * t630 + t483 * t627 + (-t479 * t620 - t481 * t618 - t483 * t616) * t547 + (t488 * t589 + t489 * t588 + t490 * t587) * t546, -t479 * t455 - t481 * t456 - t483 * t457 + (-t479 * t586 - t481 * t585 - t483 * t584) * t547 + (-t488 * t694 - t489 * t691 - t490 * t688) * t546, (t488 * t595 + t489 * t593 + t490 * t591) * t546 + (t488 * t582 + t489 * t580 + t490 * t578) * t599, 0.2e1 * (t488 * t634 + t489 * t631 + t490 * t628) * t546 + 0.2e1 * (t488 * t583 + t489 * t581 + t490 * t579) * t644, (t488 * t701 + t489 * t700 + t490 * t699) * t652, 0, 0, (t641 * t483 + t490 * t621) * t506 + (t642 * t481 + t489 * t623) * t501 + (t643 * t479 + t488 * t625) * t496, (t638 * t483 - t490 * t575) * t506 + (t639 * t481 - t489 * t576) * t501 + (t640 * t479 - t488 * t577) * t496, 0; t464 * t608 + t465 * t606 + t466 * t604, (-t531 * t594 - t533 * t592 - t535 * t590) * t546, (-t507 * t615 + (t457 * t655 - t526 * t622) * t506) * t535 + (-t502 * t617 + (t456 * t658 - t520 * t624) * t501) * t533 + (-t497 * t619 + (t455 * t661 - t514 * t626) * t496) * t531, (-t457 * t521 + (-t508 * t615 + t622) * t542) * t535 + (-t456 * t515 + (-t503 * t617 + t624) * t540) * t533 + (-t455 * t509 + (-t498 * t619 + t626) * t538) * t531, (-t494 * t594 - t499 * t592 - t504 * t590) * t546 + (-t531 * t620 - t533 * t618 - t535 * t616 + (-t609 * t668 - t611 * t674 - t613 * t680) * t637) * t548, (-t485 * t513 * t531 * t613 - t486 * t519 * t533 * t611 - t487 * t525 * t535 * t609) * t599 + ((-0.2e1 * t475 * t522 * t663 - 0.2e1 * t505 * t627 + t472) * t521 + (-0.2e1 * t474 * t516 * t669 - 0.2e1 * t500 * t630 + t471) * t515 + (-0.2e1 * t473 * t510 * t675 - 0.2e1 * t495 * t633 + t470) * t509) * t546, (t455 * t662 + t456 * t659 + t457 * t656) * t546 + ((-t504 * t525 - t654) * t665 + (-t499 * t519 - t657) * t671 + (-t494 * t513 - t660) * t677) * t652, (t455 + t456 + t457) * t546, (t512 * t531 + t518 * t533 + t524 * t535) * t652, t641 * t604 + t642 * t606 + t643 * t608 + ((-t536 - t666) * t466 * t656 + (-t534 - t672) * t465 * t659 + (-t532 - t678) * t464 * t662) * t546, t638 * t604 + t639 * t606 + t640 * t608 + ((t505 * t522 * t666 - t536) * t466 + (t500 * t516 * t672 - t534) * t465 + (t495 * t510 * t678 - t532) * t464) * t546, 0;];
tau_reg  = t1;
