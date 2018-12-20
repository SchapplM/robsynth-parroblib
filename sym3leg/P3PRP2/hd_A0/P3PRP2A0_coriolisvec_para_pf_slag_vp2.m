% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRP2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (including platform)
% mrSges [4x3]
%   first moment of all robot links (mass times center of mass in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taucX = P3PRP2A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:39:01
% EndTime: 2018-12-20 17:39:05
% DurationCPUTime: 3.80s
% Computational Cost: add. (26233->373), mult. (49823->591), div. (2250->3), fcn. (23653->14), ass. (0->258)
t612 = (qJ(3,1) ^ 2);
t620 = (pkin(2) ^ 2);
t708 = -t620 - 1;
t572 = -t612 - t708;
t603 = cos(qJ(2,1));
t589 = t603 ^ 2;
t600 = sin(qJ(2,1));
t663 = t600 * t603;
t647 = pkin(2) * t663;
t657 = t612 + t620;
t715 = 2 * qJ(3,1);
t534 = 0.1e1 / (t572 * t589 + t647 * t715 - t657 - 0.1e1);
t594 = legFrame(1,3);
t581 = cos(t594);
t578 = sin(t594);
t703 = qJ(3,1) * t581;
t652 = pkin(2) * t703;
t627 = t578 * t612 - t652;
t704 = qJ(3,1) * t578;
t562 = pkin(2) * t704;
t660 = t612 * t581 + t562;
t525 = t627 * t603 - t600 * (t581 + t660);
t607 = xP(3);
t583 = sin(t607);
t584 = cos(t607);
t615 = koppelP(1,2);
t618 = koppelP(1,1);
t558 = -t583 * t615 + t584 * t618;
t604 = xDP(3);
t605 = xDP(2);
t537 = t558 * t604 + t605;
t681 = t525 * t537;
t522 = t660 * t603 - t600 * (-t578 - t627);
t555 = t583 * t618 + t584 * t615;
t606 = xDP(1);
t540 = -t555 * t604 + t606;
t684 = t522 * t540;
t719 = (t681 + t684) * t534;
t611 = (qJ(3,2) ^ 2);
t571 = -t611 - t708;
t602 = cos(qJ(2,2));
t588 = t602 ^ 2;
t599 = sin(qJ(2,2));
t664 = t599 * t602;
t648 = pkin(2) * t664;
t658 = t611 + t620;
t714 = 2 * qJ(3,2);
t533 = 0.1e1 / (t571 * t588 + t648 * t714 - t658 - 0.1e1);
t593 = legFrame(2,3);
t580 = cos(t593);
t577 = sin(t593);
t700 = qJ(3,2) * t580;
t651 = pkin(2) * t700;
t626 = t577 * t611 - t651;
t701 = qJ(3,2) * t577;
t561 = pkin(2) * t701;
t661 = t611 * t580 + t561;
t524 = t626 * t602 - t599 * (t580 + t661);
t614 = koppelP(2,2);
t617 = koppelP(2,1);
t557 = -t583 * t614 + t584 * t617;
t536 = t557 * t604 + t605;
t682 = t524 * t536;
t521 = t661 * t602 - t599 * (-t577 - t626);
t554 = t583 * t617 + t584 * t614;
t539 = -t554 * t604 + t606;
t685 = t521 * t539;
t718 = (t682 + t685) * t533;
t610 = (qJ(3,3) ^ 2);
t570 = -t610 - t708;
t601 = cos(qJ(2,3));
t587 = t601 ^ 2;
t598 = sin(qJ(2,3));
t665 = t598 * t601;
t649 = pkin(2) * t665;
t659 = t610 + t620;
t713 = 2 * qJ(3,3);
t532 = 0.1e1 / (t570 * t587 + t649 * t713 - t659 - 0.1e1);
t592 = legFrame(3,3);
t579 = cos(t592);
t576 = sin(t592);
t697 = qJ(3,3) * t579;
t650 = pkin(2) * t697;
t625 = t576 * t610 - t650;
t698 = qJ(3,3) * t576;
t560 = pkin(2) * t698;
t662 = t610 * t579 + t560;
t523 = t625 * t601 - t598 * (t579 + t662);
t613 = koppelP(3,2);
t616 = koppelP(3,1);
t556 = -t583 * t613 + t584 * t616;
t535 = t556 * t604 + t605;
t683 = t523 * t535;
t520 = t662 * t601 - t598 * (-t576 - t625);
t553 = t583 * t616 + t584 * t613;
t538 = -t553 * t604 + t606;
t686 = t520 * t538;
t717 = (t683 + t686) * t532;
t716 = 0.2e1 * pkin(2);
t590 = t604 ^ 2;
t696 = qJ(3,3) * t601;
t527 = -0.2e1 * t579 * t696 + t598 * (pkin(2) * t579 - t698);
t679 = t527 * t532;
t526 = 0.2e1 * t576 * t696 - t598 * (pkin(2) * t576 + t697);
t680 = t526 * t532;
t502 = t535 * t679 + t538 * t680;
t712 = 0.2e1 * t502;
t699 = qJ(3,2) * t602;
t529 = -0.2e1 * t580 * t699 + t599 * (pkin(2) * t580 - t701);
t677 = t529 * t533;
t528 = 0.2e1 * t577 * t699 - t599 * (pkin(2) * t577 + t700);
t678 = t528 * t533;
t503 = t536 * t677 + t539 * t678;
t711 = 0.2e1 * t503;
t702 = qJ(3,1) * t603;
t531 = -0.2e1 * t581 * t702 + t600 * (pkin(2) * t581 - t704);
t675 = t531 * t534;
t530 = 0.2e1 * t578 * t702 - t600 * (pkin(2) * t578 + t703);
t676 = t530 * t534;
t504 = t537 * t675 + t540 * t676;
t710 = 0.2e1 * t504;
t709 = m(3) * pkin(2);
t707 = m(3) * t601;
t706 = m(3) * t602;
t705 = m(3) * t603;
t573 = m(3) * qJ(3,3) + mrSges(3,3);
t574 = m(3) * qJ(3,2) + mrSges(3,3);
t575 = m(3) * qJ(3,1) + mrSges(3,3);
t619 = pkin(2) * t620;
t490 = t610 * t502;
t636 = -t620 * t502 + t716 * t717 - t490;
t656 = -t620 + t708;
t689 = t502 * t532;
t439 = (t636 * t696 + ((t619 + (1 + t610) * pkin(2)) * t502 - t717 + t656 * t717) * t598) * t689;
t497 = pkin(2) * t502;
t461 = t497 - t717;
t442 = ((t461 - t717) * t665 + (-t587 * t502 + t502 - t636) * qJ(3,3)) * t689;
t668 = t587 * qJ(3,3);
t446 = ((0.2e1 * (t497 + (-t686 / 0.2e1 - t683 / 0.2e1) * t532) * t668 - (pkin(2) * t461 - t490) * t665 - qJ(3,3) * t717) * t532 + (-qJ(3,3) + t649 - t668) * t532 * t717) * t502;
t582 = mrSges(3,1) + t709;
t436 = t582 * t446 + (t439 * t601 - t442) * m(3);
t695 = t436 * t532;
t491 = t611 * t503;
t635 = -t620 * t503 + t716 * t718 - t491;
t688 = t503 * t533;
t440 = (t635 * t699 + ((t619 + (1 + t611) * pkin(2)) * t503 - t718 + t656 * t718) * t599) * t688;
t498 = pkin(2) * t503;
t462 = t498 - t718;
t443 = ((t462 - t718) * t664 + (-t588 * t503 + t503 - t635) * qJ(3,2)) * t688;
t667 = t588 * qJ(3,2);
t447 = ((0.2e1 * (t498 + (-t685 / 0.2e1 - t682 / 0.2e1) * t533) * t667 - (pkin(2) * t462 - t491) * t664 - qJ(3,2) * t718) * t533 + (-qJ(3,2) + t648 - t667) * t533 * t718) * t503;
t437 = t582 * t447 + (t440 * t602 - t443) * m(3);
t694 = t437 * t533;
t492 = t612 * t504;
t634 = -t620 * t504 + t716 * t719 - t492;
t687 = t504 * t534;
t441 = (t634 * t702 + ((t619 + (1 + t612) * pkin(2)) * t504 - t719 + t656 * t719) * t600) * t687;
t493 = t504 * pkin(2);
t460 = t493 - t719;
t444 = ((t460 - t719) * t663 + (-t589 * t504 + t504 - t634) * qJ(3,1)) * t687;
t666 = t589 * qJ(3,1);
t445 = ((0.2e1 * (t493 + (-t684 / 0.2e1 - t681 / 0.2e1) * t534) * t666 - (pkin(2) * t460 - t492) * t663 - qJ(3,1) * t719) * t534 + (-qJ(3,1) + t647 - t666) * t534 * t719) * t504;
t438 = t582 * t445 + (t441 * t603 - t444) * m(3);
t693 = t438 * t534;
t692 = t502 ^ 2 * t573;
t691 = t503 ^ 2 * t574;
t690 = t504 ^ 2 * t575;
t566 = -mrSges(2,2) + t573;
t569 = mrSges(2,1) + t582;
t541 = t566 * t598 + t569 * t601;
t585 = m(1) + m(2) + m(3);
t433 = -t585 * t439 + t442 * t707 - t541 * t446;
t674 = t532 * t433;
t673 = t532 * t590;
t567 = -mrSges(2,2) + t574;
t542 = t567 * t599 + t569 * t602;
t434 = -t585 * t440 + t443 * t706 - t542 * t447;
t672 = t533 * t434;
t671 = t533 * t590;
t568 = -mrSges(2,2) + t575;
t543 = t568 * t600 + t569 * t603;
t435 = -t585 * t441 + t444 * t705 - t543 * t445;
t670 = t534 * t435;
t669 = t534 * t590;
t559 = -t709 / 0.2e1 - mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1;
t655 = ((m(3) * t717 + t559 * t502) * t598 + t601 * t502 * t566 / 0.2e1) * t712;
t654 = ((m(3) * t718 + t559 * t503) * t599 + t602 * t503 * t567 / 0.2e1) * t711;
t653 = ((m(3) * t719 + t559 * t504) * t600 + t603 * t504 * t568 / 0.2e1) * t710;
t646 = t532 * t692;
t645 = t533 * t691;
t644 = t534 * t690;
t643 = t532 * t655;
t642 = t533 * t654;
t641 = t534 * t653;
t640 = t717 * t573 * t712;
t639 = t718 * t574 * t711;
t638 = t719 * t575 * t710;
t637 = mrSges(3,1) * t716 + Ifges(3,2) + Ifges(2,3);
t633 = t532 * t640;
t632 = t533 * t639;
t631 = t534 * t638;
t609 = mrSges(4,1);
t608 = mrSges(4,2);
t552 = m(3) * t657 + mrSges(3,3) * t715 + t637;
t551 = m(3) * t658 + mrSges(3,3) * t714 + t637;
t550 = m(3) * t659 + mrSges(3,3) * t713 + t637;
t549 = t572 * t578 + 0.2e1 * t652;
t548 = t571 * t577 + 0.2e1 * t651;
t547 = t570 * t576 + 0.2e1 * t650;
t546 = t581 * t572 - 0.2e1 * t562;
t545 = t580 * t571 - 0.2e1 * t561;
t544 = t579 * t570 - 0.2e1 * t560;
t519 = -t546 * t663 + t549 * t589 - t578 * t620 - t578 - t652;
t518 = -t545 * t664 + t548 * t588 - t577 * t620 - t577 - t651;
t517 = -t544 * t665 + t547 * t587 - t576 * t620 - t576 - t650;
t516 = t546 * t589 + t549 * t663 - t620 * t581 + t562 - t581;
t515 = t545 * t588 + t548 * t664 - t620 * t580 + t561 - t580;
t514 = t544 * t587 + t547 * t665 - t620 * t579 + t560 - t579;
t507 = (-t530 * t555 + t531 * t558) * t534;
t506 = (-t528 * t554 + t529 * t557) * t533;
t505 = (-t526 * t553 + t527 * t556) * t532;
t489 = (-t522 * t555 + t525 * t558) * t534;
t488 = (-t521 * t554 + t524 * t557) * t533;
t487 = (-t520 * t553 + t523 * t556) * t532;
t483 = (-t516 * t555 + t519 * t558) * t534;
t482 = (-t515 * t554 + t518 * t557) * t533;
t481 = (-t514 * t553 + t517 * t556) * t532;
t480 = (-t531 * t582 + (-t519 * t603 + t525) * m(3)) * t534;
t479 = (-t529 * t582 + (-t518 * t602 + t524) * m(3)) * t533;
t478 = (-t527 * t582 + (-t517 * t601 + t523) * m(3)) * t532;
t477 = (-t530 * t582 + (-t516 * t603 + t522) * m(3)) * t534;
t476 = (-t528 * t582 + (-t515 * t602 + t521) * m(3)) * t533;
t475 = (-t526 * t582 + (-t514 * t601 + t520) * m(3)) * t532;
t474 = (t519 * t585 - t525 * t705 + t531 * t543) * t534;
t473 = (t518 * t585 - t524 * t706 + t529 * t542) * t533;
t472 = (t517 * t585 - t523 * t707 + t527 * t541) * t532;
t471 = (t516 * t585 - t522 * t705 + t530 * t543) * t534;
t470 = (t515 * t585 - t521 * t706 + t528 * t542) * t533;
t469 = (t514 * t585 - t520 * t707 + t526 * t541) * t532;
t468 = (t519 * t543 - t525 * t582 + t531 * t552) * t534;
t467 = (t518 * t542 - t524 * t582 + t529 * t551) * t533;
t466 = (t517 * t541 - t523 * t582 + t527 * t550) * t532;
t465 = (t516 * t543 - t522 * t582 + t530 * t552) * t534;
t464 = (t515 * t542 - t521 * t582 + t528 * t551) * t533;
t463 = (t514 * t541 - t520 * t582 + t526 * t550) * t532;
t456 = -t507 * t582 + (-t483 * t603 + t489) * m(3);
t455 = -t506 * t582 + (-t482 * t602 + t488) * m(3);
t454 = -t505 * t582 + (-t481 * t601 + t487) * m(3);
t453 = t483 * t585 - t489 * t705 + t507 * t543;
t452 = t482 * t585 - t488 * t706 + t506 * t542;
t451 = t481 * t585 - t487 * t707 + t505 * t541;
t450 = t483 * t543 - t489 * t582 + t507 * t552;
t449 = t482 * t542 - t488 * t582 + t506 * t551;
t448 = t481 * t541 - t487 * t582 + t505 * t550;
t432 = -t543 * t441 + t582 * t444 - t552 * t445;
t431 = -t542 * t440 + t582 * t443 - t551 * t447;
t430 = -t541 * t439 + t582 * t442 - t550 * t446;
t1 = [(-(t465 * t530 + t471 * t516 + t477 * t522) * t558 - (t465 * t531 + t471 * t519 + t477 * t525) * t555) * t669 + t516 * t670 + t432 * t676 + t522 * t693 + t516 * t641 + t530 * t631 - t522 * t644 + (-(t464 * t528 + t470 * t515 + t476 * t521) * t557 - (t464 * t529 + t470 * t518 + t476 * t524) * t554) * t671 + t515 * t672 + t431 * t678 + t521 * t694 + t515 * t642 + t528 * t632 - t521 * t645 + (-(t463 * t526 + t469 * t514 + t475 * t520) * t556 - (t463 * t527 + t469 * t517 + t475 * t523) * t553) * t673 + t514 * t674 + t430 * t680 + t520 * t695 + t514 * t643 + t526 * t633 - t520 * t646 - t590 * (-t583 * t608 + t584 * t609); (-(t468 * t530 + t474 * t516 + t480 * t522) * t558 - (t468 * t531 + t474 * t519 + t480 * t525) * t555) * t669 + t519 * t670 + t432 * t675 + t525 * t693 + t519 * t641 + t531 * t631 - t525 * t644 + (-(t467 * t528 + t473 * t515 + t479 * t521) * t557 - (t467 * t529 + t473 * t518 + t479 * t524) * t554) * t671 + t518 * t672 + t431 * t677 + t524 * t694 + t518 * t642 + t529 * t632 - t524 * t645 + (-(t466 * t526 + t472 * t514 + t478 * t520) * t556 - (t466 * t527 + t472 * t517 + t478 * t523) * t553) * t673 + t517 * t674 + t430 * t679 + t523 * t695 + t517 * t643 + t527 * t633 - t523 * t646 - t590 * (t583 * t609 + t584 * t608); (-(t450 * t530 + t453 * t516 + t456 * t522) * t558 - (t450 * t531 + t453 * t519 + t456 * t525) * t555) * t669 + (-(t449 * t528 + t452 * t515 + t455 * t521) * t557 - (t449 * t529 + t452 * t518 + t455 * t524) * t554) * t671 + (-(t448 * t526 + t451 * t514 + t454 * t520) * t556 - (t448 * t527 + t451 * t517 + t454 * t523) * t553) * t673 + (t432 + t638) * t507 + (t431 + t639) * t506 + (t430 + t640) * t505 + (t438 - t690) * t489 + (t437 - t691) * t488 + (t436 - t692) * t487 + (t435 + t653) * t483 + (t434 + t654) * t482 + (t433 + t655) * t481;];
taucX  = t1;
