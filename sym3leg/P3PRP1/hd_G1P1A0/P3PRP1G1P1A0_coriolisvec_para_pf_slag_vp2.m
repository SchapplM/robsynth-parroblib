% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRP1G1P1A0
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
%   mass of all robot links (leg links until cut joint, platform)
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
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1G1P1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:41:44
% EndTime: 2019-05-03 14:41:48
% DurationCPUTime: 3.95s
% Computational Cost: add. (26233->373), mult. (49823->588), div. (2250->3), fcn. (23653->14), ass. (0->261)
t608 = (pkin(2) ^ 2);
t579 = 1 + t608;
t600 = (qJ(3,1) ^ 2);
t560 = -t600 + t579;
t591 = cos(qJ(2,1));
t577 = t591 ^ 2;
t588 = sin(qJ(2,1));
t660 = t588 * t591;
t643 = pkin(2) * t660;
t657 = t600 + t608;
t711 = 2 * qJ(3,1);
t531 = 0.1e1 / (t560 * t577 + t643 * t711 - t657 - 0.1e1);
t582 = legFrame(1,3);
t566 = sin(t582);
t569 = cos(t582);
t701 = qJ(3,1) * t566;
t649 = pkin(2) * t701;
t619 = t569 * t600 - t649;
t700 = qJ(3,1) * t569;
t644 = pkin(2) * t700;
t620 = -t566 * t600 - t644;
t522 = t619 * t591 - t588 * (t566 - t620);
t595 = xP(3);
t571 = sin(t595);
t572 = cos(t595);
t603 = koppelP(1,2);
t606 = koppelP(1,1);
t552 = -t571 * t603 + t572 * t606;
t592 = xDP(3);
t593 = xDP(2);
t534 = t552 * t592 + t593;
t672 = t522 * t534;
t519 = t620 * t591 + t588 * (-t569 - t619);
t549 = t571 * t606 + t572 * t603;
t594 = xDP(1);
t537 = -t549 * t592 + t594;
t675 = t519 * t537;
t715 = (t672 + t675) * t531;
t599 = (qJ(3,2) ^ 2);
t559 = -t599 + t579;
t590 = cos(qJ(2,2));
t576 = t590 ^ 2;
t587 = sin(qJ(2,2));
t661 = t587 * t590;
t641 = pkin(2) * t661;
t658 = t599 + t608;
t710 = 2 * qJ(3,2);
t530 = 0.1e1 / (t559 * t576 + t641 * t710 - t658 - 0.1e1);
t581 = legFrame(2,3);
t565 = sin(t581);
t568 = cos(t581);
t697 = qJ(3,2) * t565;
t648 = pkin(2) * t697;
t616 = t568 * t599 - t648;
t696 = qJ(3,2) * t568;
t645 = pkin(2) * t696;
t617 = -t565 * t599 - t645;
t521 = t616 * t590 - t587 * (t565 - t617);
t602 = koppelP(2,2);
t605 = koppelP(2,1);
t551 = -t571 * t602 + t572 * t605;
t533 = t551 * t592 + t593;
t673 = t521 * t533;
t518 = t617 * t590 + t587 * (-t568 - t616);
t548 = t571 * t605 + t572 * t602;
t536 = -t548 * t592 + t594;
t676 = t518 * t536;
t714 = (t673 + t676) * t530;
t598 = (qJ(3,3) ^ 2);
t558 = -t598 + t579;
t589 = cos(qJ(2,3));
t575 = t589 ^ 2;
t586 = sin(qJ(2,3));
t662 = t586 * t589;
t642 = pkin(2) * t662;
t659 = t598 + t608;
t709 = 2 * qJ(3,3);
t529 = 0.1e1 / (t558 * t575 + t642 * t709 - t659 - 0.1e1);
t580 = legFrame(3,3);
t564 = sin(t580);
t567 = cos(t580);
t693 = qJ(3,3) * t564;
t647 = pkin(2) * t693;
t613 = t567 * t598 - t647;
t692 = qJ(3,3) * t567;
t646 = pkin(2) * t692;
t614 = -t564 * t598 - t646;
t520 = t613 * t589 - t586 * (t564 - t614);
t601 = koppelP(3,2);
t604 = koppelP(3,1);
t550 = -t571 * t601 + t572 * t604;
t532 = t550 * t592 + t593;
t674 = t520 * t532;
t517 = t614 * t589 + t586 * (-t567 - t613);
t547 = t571 * t604 + t572 * t601;
t535 = -t547 * t592 + t594;
t677 = t517 * t535;
t713 = (t674 + t677) * t529;
t712 = 0.2e1 * pkin(2);
t578 = t592 ^ 2;
t690 = qJ(3,3) * t589;
t653 = -0.2e1 * t690;
t524 = t567 * t653 + t586 * (pkin(2) * t567 + t693);
t670 = t524 * t529;
t523 = t564 * t653 + t586 * (pkin(2) * t564 - t692);
t671 = t523 * t529;
t499 = t532 * t671 + t535 * t670;
t708 = 0.2e1 * t499;
t694 = qJ(3,2) * t590;
t654 = -0.2e1 * t694;
t526 = t568 * t654 + t587 * (pkin(2) * t568 + t697);
t668 = t526 * t530;
t525 = t565 * t654 + t587 * (pkin(2) * t565 - t696);
t669 = t525 * t530;
t500 = t533 * t669 + t536 * t668;
t707 = 0.2e1 * t500;
t698 = qJ(3,1) * t591;
t655 = -0.2e1 * t698;
t528 = t569 * t655 + t588 * (pkin(2) * t569 + t701);
t666 = t528 * t531;
t527 = t566 * t655 + t588 * (pkin(2) * t566 - t700);
t667 = t527 * t531;
t501 = t534 * t667 + t537 * t666;
t706 = 0.2e1 * t501;
t705 = m(3) * pkin(2);
t704 = m(3) * t589;
t703 = m(3) * t590;
t702 = m(3) * t591;
t561 = m(3) * qJ(3,3) + mrSges(3,3);
t562 = m(3) * qJ(3,2) + mrSges(3,3);
t563 = m(3) * qJ(3,1) + mrSges(3,3);
t699 = qJ(3,1) * t577;
t695 = qJ(3,2) * t576;
t691 = qJ(3,3) * t575;
t607 = pkin(2) * t608;
t487 = t598 * t499;
t630 = -t608 * t499 + t712 * t713 - t487;
t656 = -t608 - t579;
t680 = t499 * t529;
t436 = (t630 * t690 + ((t607 + (1 + t598) * pkin(2)) * t499 - t713 + t656 * t713) * t586) * t680;
t494 = pkin(2) * t499;
t458 = t494 - t713;
t439 = ((t458 - t713) * t662 + (-t499 * t575 + t499 - t630) * qJ(3,3)) * t680;
t443 = ((0.2e1 * (t494 + (-t677 / 0.2e1 - t674 / 0.2e1) * t529) * t691 - (pkin(2) * t458 - t487) * t662 - qJ(3,3) * t713) * t529 + (-qJ(3,3) + t642 - t691) * t529 * t713) * t499;
t554 = -mrSges(2,2) + t561;
t570 = mrSges(3,1) + t705;
t557 = mrSges(2,1) + t570;
t538 = t554 * t586 + t557 * t589;
t573 = m(1) + m(2) + m(3);
t430 = -t436 * t573 + t439 * t704 - t443 * t538;
t689 = t430 * t529;
t488 = t599 * t500;
t629 = -t608 * t500 + t712 * t714 - t488;
t679 = t500 * t530;
t437 = (t629 * t694 + ((t607 + (1 + t599) * pkin(2)) * t500 - t714 + t656 * t714) * t587) * t679;
t495 = pkin(2) * t500;
t459 = t495 - t714;
t440 = ((t459 - t714) * t661 + (-t500 * t576 + t500 - t629) * qJ(3,2)) * t679;
t444 = ((0.2e1 * (t495 + (-t676 / 0.2e1 - t673 / 0.2e1) * t530) * t695 - (pkin(2) * t459 - t488) * t661 - qJ(3,2) * t714) * t530 + (-qJ(3,2) + t641 - t695) * t530 * t714) * t500;
t555 = -mrSges(2,2) + t562;
t539 = t555 * t587 + t557 * t590;
t431 = -t437 * t573 + t440 * t703 - t444 * t539;
t688 = t431 * t530;
t489 = t600 * t501;
t628 = -t608 * t501 + t712 * t715 - t489;
t678 = t501 * t531;
t438 = (t628 * t698 + ((t607 + (1 + t600) * pkin(2)) * t501 - t715 + t656 * t715) * t588) * t678;
t490 = t501 * pkin(2);
t457 = t490 - t715;
t441 = ((t457 - t715) * t660 + (-t501 * t577 + t501 - t628) * qJ(3,1)) * t678;
t442 = ((0.2e1 * (t490 + (-t675 / 0.2e1 - t672 / 0.2e1) * t531) * t699 - (pkin(2) * t457 - t489) * t660 - qJ(3,1) * t715) * t531 + (-qJ(3,1) + t643 - t699) * t531 * t715) * t501;
t556 = -mrSges(2,2) + t563;
t540 = t556 * t588 + t557 * t591;
t432 = -t438 * t573 + t441 * t702 - t442 * t540;
t687 = t432 * t531;
t433 = t570 * t443 + (t436 * t589 - t439) * m(3);
t686 = t433 * t529;
t434 = t570 * t444 + (t437 * t590 - t440) * m(3);
t685 = t434 * t530;
t435 = t570 * t442 + (t438 * t591 - t441) * m(3);
t684 = t435 * t531;
t683 = t499 ^ 2 * t561;
t682 = t500 ^ 2 * t562;
t681 = t501 ^ 2 * t563;
t665 = t529 * t578;
t664 = t530 * t578;
t663 = t531 * t578;
t553 = -t705 / 0.2e1 - mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1;
t652 = ((m(3) * t713 + t499 * t553) * t586 + t589 * t499 * t554 / 0.2e1) * t708;
t651 = ((m(3) * t714 + t500 * t553) * t587 + t590 * t500 * t555 / 0.2e1) * t707;
t650 = ((m(3) * t715 + t501 * t553) * t588 + t591 * t501 * t556 / 0.2e1) * t706;
t640 = t529 * t683;
t639 = t530 * t682;
t638 = t531 * t681;
t637 = t529 * t652;
t636 = t530 * t651;
t635 = t531 * t650;
t634 = t713 * t561 * t708;
t633 = t714 * t562 * t707;
t632 = t715 * t563 * t706;
t631 = mrSges(3,1) * t712 + Ifges(3,2) + Ifges(2,3);
t627 = t529 * t634;
t626 = t530 * t633;
t625 = t531 * t632;
t621 = -t560 * t566 + 0.2e1 * t644;
t618 = -t559 * t565 + 0.2e1 * t645;
t615 = -t558 * t564 + 0.2e1 * t646;
t597 = mrSges(4,1);
t596 = mrSges(4,2);
t546 = m(3) * t657 + mrSges(3,3) * t711 + t631;
t545 = m(3) * t658 + mrSges(3,3) * t710 + t631;
t544 = m(3) * t659 + mrSges(3,3) * t709 + t631;
t543 = t560 * t569 + 0.2e1 * t649;
t542 = t559 * t568 + 0.2e1 * t648;
t541 = t558 * t567 + 0.2e1 * t647;
t516 = -t543 * t660 + t608 * t566 + t621 * t577 + t566 - t644;
t515 = t543 * t577 - t569 * t608 + t621 * t660 - t569 - t649;
t514 = -t542 * t661 + t608 * t565 + t618 * t576 + t565 - t645;
t513 = t542 * t576 - t568 * t608 + t618 * t661 - t568 - t648;
t512 = -t541 * t662 + t608 * t564 + t615 * t575 + t564 - t646;
t511 = t541 * t575 - t567 * t608 + t615 * t662 - t567 - t647;
t504 = (t527 * t552 - t528 * t549) * t531;
t503 = (t525 * t551 - t526 * t548) * t530;
t502 = (t523 * t550 - t524 * t547) * t529;
t486 = (-t519 * t549 + t522 * t552) * t531;
t485 = (-t518 * t548 + t521 * t551) * t530;
t484 = (-t517 * t547 + t520 * t550) * t529;
t480 = (t515 * t552 - t516 * t549) * t531;
t479 = (t513 * t551 - t514 * t548) * t530;
t478 = (t511 * t550 - t512 * t547) * t529;
t477 = (-t528 * t570 + (-t516 * t591 + t519) * m(3)) * t531;
t476 = (-t527 * t570 + (-t515 * t591 + t522) * m(3)) * t531;
t475 = (-t526 * t570 + (-t514 * t590 + t518) * m(3)) * t530;
t474 = (-t525 * t570 + (-t513 * t590 + t521) * m(3)) * t530;
t473 = (-t524 * t570 + (-t512 * t589 + t517) * m(3)) * t529;
t472 = (-t523 * t570 + (-t511 * t589 + t520) * m(3)) * t529;
t471 = (t516 * t573 - t519 * t702 + t528 * t540) * t531;
t470 = (t515 * t573 - t522 * t702 + t527 * t540) * t531;
t469 = (t514 * t573 - t518 * t703 + t526 * t539) * t530;
t468 = (t513 * t573 - t521 * t703 + t525 * t539) * t530;
t467 = (t512 * t573 - t517 * t704 + t524 * t538) * t529;
t466 = (t511 * t573 - t520 * t704 + t523 * t538) * t529;
t465 = (t516 * t540 - t519 * t570 + t528 * t546) * t531;
t464 = (t515 * t540 - t522 * t570 + t527 * t546) * t531;
t463 = (t514 * t539 - t518 * t570 + t526 * t545) * t530;
t462 = (t513 * t539 - t521 * t570 + t525 * t545) * t530;
t461 = (t512 * t538 - t517 * t570 + t524 * t544) * t529;
t460 = (t511 * t538 - t520 * t570 + t523 * t544) * t529;
t453 = -t504 * t570 + (-t480 * t591 + t486) * m(3);
t452 = -t503 * t570 + (-t479 * t590 + t485) * m(3);
t451 = -t502 * t570 + (-t478 * t589 + t484) * m(3);
t450 = t480 * t573 - t486 * t702 + t504 * t540;
t449 = t479 * t573 - t485 * t703 + t503 * t539;
t448 = t478 * t573 - t484 * t704 + t502 * t538;
t447 = t480 * t540 - t486 * t570 + t504 * t546;
t446 = t479 * t539 - t485 * t570 + t503 * t545;
t445 = t478 * t538 - t484 * t570 + t502 * t544;
t429 = -t438 * t540 + t441 * t570 - t442 * t546;
t428 = -t437 * t539 + t440 * t570 - t444 * t545;
t427 = -t436 * t538 + t439 * t570 - t443 * t544;
t1 = [(-(t465 * t528 + t471 * t516 + t477 * t519) * t552 - (t465 * t527 + t471 * t515 + t477 * t522) * t549) * t663 + t516 * t687 + t429 * t666 + t519 * t684 + t516 * t635 + t528 * t625 - t519 * t638 + (-(t463 * t526 + t469 * t514 + t475 * t518) * t551 - (t463 * t525 + t469 * t513 + t475 * t521) * t548) * t664 + t514 * t688 + t428 * t668 + t518 * t685 + t514 * t636 + t526 * t626 - t518 * t639 + (-(t461 * t524 + t467 * t512 + t473 * t517) * t550 - (t461 * t523 + t467 * t511 + t473 * t520) * t547) * t665 + t512 * t689 + t427 * t670 + t517 * t686 + t512 * t637 + t524 * t627 - t517 * t640 - t578 * (-t571 * t596 + t572 * t597); (-(t464 * t528 + t470 * t516 + t476 * t519) * t552 - (t464 * t527 + t470 * t515 + t476 * t522) * t549) * t663 + t515 * t687 + t429 * t667 + t522 * t684 + t515 * t635 + t527 * t625 - t522 * t638 + (-(t462 * t526 + t468 * t514 + t474 * t518) * t551 - (t462 * t525 + t468 * t513 + t474 * t521) * t548) * t664 + t513 * t688 + t428 * t669 + t521 * t685 + t513 * t636 + t525 * t626 - t521 * t639 + (-(t460 * t524 + t466 * t512 + t472 * t517) * t550 - (t460 * t523 + t466 * t511 + t472 * t520) * t547) * t665 + t511 * t689 + t427 * t671 + t520 * t686 + t511 * t637 + t523 * t627 - t520 * t640 - t578 * (t571 * t597 + t572 * t596); (-(t447 * t528 + t450 * t516 + t453 * t519) * t552 - (t447 * t527 + t450 * t515 + t453 * t522) * t549) * t663 + (-(t446 * t526 + t449 * t514 + t452 * t518) * t551 - (t446 * t525 + t449 * t513 + t452 * t521) * t548) * t664 + (-(t445 * t524 + t448 * t512 + t451 * t517) * t550 - (t445 * t523 + t448 * t511 + t451 * t520) * t547) * t665 + (t429 + t632) * t504 + (t428 + t633) * t503 + (t427 + t634) * t502 + (t435 - t681) * t486 + (t434 - t682) * t485 + (t433 - t683) * t484 + (t432 + t650) * t480 + (t431 + t651) * t479 + (t430 + t652) * t478;];
taucX  = t1;
