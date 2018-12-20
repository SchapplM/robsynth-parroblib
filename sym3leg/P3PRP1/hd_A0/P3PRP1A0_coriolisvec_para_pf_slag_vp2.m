% Calculate vector of centrifugal and coriolis load on the joints for
% P3PRP1A0
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
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taucX = P3PRP1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:35:00
% EndTime: 2018-12-20 17:35:05
% DurationCPUTime: 4.35s
% Computational Cost: add. (26233->373), mult. (49823->588), div. (2250->3), fcn. (23653->14), ass. (0->258)
t606 = (pkin(2) ^ 2);
t577 = 1 + t606;
t598 = (qJ(3,1) ^ 2);
t558 = -t598 + t577;
t589 = cos(qJ(2,1));
t575 = t589 ^ 2;
t586 = sin(qJ(2,1));
t658 = t586 * t589;
t645 = pkin(2) * t658;
t655 = t598 + t606;
t706 = 2 * qJ(3,1);
t529 = 0.1e1 / (t558 * t575 + t645 * t706 - t655 - 0.1e1);
t580 = legFrame(1,3);
t564 = sin(t580);
t567 = cos(t580);
t696 = qJ(3,1) * t564;
t644 = pkin(2) * t696;
t617 = t567 * t598 - t644;
t695 = qJ(3,1) * t567;
t639 = pkin(2) * t695;
t618 = -t564 * t598 - t639;
t520 = t617 * t589 - t586 * (t564 - t618);
t593 = xP(3);
t569 = sin(t593);
t570 = cos(t593);
t601 = koppelP(1,2);
t604 = koppelP(1,1);
t550 = -t569 * t601 + t570 * t604;
t590 = xDP(3);
t591 = xDP(2);
t532 = t550 * t590 + t591;
t673 = t520 * t532;
t517 = t618 * t589 + t586 * (-t567 - t617);
t547 = t569 * t604 + t570 * t601;
t592 = xDP(1);
t535 = -t547 * t590 + t592;
t676 = t517 * t535;
t710 = (t673 + t676) * t529;
t597 = (qJ(3,2) ^ 2);
t557 = -t597 + t577;
t588 = cos(qJ(2,2));
t574 = t588 ^ 2;
t585 = sin(qJ(2,2));
t659 = t585 * t588;
t646 = pkin(2) * t659;
t656 = t597 + t606;
t705 = 2 * qJ(3,2);
t528 = 0.1e1 / (t557 * t574 + t646 * t705 - t656 - 0.1e1);
t579 = legFrame(2,3);
t563 = sin(t579);
t566 = cos(t579);
t692 = qJ(3,2) * t563;
t643 = pkin(2) * t692;
t614 = t566 * t597 - t643;
t691 = qJ(3,2) * t566;
t640 = pkin(2) * t691;
t615 = -t563 * t597 - t640;
t519 = t614 * t588 - t585 * (t563 - t615);
t600 = koppelP(2,2);
t603 = koppelP(2,1);
t549 = -t569 * t600 + t570 * t603;
t531 = t549 * t590 + t591;
t674 = t519 * t531;
t516 = t615 * t588 + t585 * (-t566 - t614);
t546 = t569 * t603 + t570 * t600;
t534 = -t546 * t590 + t592;
t677 = t516 * t534;
t709 = (t674 + t677) * t528;
t596 = (qJ(3,3) ^ 2);
t556 = -t596 + t577;
t587 = cos(qJ(2,3));
t573 = t587 ^ 2;
t584 = sin(qJ(2,3));
t660 = t584 * t587;
t647 = pkin(2) * t660;
t657 = t596 + t606;
t704 = 2 * qJ(3,3);
t527 = 0.1e1 / (t556 * t573 + t647 * t704 - t657 - 0.1e1);
t578 = legFrame(3,3);
t562 = sin(t578);
t565 = cos(t578);
t688 = qJ(3,3) * t562;
t642 = pkin(2) * t688;
t611 = t565 * t596 - t642;
t687 = qJ(3,3) * t565;
t641 = pkin(2) * t687;
t612 = -t562 * t596 - t641;
t518 = t611 * t587 - t584 * (t562 - t612);
t599 = koppelP(3,2);
t602 = koppelP(3,1);
t548 = -t569 * t599 + t570 * t602;
t530 = t548 * t590 + t591;
t675 = t518 * t530;
t515 = t612 * t587 + t584 * (-t565 - t611);
t545 = t569 * t602 + t570 * t599;
t533 = -t545 * t590 + t592;
t678 = t515 * t533;
t708 = (t675 + t678) * t527;
t707 = 0.2e1 * pkin(2);
t576 = t590 ^ 2;
t685 = qJ(3,3) * t587;
t651 = -0.2e1 * t685;
t521 = t562 * t651 + t584 * (pkin(2) * t562 - t687);
t522 = t565 * t651 + t584 * (pkin(2) * t565 + t688);
t497 = (t521 * t530 + t522 * t533) * t527;
t703 = 0.2e1 * t497;
t689 = qJ(3,2) * t588;
t652 = -0.2e1 * t689;
t523 = t563 * t652 + t585 * (pkin(2) * t563 - t691);
t524 = t566 * t652 + t585 * (pkin(2) * t566 + t692);
t498 = (t523 * t531 + t524 * t534) * t528;
t702 = 0.2e1 * t498;
t693 = qJ(3,1) * t589;
t653 = -0.2e1 * t693;
t525 = t564 * t653 + t586 * (pkin(2) * t564 - t695);
t526 = t567 * t653 + t586 * (pkin(2) * t567 + t696);
t499 = (t525 * t532 + t526 * t535) * t529;
t701 = 0.2e1 * t499;
t700 = m(3) * pkin(2);
t699 = m(3) * t587;
t698 = m(3) * t588;
t697 = m(3) * t589;
t559 = m(3) * qJ(3,3) + mrSges(3,3);
t560 = m(3) * qJ(3,2) + mrSges(3,3);
t561 = m(3) * qJ(3,1) + mrSges(3,3);
t694 = qJ(3,1) * t575;
t690 = qJ(3,2) * t574;
t686 = qJ(3,3) * t573;
t684 = t497 ^ 2 * t559;
t683 = t498 ^ 2 * t560;
t682 = t499 ^ 2 * t561;
t681 = t497 * t527;
t680 = t498 * t528;
t679 = t499 * t529;
t605 = pkin(2) * t606;
t485 = t596 * t497;
t628 = -t606 * t497 + t707 * t708 - t485;
t654 = -t606 - t577;
t434 = (t628 * t685 + ((t605 + (1 + t596) * pkin(2)) * t497 - t708 + t654 * t708) * t584) * t681;
t492 = pkin(2) * t497;
t456 = t492 - t708;
t437 = ((t456 - t708) * t660 + (-t497 * t573 + t497 - t628) * qJ(3,3)) * t681;
t441 = ((0.2e1 * (t492 + (-t678 / 0.2e1 - t675 / 0.2e1) * t527) * t686 - (pkin(2) * t456 - t485) * t660 - qJ(3,3) * t708) * t527 + (-qJ(3,3) + t647 - t686) * t527 * t708) * t497;
t552 = -mrSges(2,2) + t559;
t568 = mrSges(3,1) + t700;
t555 = mrSges(2,1) + t568;
t536 = t552 * t584 + t555 * t587;
t629 = mrSges(3,1) * t707 + Ifges(3,2) + Ifges(2,3);
t542 = m(3) * t657 + mrSges(3,3) * t704 + t629;
t425 = -t434 * t536 + t437 * t568 - t441 * t542;
t672 = t527 * t425;
t571 = m(1) + m(2) + m(3);
t428 = -t434 * t571 + t437 * t699 - t441 * t536;
t671 = t527 * t428;
t431 = t441 * t568 + (t434 * t587 - t437) * m(3);
t670 = t527 * t431;
t669 = t527 * t576;
t486 = t597 * t498;
t627 = -t606 * t498 + t707 * t709 - t486;
t435 = (t627 * t689 + ((t605 + (1 + t597) * pkin(2)) * t498 - t709 + t654 * t709) * t585) * t680;
t493 = pkin(2) * t498;
t457 = t493 - t709;
t438 = ((t457 - t709) * t659 + (-t498 * t574 + t498 - t627) * qJ(3,2)) * t680;
t442 = ((0.2e1 * (t493 + (-t677 / 0.2e1 - t674 / 0.2e1) * t528) * t690 - (pkin(2) * t457 - t486) * t659 - qJ(3,2) * t709) * t528 + (-qJ(3,2) + t646 - t690) * t528 * t709) * t498;
t553 = -mrSges(2,2) + t560;
t537 = t553 * t585 + t555 * t588;
t543 = m(3) * t656 + mrSges(3,3) * t705 + t629;
t426 = -t435 * t537 + t438 * t568 - t442 * t543;
t668 = t528 * t426;
t429 = -t435 * t571 + t438 * t698 - t442 * t537;
t667 = t528 * t429;
t432 = t442 * t568 + (t435 * t588 - t438) * m(3);
t666 = t528 * t432;
t665 = t528 * t576;
t487 = t598 * t499;
t626 = -t606 * t499 + t707 * t710 - t487;
t436 = (t626 * t693 + ((t605 + (1 + t598) * pkin(2)) * t499 - t710 + t654 * t710) * t586) * t679;
t488 = t499 * pkin(2);
t455 = t488 - t710;
t439 = ((t455 - t710) * t658 + (-t499 * t575 + t499 - t626) * qJ(3,1)) * t679;
t440 = ((0.2e1 * (t488 + (-t676 / 0.2e1 - t673 / 0.2e1) * t529) * t694 - (pkin(2) * t455 - t487) * t658 - qJ(3,1) * t710) * t529 + (-qJ(3,1) + t645 - t694) * t529 * t710) * t499;
t554 = -mrSges(2,2) + t561;
t538 = t554 * t586 + t555 * t589;
t544 = m(3) * t655 + mrSges(3,3) * t706 + t629;
t427 = -t436 * t538 + t439 * t568 - t440 * t544;
t664 = t529 * t427;
t430 = -t436 * t571 + t439 * t697 - t440 * t538;
t663 = t529 * t430;
t433 = t440 * t568 + (t436 * t589 - t439) * m(3);
t662 = t529 * t433;
t661 = t529 * t576;
t551 = -t700 / 0.2e1 - mrSges(2,1) / 0.2e1 - mrSges(3,1) / 0.2e1;
t650 = ((m(3) * t708 + t497 * t551) * t584 + t587 * t497 * t552 / 0.2e1) * t703;
t649 = ((m(3) * t709 + t498 * t551) * t585 + t588 * t498 * t553 / 0.2e1) * t702;
t648 = ((m(3) * t710 + t499 * t551) * t586 + t589 * t499 * t554 / 0.2e1) * t701;
t638 = t527 * t684;
t637 = t528 * t683;
t636 = t529 * t682;
t635 = t708 * t559 * t703;
t634 = t709 * t560 * t702;
t633 = t710 * t561 * t701;
t632 = t527 * t650;
t631 = t528 * t649;
t630 = t529 * t648;
t625 = t527 * t635;
t624 = t528 * t634;
t623 = t529 * t633;
t619 = -t558 * t564 + 0.2e1 * t639;
t616 = -t557 * t563 + 0.2e1 * t640;
t613 = -t556 * t562 + 0.2e1 * t641;
t595 = mrSges(4,1);
t594 = mrSges(4,2);
t541 = t558 * t567 + 0.2e1 * t644;
t540 = t557 * t566 + 0.2e1 * t643;
t539 = t556 * t565 + 0.2e1 * t642;
t514 = -t541 * t658 + t606 * t564 + t575 * t619 + t564 - t639;
t513 = t541 * t575 - t567 * t606 + t619 * t658 - t567 - t644;
t512 = -t540 * t659 + t606 * t563 + t574 * t616 + t563 - t640;
t511 = t540 * t574 - t566 * t606 + t616 * t659 - t566 - t643;
t510 = -t539 * t660 + t606 * t562 + t573 * t613 + t562 - t641;
t509 = t539 * t573 - t565 * t606 + t613 * t660 - t565 - t642;
t502 = (t525 * t550 - t526 * t547) * t529;
t501 = (t523 * t549 - t524 * t546) * t528;
t500 = (t521 * t548 - t522 * t545) * t527;
t484 = (-t517 * t547 + t520 * t550) * t529;
t483 = (-t516 * t546 + t519 * t549) * t528;
t482 = (-t515 * t545 + t518 * t548) * t527;
t478 = (t513 * t550 - t514 * t547) * t529;
t477 = (t511 * t549 - t512 * t546) * t528;
t476 = (t509 * t548 - t510 * t545) * t527;
t475 = (-t526 * t568 + (-t514 * t589 + t517) * m(3)) * t529;
t474 = (-t525 * t568 + (-t513 * t589 + t520) * m(3)) * t529;
t473 = (-t524 * t568 + (-t512 * t588 + t516) * m(3)) * t528;
t472 = (-t523 * t568 + (-t511 * t588 + t519) * m(3)) * t528;
t471 = (-t522 * t568 + (-t510 * t587 + t515) * m(3)) * t527;
t470 = (-t521 * t568 + (-t509 * t587 + t518) * m(3)) * t527;
t469 = (t514 * t571 - t517 * t697 + t526 * t538) * t529;
t468 = (t513 * t571 - t520 * t697 + t525 * t538) * t529;
t467 = (t512 * t571 - t516 * t698 + t524 * t537) * t528;
t466 = (t511 * t571 - t519 * t698 + t523 * t537) * t528;
t465 = (t510 * t571 - t515 * t699 + t522 * t536) * t527;
t464 = (t509 * t571 - t518 * t699 + t521 * t536) * t527;
t463 = (t514 * t538 - t517 * t568 + t526 * t544) * t529;
t462 = (t513 * t538 - t520 * t568 + t525 * t544) * t529;
t461 = (t512 * t537 - t516 * t568 + t524 * t543) * t528;
t460 = (t511 * t537 - t519 * t568 + t523 * t543) * t528;
t459 = (t510 * t536 - t515 * t568 + t522 * t542) * t527;
t458 = (t509 * t536 - t518 * t568 + t521 * t542) * t527;
t451 = -t502 * t568 + (-t478 * t589 + t484) * m(3);
t450 = -t501 * t568 + (-t477 * t588 + t483) * m(3);
t449 = -t500 * t568 + (-t476 * t587 + t482) * m(3);
t448 = t478 * t571 - t484 * t697 + t502 * t538;
t447 = t477 * t571 - t483 * t698 + t501 * t537;
t446 = t476 * t571 - t482 * t699 + t500 * t536;
t445 = t478 * t538 - t484 * t568 + t502 * t544;
t444 = t477 * t537 - t483 * t568 + t501 * t543;
t443 = t476 * t536 - t482 * t568 + t500 * t542;
t1 = [(-(t463 * t526 + t469 * t514 + t475 * t517) * t550 - (t463 * t525 + t469 * t513 + t475 * t520) * t547) * t661 + t514 * t663 + t526 * t664 + t517 * t662 + t514 * t630 + t526 * t623 - t517 * t636 + (-(t461 * t524 + t467 * t512 + t473 * t516) * t549 - (t461 * t523 + t467 * t511 + t473 * t519) * t546) * t665 + t512 * t667 + t524 * t668 + t516 * t666 + t512 * t631 + t524 * t624 - t516 * t637 + (-(t459 * t522 + t465 * t510 + t471 * t515) * t548 - (t459 * t521 + t465 * t509 + t471 * t518) * t545) * t669 + t510 * t671 + t522 * t672 + t515 * t670 + t510 * t632 + t522 * t625 - t515 * t638 - t576 * (-t569 * t594 + t570 * t595); (-(t462 * t526 + t468 * t514 + t474 * t517) * t550 - (t462 * t525 + t468 * t513 + t474 * t520) * t547) * t661 + t513 * t663 + t525 * t664 + t520 * t662 + t513 * t630 + t525 * t623 - t520 * t636 + (-(t460 * t524 + t466 * t512 + t472 * t516) * t549 - (t460 * t523 + t466 * t511 + t472 * t519) * t546) * t665 + t511 * t667 + t523 * t668 + t519 * t666 + t511 * t631 + t523 * t624 - t519 * t637 + (-(t458 * t522 + t464 * t510 + t470 * t515) * t548 - (t458 * t521 + t464 * t509 + t470 * t518) * t545) * t669 + t509 * t671 + t521 * t672 + t518 * t670 + t509 * t632 + t521 * t625 - t518 * t638 - t576 * (t569 * t595 + t570 * t594); (-(t445 * t526 + t448 * t514 + t451 * t517) * t550 - (t445 * t525 + t448 * t513 + t451 * t520) * t547) * t661 + (-(t444 * t524 + t447 * t512 + t450 * t516) * t549 - (t444 * t523 + t447 * t511 + t450 * t519) * t546) * t665 + (-(t443 * t522 + t446 * t510 + t449 * t515) * t548 - (t443 * t521 + t446 * t509 + t449 * t518) * t545) * t669 + (t427 + t633) * t502 + (t426 + t634) * t501 + (t425 + t635) * t500 + (t433 - t682) * t484 + (t432 - t683) * t483 + (t431 - t684) * t482 + (t430 + t648) * t478 + (t429 + t649) * t477 + (t428 + t650) * t476;];
taucX  = t1;
