% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPP1A0
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
%   pkin=[a2,a3,d1]';
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
% Datum: 2018-12-20 17:50
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function taucX = P3RPP1A0_coriolisvec_para_pf_slag_vp2(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1A0_coriolisvec_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:50:33
% EndTime: 2018-12-20 17:50:37
% DurationCPUTime: 4.14s
% Computational Cost: add. (27163->325), mult. (37037->531), div. (2067->3), fcn. (14716->14), ass. (0->246)
t615 = (qJ(3,3) ^ 2);
t628 = (pkin(1) ^ 2);
t690 = 1 + t628;
t693 = 2 * qJ(3,3);
t570 = pkin(1) * t693 + t615 + t690;
t601 = sin(qJ(1,3));
t597 = -qJ(3,3) - pkin(1);
t604 = cos(qJ(1,3));
t665 = t597 * t604;
t654 = qJ(2,3) * t665;
t540 = t570 * t601 + t654;
t666 = t597 * t601;
t573 = qJ(2,3) * t666;
t543 = t570 * t604 - t573;
t594 = legFrame(3,3);
t582 = sin(t594);
t585 = cos(t594);
t513 = t540 * t585 + t543 * t582;
t516 = -t540 * t582 + t543 * t585;
t612 = xP(3);
t588 = sin(t612);
t589 = cos(t612);
t621 = koppelP(3,2);
t624 = koppelP(3,1);
t555 = -t588 * t621 + t589 * t624;
t607 = xDP(3);
t608 = xDP(2);
t534 = t555 * t607 + t608;
t552 = t588 * t624 + t589 * t621;
t609 = xDP(1);
t537 = -t552 * t607 + t609;
t616 = qJ(2,3) ^ 2;
t704 = t570 + t616;
t564 = 0.1e1 / t704;
t477 = (t513 * t534 + t516 * t537) * t564;
t714 = t477 * t597;
t617 = (qJ(3,2) ^ 2);
t694 = 2 * qJ(3,2);
t571 = pkin(1) * t694 + t617 + t690;
t602 = sin(qJ(1,2));
t598 = -qJ(3,2) - pkin(1);
t605 = cos(qJ(1,2));
t663 = t598 * t605;
t655 = qJ(2,2) * t663;
t541 = t571 * t602 + t655;
t664 = t598 * t602;
t574 = qJ(2,2) * t664;
t544 = t571 * t605 - t574;
t595 = legFrame(2,3);
t583 = sin(t595);
t586 = cos(t595);
t514 = t541 * t586 + t544 * t583;
t517 = -t541 * t583 + t544 * t586;
t622 = koppelP(2,2);
t625 = koppelP(2,1);
t556 = -t588 * t622 + t589 * t625;
t535 = t556 * t607 + t608;
t553 = t588 * t625 + t589 * t622;
t538 = -t553 * t607 + t609;
t618 = qJ(2,2) ^ 2;
t703 = t571 + t618;
t565 = 0.1e1 / t703;
t478 = (t514 * t535 + t517 * t538) * t565;
t713 = t478 * t598;
t619 = (qJ(3,1) ^ 2);
t695 = 2 * qJ(3,1);
t572 = pkin(1) * t695 + t619 + t690;
t603 = sin(qJ(1,1));
t599 = -qJ(3,1) - pkin(1);
t606 = cos(qJ(1,1));
t661 = t599 * t606;
t656 = qJ(2,1) * t661;
t542 = t572 * t603 + t656;
t662 = t599 * t603;
t575 = qJ(2,1) * t662;
t545 = t572 * t606 - t575;
t596 = legFrame(1,3);
t584 = sin(t596);
t587 = cos(t596);
t515 = t542 * t587 + t545 * t584;
t518 = -t542 * t584 + t545 * t587;
t623 = koppelP(1,2);
t626 = koppelP(1,1);
t557 = -t588 * t623 + t589 * t626;
t536 = t557 * t607 + t608;
t554 = t588 * t626 + t589 * t623;
t539 = -t554 * t607 + t609;
t620 = qJ(2,1) ^ 2;
t702 = t572 + t620;
t566 = 0.1e1 / t702;
t479 = (t515 * t536 + t518 * t539) * t566;
t712 = t479 * t599;
t590 = 0.1e1 + t616;
t591 = 0.1e1 + t618;
t592 = 0.1e1 + t620;
t708 = -m(2) * pkin(1) + mrSges(2,2);
t560 = qJ(2,1) * t606 + t662;
t563 = qJ(2,1) * t603 - t661;
t530 = t560 * t587 - t563 * t584;
t533 = t560 * t584 + t563 * t587;
t707 = (t530 * t539 + t533 * t536) * t566;
t559 = qJ(2,2) * t605 + t664;
t562 = qJ(2,2) * t602 - t663;
t529 = t559 * t586 - t562 * t583;
t532 = t559 * t583 + t562 * t586;
t706 = (t529 * t538 + t532 * t535) * t565;
t558 = qJ(2,3) * t604 + t666;
t561 = qJ(2,3) * t601 - t665;
t528 = t558 * t585 - t561 * t582;
t531 = t558 * t582 + t561 * t585;
t705 = (t528 * t537 + t531 * t534) * t564;
t600 = mrSges(3,2) + mrSges(2,3);
t610 = m(2) + m(3);
t701 = qJ(2,1) * t610 + t600;
t700 = qJ(2,2) * t610 + t600;
t699 = qJ(2,3) * t610 + t600;
t698 = -2 * m(3);
t697 = 2 * m(3);
t696 = 2 * pkin(1);
t593 = t607 ^ 2;
t692 = 0.2e1 * t600;
t691 = -3 * t628;
t689 = -mrSges(2,2) + mrSges(3,3);
t548 = t592 * t603 - t656;
t551 = -t592 * t606 - t575;
t521 = t548 * t587 - t551 * t584;
t524 = t548 * t584 + t551 * t587;
t485 = (t521 * t539 + t524 * t536) * t566;
t687 = qJ(2,1) * t485;
t547 = t591 * t602 - t655;
t550 = -t591 * t605 - t574;
t520 = t547 * t586 - t550 * t583;
t523 = t547 * t583 + t550 * t586;
t484 = (t520 * t538 + t523 * t535) * t565;
t685 = qJ(2,2) * t484;
t546 = t590 * t601 - t654;
t549 = -t590 * t604 - t573;
t522 = t546 * t582 + t549 * t585;
t676 = t522 * t564;
t519 = t546 * t585 - t549 * t582;
t677 = t519 * t564;
t483 = t534 * t676 + t537 * t677;
t683 = qJ(2,3) * t483;
t627 = pkin(1) * t628;
t680 = t705 * t564;
t441 = (0.2e1 * t570 * t483 + 0.2e1 * qJ(2,3) * t714 + (-t627 + (-t615 - t616 + t691) * qJ(3,3) - qJ(3,3) + (-(3 * t615) - t590) * pkin(1)) * t705) * t680;
t447 = (-0.2e1 * t597 * t683 - t477 + (-t616 - t590) * t477 - t704 * t705 * qJ(2,3)) * t680;
t450 = 0.2e1 * (t683 - t714) * t680;
t657 = Ifges(2,1) + Ifges(3,1) + Ifges(1,3);
t525 = (m(3) * t615) + (mrSges(3,3) * t693) + (t628 + t616) * t610 + ((m(3) * qJ(3,3) + t689) * t696) + qJ(2,3) * t692 + t657;
t576 = -m(3) * t597 + mrSges(3,3);
t567 = -t576 + t708;
t579 = m(3) * qJ(2,3) + mrSges(3,2);
t435 = -t441 * t579 - t447 * t567 - t450 * t525;
t682 = t435 * t564;
t438 = -m(3) * t441 - t450 * t579;
t681 = t438 * t564;
t679 = t706 * t565;
t678 = t707 * t566;
t675 = t564 * t593;
t442 = (0.2e1 * t571 * t484 + 0.2e1 * qJ(2,2) * t713 + (-t627 + (-t617 - t618 + t691) * qJ(3,2) - qJ(3,2) + (-(3 * t617) - t591) * pkin(1)) * t706) * t679;
t448 = (-0.2e1 * t598 * t685 - t478 + (-t618 - t591) * t478 - t703 * t706 * qJ(2,2)) * t679;
t451 = 0.2e1 * (t685 - t713) * t679;
t526 = (m(3) * t617) + (mrSges(3,3) * t694) + (t628 + t618) * t610 + ((m(3) * qJ(3,2) + t689) * t696) + qJ(2,2) * t692 + t657;
t577 = -m(3) * t598 + mrSges(3,3);
t568 = -t577 + t708;
t580 = m(3) * qJ(2,2) + mrSges(3,2);
t436 = -t442 * t580 - t448 * t568 - t451 * t526;
t674 = t565 * t436;
t439 = -m(3) * t442 - t451 * t580;
t673 = t565 * t439;
t445 = -t448 * t610 - t451 * t568;
t672 = t565 * t445;
t671 = t565 * t593;
t443 = (0.2e1 * t572 * t485 + 0.2e1 * qJ(2,1) * t712 + (-t627 + (-t619 - t620 + t691) * qJ(3,1) - qJ(3,1) + (-(3 * t619) - t592) * pkin(1)) * t707) * t678;
t449 = (-0.2e1 * t599 * t687 - t479 + (-t620 - t592) * t479 - t702 * t707 * qJ(2,1)) * t678;
t452 = 0.2e1 * (t687 - t712) * t678;
t527 = (m(3) * t619) + (mrSges(3,3) * t695) + (t628 + t620) * t610 + ((m(3) * qJ(3,1) + t689) * t696) + qJ(2,1) * t692 + t657;
t578 = -m(3) * t599 + mrSges(3,3);
t569 = -t578 + t708;
t581 = m(3) * qJ(2,1) + mrSges(3,2);
t437 = -t443 * t581 - t449 * t569 - t452 * t527;
t670 = t566 * t437;
t440 = -m(3) * t443 - t452 * t581;
t669 = t566 * t440;
t446 = -t449 * t610 - t452 * t569;
t668 = t566 * t446;
t667 = t566 * t593;
t660 = 0.2e1 * t705 * (t576 * t477 + t483 * t699);
t659 = 0.2e1 * t706 * (t577 * t478 + t484 * t700);
t658 = 0.2e1 * t707 * (t578 * t479 + t485 * t701);
t459 = t477 * t697 + t699 * t705;
t653 = t459 * t680;
t462 = t483 * t698 + t576 * t705;
t652 = t462 * t680;
t463 = t484 * t698 + t577 * t706;
t651 = t463 * t679;
t460 = t478 * t697 + t700 * t706;
t650 = t460 * t679;
t464 = t485 * t698 + t578 * t707;
t649 = t464 * t678;
t461 = t479 * t697 + t701 * t707;
t648 = t461 * t678;
t647 = t564 * t660;
t646 = t565 * t659;
t645 = t566 * t658;
t614 = mrSges(4,1);
t613 = mrSges(4,2);
t506 = (t524 * t610 + t533 * t569) * t566;
t505 = (t523 * t610 + t532 * t568) * t565;
t504 = (t522 * t610 + t531 * t567) * t564;
t503 = (t521 * t610 + t530 * t569) * t566;
t502 = (t520 * t610 + t529 * t568) * t565;
t501 = (t519 * t610 + t528 * t567) * t564;
t500 = (m(3) * t515 + t533 * t581) * t566;
t499 = (m(3) * t514 + t532 * t580) * t565;
t498 = (m(3) * t513 + t531 * t579) * t564;
t497 = (m(3) * t518 + t530 * t581) * t566;
t496 = (m(3) * t517 + t529 * t580) * t565;
t495 = (m(3) * t516 + t528 * t579) * t564;
t494 = (-t530 * t554 + t533 * t557) * t566;
t493 = (-t529 * t553 + t532 * t556) * t565;
t492 = (-t528 * t552 + t531 * t555) * t564;
t488 = (-t521 * t554 + t524 * t557) * t566;
t487 = (-t520 * t553 + t523 * t556) * t565;
t486 = (-t519 * t552 + t522 * t555) * t564;
t482 = (t515 * t557 - t518 * t554) * t566;
t481 = (t514 * t556 - t517 * t553) * t565;
t480 = (t513 * t555 - t516 * t552) * t564;
t476 = (t515 * t581 + t524 * t569 + t527 * t533) * t566;
t475 = (t514 * t580 + t523 * t568 + t526 * t532) * t565;
t474 = (t513 * t579 + t522 * t567 + t525 * t531) * t564;
t473 = (t518 * t581 + t521 * t569 + t527 * t530) * t566;
t472 = (t517 * t580 + t520 * t568 + t526 * t529) * t565;
t471 = (t516 * t579 + t519 * t567 + t525 * t528) * t564;
t470 = t488 * t610 + t494 * t569;
t469 = t487 * t610 + t493 * t568;
t468 = t486 * t610 + t492 * t567;
t467 = m(3) * t482 + t494 * t581;
t466 = m(3) * t481 + t493 * t580;
t465 = m(3) * t480 + t492 * t579;
t455 = t482 * t581 + t488 * t569 + t494 * t527;
t454 = t481 * t580 + t487 * t568 + t493 * t526;
t453 = t480 * t579 + t486 * t567 + t492 * t525;
t444 = -t447 * t610 - t450 * t567;
t1 = [(-(t473 * t530 + t497 * t518 + t503 * t521) * t557 - (t473 * t533 + t497 * t515 + t503 * t524) * t554) * t667 + t530 * t670 + t521 * t668 + t518 * t669 + t530 * t645 - t521 * t648 - t518 * t649 + (-(t472 * t529 + t496 * t517 + t502 * t520) * t556 - (t472 * t532 + t496 * t514 + t502 * t523) * t553) * t671 + t529 * t674 + t520 * t672 + t517 * t673 + t529 * t646 - t520 * t650 - t517 * t651 + (-(t471 * t528 + t495 * t516 + t501 * t519) * t555 - (t471 * t531 + t495 * t513 + t501 * t522) * t552) * t675 + t528 * t682 + t444 * t677 + t516 * t681 + t528 * t647 - t519 * t653 - t516 * t652 - t593 * (-t588 * t613 + t589 * t614); (-(t476 * t530 + t500 * t518 + t506 * t521) * t557 - (t476 * t533 + t500 * t515 + t506 * t524) * t554) * t667 + t533 * t670 + t524 * t668 + t515 * t669 + t533 * t645 - t524 * t648 - t515 * t649 + (-(t475 * t529 + t499 * t517 + t505 * t520) * t556 - (t475 * t532 + t499 * t514 + t505 * t523) * t553) * t671 + t532 * t674 + t523 * t672 + t514 * t673 + t532 * t646 - t523 * t650 - t514 * t651 + (-(t474 * t528 + t498 * t516 + t504 * t519) * t555 - (t474 * t531 + t498 * t513 + t504 * t522) * t552) * t675 + t531 * t682 + t444 * t676 + t513 * t681 + t531 * t647 - t522 * t653 - t513 * t652 - t593 * (t588 * t614 + t589 * t613); (-(t455 * t530 + t467 * t518 + t470 * t521) * t557 - (t455 * t533 + t467 * t515 + t470 * t524) * t554) * t667 + (-(t454 * t529 + t466 * t517 + t469 * t520) * t556 - (t454 * t532 + t466 * t514 + t469 * t523) * t553) * t671 + (-(t453 * t528 + t465 * t516 + t468 * t519) * t555 - (t453 * t531 + t465 * t513 + t468 * t522) * t552) * t675 + (t437 + t658) * t494 + (t436 + t659) * t493 + (t435 + t660) * t492 + (-t461 * t707 + t446) * t488 + (-t460 * t706 + t445) * t487 + (-t459 * t705 + t444) * t486 + (-t464 * t707 + t440) * t482 + (-t463 * t706 + t439) * t481 + (-t462 * t705 + t438) * t480;];
taucX  = t1;
