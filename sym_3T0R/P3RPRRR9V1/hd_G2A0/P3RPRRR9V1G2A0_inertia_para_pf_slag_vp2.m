% Calculate inertia matrix for parallel robot
% P3RPRRR9V1G2A0
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
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:53
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:51:38
% EndTime: 2020-08-06 18:51:39
% DurationCPUTime: 1.38s
% Computational Cost: add. (4080->289), mult. (5952->508), div. (396->10), fcn. (4008->26), ass. (0->222)
t653 = 2 * pkin(1);
t652 = 2 * pkin(3);
t651 = 4 * Ifges(3,4);
t530 = cos(pkin(7));
t518 = t530 ^ 2;
t650 = 0.2e1 * t518;
t649 = 2 * mrSges(2,3) + 2 * mrSges(3,3);
t648 = (pkin(2) * mrSges(3,1));
t529 = sin(pkin(7));
t498 = pkin(1) * t529;
t551 = (m(2) + m(3));
t647 = pkin(1) * t551;
t545 = cos(qJ(3,3));
t526 = t545 ^ 2;
t646 = t526 * pkin(3);
t547 = cos(qJ(3,2));
t527 = t547 ^ 2;
t645 = t527 * pkin(3);
t549 = cos(qJ(3,1));
t528 = t549 ^ 2;
t644 = t528 * pkin(3);
t643 = t545 * pkin(2);
t642 = t547 * pkin(2);
t641 = t549 * pkin(2);
t534 = Ifges(3,1) - Ifges(3,2);
t640 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t533 = pkin(5) + qJ(2,1);
t532 = pkin(5) + qJ(2,2);
t531 = pkin(5) + qJ(2,3);
t539 = sin(qJ(3,3));
t639 = Ifges(3,4) * t539;
t541 = sin(qJ(3,2));
t638 = Ifges(3,4) * t541;
t543 = sin(qJ(3,1));
t637 = Ifges(3,4) * t543;
t636 = m(3) * pkin(2) + mrSges(2,1);
t635 = t539 * mrSges(3,1) + mrSges(2,2);
t634 = t541 * mrSges(3,1) + mrSges(2,2);
t633 = t543 * mrSges(3,1) + mrSges(2,2);
t540 = sin(qJ(1,3));
t515 = pkin(6) + t531;
t546 = cos(qJ(1,3));
t567 = pkin(1) * t540 - t546 * t515;
t587 = t529 * t539;
t452 = t567 * t587 + (t526 - 0.1e1) * t540 * pkin(3);
t458 = pkin(1) * t539 + (-pkin(3) + t643 + 0.2e1 * t646) * t529;
t554 = -pkin(3) / 0.2e1;
t473 = t646 + t643 / 0.2e1 + t554;
t536 = legFrame(3,2);
t505 = sin(t536);
t508 = cos(t536);
t570 = t540 * t587;
t559 = pkin(2) * t570 + (t570 * t652 - t567) * t545;
t581 = t545 * (-t539 * pkin(3) + t498);
t599 = t505 * t540;
t555 = pkin(2) / 0.2e1;
t608 = (t545 * pkin(3) + t555) * t539;
t428 = (-t473 * t599 + t508 * t608) * t650 + (t508 * t458 + t559 * t505) * t530 + t452 * t505 + t508 * t581;
t470 = 0.1e1 / (t530 * t545 - t587);
t632 = t428 * t470;
t593 = t508 * t540;
t429 = (t473 * t593 + t505 * t608) * t650 + (t505 * t458 - t559 * t508) * t530 - t452 * t508 + t505 * t581;
t631 = t429 * t470;
t542 = sin(qJ(1,2));
t516 = pkin(6) + t532;
t548 = cos(qJ(1,2));
t566 = pkin(1) * t542 - t548 * t516;
t586 = t529 * t541;
t453 = t566 * t586 + (t527 - 0.1e1) * t542 * pkin(3);
t459 = pkin(1) * t541 + (-pkin(3) + t642 + 0.2e1 * t645) * t529;
t474 = t645 + t642 / 0.2e1 + t554;
t537 = legFrame(2,2);
t506 = sin(t537);
t509 = cos(t537);
t569 = t542 * t586;
t558 = pkin(2) * t569 + (t569 * t652 - t566) * t547;
t580 = t547 * (-t541 * pkin(3) + t498);
t597 = t506 * t542;
t607 = (t547 * pkin(3) + t555) * t541;
t430 = (-t474 * t597 + t509 * t607) * t650 + (t509 * t459 + t558 * t506) * t530 + t453 * t506 + t509 * t580;
t471 = 0.1e1 / (t530 * t547 - t586);
t630 = t430 * t471;
t591 = t509 * t542;
t431 = (t474 * t591 + t506 * t607) * t650 + (t506 * t459 - t558 * t509) * t530 - t453 * t509 + t506 * t580;
t629 = t431 * t471;
t544 = sin(qJ(1,1));
t517 = pkin(6) + t533;
t550 = cos(qJ(1,1));
t565 = pkin(1) * t544 - t550 * t517;
t585 = t529 * t543;
t454 = t565 * t585 + (t528 - 0.1e1) * t544 * pkin(3);
t460 = pkin(1) * t543 + (-pkin(3) + t641 + 0.2e1 * t644) * t529;
t475 = t644 + t641 / 0.2e1 + t554;
t538 = legFrame(1,2);
t507 = sin(t538);
t510 = cos(t538);
t568 = t544 * t585;
t557 = pkin(2) * t568 + (t568 * t652 - t565) * t549;
t579 = t549 * (-t543 * pkin(3) + t498);
t595 = t507 * t544;
t606 = (t549 * pkin(3) + t555) * t543;
t432 = (-t475 * t595 + t510 * t606) * t650 + (t510 * t460 + t557 * t507) * t530 + t454 * t507 + t510 * t579;
t472 = 0.1e1 / (t530 * t549 - t585);
t628 = t432 * t472;
t589 = t510 * t544;
t433 = (t475 * t589 + t507 * t606) * t650 + (t507 * t460 - t557 * t510) * t530 - t454 * t510 + t507 * t579;
t627 = t433 * t472;
t564 = -mrSges(3,2) * t539 + t636;
t449 = (-mrSges(3,1) * t545 - t564) * t530 + (mrSges(3,2) * t545 + t635) * t529 - t647;
t488 = t530 * pkin(2) + pkin(1);
t519 = pkin(7) + qJ(3,3);
t495 = cos(t519);
t455 = t540 * t515 + (pkin(3) * t495 + t488) * t546;
t502 = 0.1e1 / t515;
t443 = (t449 * t546 + t455 * t551) * t502;
t626 = t443 * t470;
t563 = -mrSges(3,2) * t541 + t636;
t450 = (-mrSges(3,1) * t547 - t563) * t530 + (mrSges(3,2) * t547 + t634) * t529 - t647;
t520 = pkin(7) + qJ(3,2);
t496 = cos(t520);
t456 = t542 * t516 + (pkin(3) * t496 + t488) * t548;
t503 = 0.1e1 / t516;
t444 = (t450 * t548 + t456 * t551) * t503;
t625 = t444 * t471;
t562 = -mrSges(3,2) * t543 + t636;
t451 = (-mrSges(3,1) * t549 - t562) * t530 + (mrSges(3,2) * t549 + t633) * t529 - t647;
t521 = pkin(7) + qJ(3,1);
t497 = cos(t521);
t457 = t544 * t517 + (pkin(3) * t497 + t488) * t550;
t504 = 0.1e1 / t517;
t445 = (t451 * t550 + t457 * t551) * t504;
t624 = t445 * t472;
t492 = sin(t519);
t461 = t508 * t492 - t495 * t599;
t489 = 0.1e1 / t495;
t623 = t461 * t489;
t622 = t461 * t502;
t462 = t505 * t492 + t495 * t593;
t621 = t462 * t489;
t620 = t462 * t502;
t493 = sin(t520);
t463 = t509 * t493 - t496 * t597;
t490 = 0.1e1 / t496;
t619 = t463 * t490;
t618 = t463 * t503;
t464 = t506 * t493 + t496 * t591;
t617 = t464 * t490;
t616 = t464 * t503;
t494 = sin(t521);
t465 = t510 * t494 - t497 * t595;
t491 = 0.1e1 / t497;
t615 = t465 * t491;
t614 = t465 * t504;
t466 = t507 * t494 + t497 * t589;
t613 = t466 * t491;
t612 = t466 * t504;
t611 = t470 * t551;
t610 = t471 * t551;
t609 = t472 * t551;
t605 = t489 * t505;
t604 = t489 * t508;
t603 = t490 * t506;
t602 = t490 * t509;
t601 = t491 * t507;
t600 = t491 * t510;
t556 = 1 / pkin(3);
t598 = t505 * t556;
t596 = t506 * t556;
t594 = t507 * t556;
t592 = t508 * t556;
t590 = t509 * t556;
t588 = t510 * t556;
t584 = t534 * t526;
t583 = t534 * t527;
t582 = t534 * t528;
t522 = -0.2e1 * pkin(2) * mrSges(3,2);
t578 = -0.2e1 * t498;
t577 = mrSges(3,2) * t498;
t479 = t531 * mrSges(3,2) - Ifges(3,6);
t482 = t531 * mrSges(3,1) - Ifges(3,5);
t446 = (-t479 * t545 - t539 * t482) * t530 + t529 * (t479 * t539 - t482 * t545);
t576 = t446 * t546 * t556;
t480 = t532 * mrSges(3,2) - Ifges(3,6);
t483 = t532 * mrSges(3,1) - Ifges(3,5);
t447 = (-t480 * t547 - t541 * t483) * t530 + t529 * (t480 * t541 - t483 * t547);
t575 = t447 * t548 * t556;
t481 = t533 * mrSges(3,2) - Ifges(3,6);
t484 = t533 * mrSges(3,1) - Ifges(3,5);
t448 = (-t481 * t549 - t543 * t484) * t530 + t529 * (t481 * t543 - t484 * t549);
t574 = t448 * t550 * t556;
t573 = t449 * t470 * t502;
t572 = t450 * t471 * t503;
t571 = t451 * t472 * t504;
t561 = pkin(2) ^ 2 * m(3) - Ifges(2,1) + Ifges(2,2) + t534;
t560 = (t551 * pkin(1) ^ 2) + (2 * mrSges(3,3) * pkin(5)) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t525 = mrSges(3,1) * t653;
t523 = 2 * t648;
t442 = (Ifges(3,3) * t594 + t448 * t612) * t491;
t441 = (Ifges(3,3) * t588 + t448 * t614) * t491;
t440 = (Ifges(3,3) * t596 + t447 * t616) * t490;
t439 = (Ifges(3,3) * t590 + t447 * t618) * t490;
t438 = (Ifges(3,3) * t598 + t446 * t620) * t489;
t437 = (Ifges(3,3) * t592 + t446 * t622) * t489;
t436 = (-0.2e1 * t582 + (t523 + 0.4e1 * t637) * t549 + t543 * t522 + t561) * t518 + (t525 * t549 + t562 * t653 + (t528 * t651 + t522 * t549 + 0.2e1 * (t534 * t549 - t648) * t543 + t640) * t529) * t530 + t582 + 0.2e1 * (-t577 - t637) * t549 + t633 * t578 + t533 ^ 2 * m(3) + t560 + ((m(2) * qJ(2,1) + t649) * qJ(2,1));
t435 = (-0.2e1 * t583 + (t523 + 0.4e1 * t638) * t547 + t541 * t522 + t561) * t518 + (t525 * t547 + t563 * t653 + (t527 * t651 + t522 * t547 + 0.2e1 * (t534 * t547 - t648) * t541 + t640) * t529) * t530 + t583 + 0.2e1 * (-t577 - t638) * t547 + t634 * t578 + t532 ^ 2 * m(3) + t560 + ((m(2) * qJ(2,2) + t649) * qJ(2,2));
t434 = (-0.2e1 * t584 + (t523 + 0.4e1 * t639) * t545 + t539 * t522 + t561) * t518 + (t525 * t545 + t564 * t653 + (t526 * t651 + t522 * t545 + 0.2e1 * (t534 * t545 - t648) * t539 + t640) * t529) * t530 + t584 + 0.2e1 * (-t577 - t639) * t545 + t635 * t578 + t531 ^ 2 * m(3) + t560 + ((m(2) * qJ(2,3) + t649) * qJ(2,3));
t427 = (t436 * t550 + t451 * t457) * t504;
t426 = (t435 * t548 + t450 * t456) * t503;
t425 = (t434 * t546 + t449 * t455) * t502;
t424 = (t433 * t609 + t451 * t613) * t504;
t423 = (t432 * t609 + t451 * t615) * t504;
t422 = (t431 * t610 + t450 * t617) * t503;
t421 = (t430 * t610 + t450 * t619) * t503;
t420 = (t429 * t611 + t449 * t621) * t502;
t419 = (t428 * t611 + t449 * t623) * t502;
t418 = t433 * t571 + (t436 * t612 + t448 * t594) * t491;
t417 = t432 * t571 + (t436 * t614 + t448 * t588) * t491;
t416 = t431 * t572 + (t435 * t616 + t447 * t596) * t490;
t415 = t430 * t572 + (t435 * t618 + t447 * t590) * t490;
t414 = t429 * t573 + (t434 * t620 + t446 * t598) * t489;
t413 = t428 * t573 + (t434 * t622 + t446 * t592) * t489;
t1 = [m(4) + (t418 * t613 + t424 * t627) * t504 + (t416 * t617 + t422 * t629) * t503 + (t414 * t621 + t420 * t631) * t502 + (t438 * t605 + t440 * t603 + t442 * t601) * t556, (t418 * t615 + t424 * t628) * t504 + (t416 * t619 + t422 * t630) * t503 + (t414 * t623 + t420 * t632) * t502 + (t438 * t604 + t440 * t602 + t442 * t600) * t556, (t418 * t550 + t424 * t457) * t504 + (t416 * t548 + t422 * t456) * t503 + (t414 * t546 + t420 * t455) * t502; (t417 * t613 + t423 * t627) * t504 + (t415 * t617 + t421 * t629) * t503 + (t413 * t621 + t419 * t631) * t502 + (t437 * t605 + t439 * t603 + t441 * t601) * t556, m(4) + (t417 * t615 + t423 * t628) * t504 + (t415 * t619 + t421 * t630) * t503 + (t413 * t623 + t419 * t632) * t502 + (t437 * t604 + t439 * t602 + t441 * t600) * t556, (t417 * t550 + t423 * t457) * t504 + (t415 * t548 + t421 * t456) * t503 + (t413 * t546 + t419 * t455) * t502; (t433 * t624 + (t427 * t466 + t507 * t574) * t491) * t504 + (t431 * t625 + (t426 * t464 + t506 * t575) * t490) * t503 + (t429 * t626 + (t425 * t462 + t505 * t576) * t489) * t502, (t432 * t624 + (t427 * t465 + t510 * t574) * t491) * t504 + (t430 * t625 + (t426 * t463 + t509 * t575) * t490) * t503 + (t428 * t626 + (t425 * t461 + t508 * t576) * t489) * t502, m(4) + (t427 * t550 + t445 * t457) * t504 + (t426 * t548 + t444 * t456) * t503 + (t425 * t546 + t443 * t455) * t502;];
MX  = t1;
