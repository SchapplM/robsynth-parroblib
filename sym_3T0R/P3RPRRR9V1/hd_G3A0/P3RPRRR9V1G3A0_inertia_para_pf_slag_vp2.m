% Calculate inertia matrix for parallel robot
% P3RPRRR9V1G3A0
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
% Datum: 2020-08-06 18:58
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:56:36
% EndTime: 2020-08-06 18:56:37
% DurationCPUTime: 1.34s
% Computational Cost: add. (4080->289), mult. (5952->508), div. (396->10), fcn. (4008->26), ass. (0->222)
t659 = 2 * pkin(1);
t658 = 2 * pkin(3);
t657 = 4 * Ifges(3,4);
t536 = cos(pkin(7));
t524 = t536 ^ 2;
t656 = 0.2e1 * t524;
t655 = 2 * mrSges(2,3) + 2 * mrSges(3,3);
t654 = (pkin(2) * mrSges(3,1));
t535 = sin(pkin(7));
t504 = pkin(1) * t535;
t557 = (m(2) + m(3));
t653 = pkin(1) * t557;
t551 = cos(qJ(3,3));
t532 = t551 ^ 2;
t652 = t532 * pkin(3);
t553 = cos(qJ(3,2));
t533 = t553 ^ 2;
t651 = t533 * pkin(3);
t555 = cos(qJ(3,1));
t534 = t555 ^ 2;
t650 = t534 * pkin(3);
t649 = t551 * pkin(2);
t648 = t553 * pkin(2);
t647 = t555 * pkin(2);
t540 = Ifges(3,1) - Ifges(3,2);
t646 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t539 = pkin(5) + qJ(2,1);
t538 = pkin(5) + qJ(2,2);
t537 = pkin(5) + qJ(2,3);
t545 = sin(qJ(3,3));
t645 = Ifges(3,4) * t545;
t547 = sin(qJ(3,2));
t644 = Ifges(3,4) * t547;
t549 = sin(qJ(3,1));
t643 = Ifges(3,4) * t549;
t642 = m(3) * pkin(2) + mrSges(2,1);
t641 = t545 * mrSges(3,1) + mrSges(2,2);
t640 = t547 * mrSges(3,1) + mrSges(2,2);
t639 = t549 * mrSges(3,1) + mrSges(2,2);
t552 = cos(qJ(1,3));
t521 = pkin(6) + t537;
t546 = sin(qJ(1,3));
t584 = pkin(1) * t552 + t546 * t521;
t593 = t535 * t545;
t455 = t584 * t593 + (t532 - 0.1e1) * t552 * pkin(3);
t461 = pkin(1) * t545 + (-pkin(3) + t649 + 0.2e1 * t652) * t535;
t560 = -pkin(3) / 0.2e1;
t476 = t652 + t649 / 0.2e1 + t560;
t542 = legFrame(3,2);
t511 = sin(t542);
t514 = cos(t542);
t573 = t552 * t593;
t565 = pkin(2) * t573 + (t573 * t658 - t584) * t551;
t587 = t551 * (-t545 * pkin(3) + t504);
t605 = t511 * t552;
t561 = pkin(2) / 0.2e1;
t614 = (t551 * pkin(3) + t561) * t545;
t431 = (-t476 * t605 + t514 * t614) * t656 + (t514 * t461 + t565 * t511) * t536 + t455 * t511 + t514 * t587;
t473 = 0.1e1 / (t536 * t551 - t593);
t638 = t431 * t473;
t599 = t514 * t552;
t432 = (t476 * t599 + t511 * t614) * t656 + (t511 * t461 - t565 * t514) * t536 - t455 * t514 + t511 * t587;
t637 = t432 * t473;
t554 = cos(qJ(1,2));
t522 = pkin(6) + t538;
t548 = sin(qJ(1,2));
t583 = pkin(1) * t554 + t548 * t522;
t592 = t535 * t547;
t456 = t583 * t592 + (t533 - 0.1e1) * t554 * pkin(3);
t462 = pkin(1) * t547 + (-pkin(3) + t648 + 0.2e1 * t651) * t535;
t477 = t651 + t648 / 0.2e1 + t560;
t543 = legFrame(2,2);
t512 = sin(t543);
t515 = cos(t543);
t572 = t554 * t592;
t564 = pkin(2) * t572 + (t572 * t658 - t583) * t553;
t586 = t553 * (-t547 * pkin(3) + t504);
t603 = t512 * t554;
t613 = (t553 * pkin(3) + t561) * t547;
t433 = (-t477 * t603 + t515 * t613) * t656 + (t515 * t462 + t564 * t512) * t536 + t456 * t512 + t515 * t586;
t474 = 0.1e1 / (t536 * t553 - t592);
t636 = t433 * t474;
t597 = t515 * t554;
t434 = (t477 * t597 + t512 * t613) * t656 + (t512 * t462 - t564 * t515) * t536 - t456 * t515 + t512 * t586;
t635 = t434 * t474;
t556 = cos(qJ(1,1));
t523 = pkin(6) + t539;
t550 = sin(qJ(1,1));
t582 = pkin(1) * t556 + t550 * t523;
t591 = t535 * t549;
t457 = t582 * t591 + (t534 - 0.1e1) * t556 * pkin(3);
t463 = pkin(1) * t549 + (-pkin(3) + t647 + 0.2e1 * t650) * t535;
t478 = t650 + t647 / 0.2e1 + t560;
t544 = legFrame(1,2);
t513 = sin(t544);
t516 = cos(t544);
t571 = t556 * t591;
t563 = pkin(2) * t571 + (t571 * t658 - t582) * t555;
t585 = t555 * (-t549 * pkin(3) + t504);
t601 = t513 * t556;
t612 = (t555 * pkin(3) + t561) * t549;
t435 = (-t478 * t601 + t516 * t612) * t656 + (t516 * t463 + t563 * t513) * t536 + t457 * t513 + t516 * t585;
t475 = 0.1e1 / (t536 * t555 - t591);
t634 = t435 * t475;
t595 = t516 * t556;
t436 = (t478 * t595 + t513 * t612) * t656 + (t513 * t463 - t563 * t516) * t536 - t457 * t516 + t513 * t585;
t633 = t436 * t475;
t570 = -mrSges(3,2) * t545 + t642;
t452 = (-mrSges(3,1) * t551 - t570) * t536 + (mrSges(3,2) * t551 + t641) * t535 - t653;
t494 = t536 * pkin(2) + pkin(1);
t525 = pkin(7) + qJ(3,3);
t501 = cos(t525);
t458 = t521 * t552 + (-pkin(3) * t501 - t494) * t546;
t508 = 0.1e1 / t521;
t446 = (-t452 * t546 + t458 * t557) * t508;
t632 = t446 * t473;
t569 = -mrSges(3,2) * t547 + t642;
t453 = (-mrSges(3,1) * t553 - t569) * t536 + (mrSges(3,2) * t553 + t640) * t535 - t653;
t526 = pkin(7) + qJ(3,2);
t502 = cos(t526);
t459 = t522 * t554 + (-pkin(3) * t502 - t494) * t548;
t509 = 0.1e1 / t522;
t447 = (-t453 * t548 + t459 * t557) * t509;
t631 = t447 * t474;
t568 = -mrSges(3,2) * t549 + t642;
t454 = (-mrSges(3,1) * t555 - t568) * t536 + (mrSges(3,2) * t555 + t639) * t535 - t653;
t527 = pkin(7) + qJ(3,1);
t503 = cos(t527);
t460 = t523 * t556 + (-pkin(3) * t503 - t494) * t550;
t510 = 0.1e1 / t523;
t448 = (-t454 * t550 + t460 * t557) * t510;
t630 = t448 * t475;
t498 = sin(t525);
t464 = t514 * t498 - t501 * t605;
t495 = 0.1e1 / t501;
t629 = t464 * t495;
t628 = t464 * t508;
t465 = t511 * t498 + t501 * t599;
t627 = t465 * t495;
t626 = t465 * t508;
t499 = sin(t526);
t466 = t515 * t499 - t502 * t603;
t496 = 0.1e1 / t502;
t625 = t466 * t496;
t624 = t466 * t509;
t467 = t512 * t499 + t502 * t597;
t623 = t467 * t496;
t622 = t467 * t509;
t500 = sin(t527);
t468 = t516 * t500 - t503 * t601;
t497 = 0.1e1 / t503;
t621 = t468 * t497;
t620 = t468 * t510;
t469 = t513 * t500 + t503 * t595;
t619 = t469 * t497;
t618 = t469 * t510;
t617 = t473 * t557;
t616 = t474 * t557;
t615 = t475 * t557;
t611 = t495 * t511;
t610 = t495 * t514;
t609 = t496 * t512;
t608 = t496 * t515;
t607 = t497 * t513;
t606 = t497 * t516;
t562 = 1 / pkin(3);
t604 = t511 * t562;
t602 = t512 * t562;
t600 = t513 * t562;
t598 = t514 * t562;
t596 = t515 * t562;
t594 = t516 * t562;
t590 = t540 * t532;
t589 = t540 * t533;
t588 = t540 * t534;
t528 = -0.2e1 * pkin(2) * mrSges(3,2);
t581 = -0.2e1 * t504;
t580 = mrSges(3,2) * t504;
t482 = t537 * mrSges(3,2) - Ifges(3,6);
t485 = t537 * mrSges(3,1) - Ifges(3,5);
t449 = (-t482 * t551 - t485 * t545) * t536 + t535 * (t545 * t482 - t485 * t551);
t579 = t449 * t546 * t562;
t483 = t538 * mrSges(3,2) - Ifges(3,6);
t486 = t538 * mrSges(3,1) - Ifges(3,5);
t450 = (-t483 * t553 - t486 * t547) * t536 + t535 * (t547 * t483 - t486 * t553);
t578 = t450 * t548 * t562;
t484 = t539 * mrSges(3,2) - Ifges(3,6);
t487 = t539 * mrSges(3,1) - Ifges(3,5);
t451 = (-t484 * t555 - t487 * t549) * t536 + t535 * (t549 * t484 - t487 * t555);
t577 = t451 * t550 * t562;
t576 = t452 * t473 * t508;
t575 = t453 * t474 * t509;
t574 = t454 * t475 * t510;
t567 = pkin(2) ^ 2 * m(3) - Ifges(2,1) + Ifges(2,2) + t540;
t566 = (t557 * pkin(1) ^ 2) + (2 * mrSges(3,3) * pkin(5)) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t530 = mrSges(3,1) * t659;
t529 = 2 * t654;
t445 = (Ifges(3,3) * t600 + t451 * t618) * t497;
t444 = (Ifges(3,3) * t594 + t451 * t620) * t497;
t443 = (Ifges(3,3) * t602 + t450 * t622) * t496;
t442 = (Ifges(3,3) * t596 + t450 * t624) * t496;
t441 = (Ifges(3,3) * t604 + t449 * t626) * t495;
t440 = (Ifges(3,3) * t598 + t449 * t628) * t495;
t439 = (-0.2e1 * t588 + (t529 + 0.4e1 * t643) * t555 + t549 * t528 + t567) * t524 + (t530 * t555 + t568 * t659 + (t534 * t657 + t528 * t555 + 0.2e1 * (t540 * t555 - t654) * t549 + t646) * t535) * t536 + t588 + 0.2e1 * (-t580 - t643) * t555 + t639 * t581 + t539 ^ 2 * m(3) + t566 + ((m(2) * qJ(2,1) + t655) * qJ(2,1));
t438 = (-0.2e1 * t589 + (t529 + 0.4e1 * t644) * t553 + t547 * t528 + t567) * t524 + (t530 * t553 + t569 * t659 + (t533 * t657 + t528 * t553 + 0.2e1 * (t540 * t553 - t654) * t547 + t646) * t535) * t536 + t589 + 0.2e1 * (-t580 - t644) * t553 + t640 * t581 + t538 ^ 2 * m(3) + t566 + ((m(2) * qJ(2,2) + t655) * qJ(2,2));
t437 = (-0.2e1 * t590 + (t529 + 0.4e1 * t645) * t551 + t545 * t528 + t567) * t524 + (t530 * t551 + t570 * t659 + (t532 * t657 + t528 * t551 + 0.2e1 * (t540 * t551 - t654) * t545 + t646) * t535) * t536 + t590 + 0.2e1 * (-t580 - t645) * t551 + t641 * t581 + t537 ^ 2 * m(3) + t566 + ((m(2) * qJ(2,3) + t655) * qJ(2,3));
t430 = (-t439 * t550 + t454 * t460) * t510;
t429 = (-t438 * t548 + t453 * t459) * t509;
t428 = (-t437 * t546 + t452 * t458) * t508;
t427 = (t436 * t615 + t454 * t619) * t510;
t426 = (t435 * t615 + t454 * t621) * t510;
t425 = (t434 * t616 + t453 * t623) * t509;
t424 = (t433 * t616 + t453 * t625) * t509;
t423 = (t432 * t617 + t452 * t627) * t508;
t422 = (t431 * t617 + t452 * t629) * t508;
t421 = t436 * t574 + (t439 * t618 + t451 * t600) * t497;
t420 = t435 * t574 + (t439 * t620 + t451 * t594) * t497;
t419 = t434 * t575 + (t438 * t622 + t450 * t602) * t496;
t418 = t433 * t575 + (t438 * t624 + t450 * t596) * t496;
t417 = t432 * t576 + (t437 * t626 + t449 * t604) * t495;
t416 = t431 * t576 + (t437 * t628 + t449 * t598) * t495;
t1 = [m(4) + (t421 * t619 + t427 * t633) * t510 + (t419 * t623 + t425 * t635) * t509 + (t417 * t627 + t423 * t637) * t508 + (t441 * t611 + t443 * t609 + t445 * t607) * t562, (t421 * t621 + t427 * t634) * t510 + (t419 * t625 + t425 * t636) * t509 + (t417 * t629 + t423 * t638) * t508 + (t441 * t610 + t443 * t608 + t445 * t606) * t562, (-t421 * t550 + t427 * t460) * t510 + (-t419 * t548 + t425 * t459) * t509 + (-t417 * t546 + t423 * t458) * t508; (t420 * t619 + t426 * t633) * t510 + (t418 * t623 + t424 * t635) * t509 + (t416 * t627 + t422 * t637) * t508 + (t440 * t611 + t442 * t609 + t444 * t607) * t562, m(4) + (t420 * t621 + t426 * t634) * t510 + (t418 * t625 + t424 * t636) * t509 + (t416 * t629 + t422 * t638) * t508 + (t440 * t610 + t442 * t608 + t444 * t606) * t562, (-t420 * t550 + t426 * t460) * t510 + (-t418 * t548 + t424 * t459) * t509 + (-t416 * t546 + t422 * t458) * t508; (t436 * t630 + (t430 * t469 - t513 * t577) * t497) * t510 + (t434 * t631 + (t429 * t467 - t512 * t578) * t496) * t509 + (t432 * t632 + (t428 * t465 - t511 * t579) * t495) * t508, (t435 * t630 + (t430 * t468 - t516 * t577) * t497) * t510 + (t433 * t631 + (t429 * t466 - t515 * t578) * t496) * t509 + (t431 * t632 + (t428 * t464 - t514 * t579) * t495) * t508, m(4) + (-t430 * t550 + t448 * t460) * t510 + (-t429 * t548 + t447 * t459) * t509 + (-t428 * t546 + t446 * t458) * t508;];
MX  = t1;
