% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G2A0
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
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2022-11-07 13:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:06:00
% EndTime: 2022-11-07 13:06:02
% DurationCPUTime: 1.30s
% Computational Cost: add. (4653->278), mult. (7920->452), div. (348->9), fcn. (4632->23), ass. (0->197)
t531 = cos(pkin(7));
t525 = t531 ^ 2;
t635 = 0.2e1 * t525;
t634 = 2 * mrSges(3,1);
t633 = 2 * mrSges(3,3);
t545 = cos(qJ(2,3));
t527 = t545 ^ 2;
t632 = 0.2e1 * t527;
t547 = cos(qJ(2,2));
t528 = t547 ^ 2;
t631 = 0.2e1 * t528;
t549 = cos(qJ(2,1));
t529 = t549 ^ 2;
t630 = 0.2e1 * t529;
t629 = 0.2e1 * t531;
t628 = m(3) * pkin(1);
t627 = m(3) * pkin(2);
t539 = sin(qJ(2,3));
t626 = pkin(1) * t539;
t541 = sin(qJ(2,2));
t625 = pkin(1) * t541;
t543 = sin(qJ(2,1));
t624 = pkin(1) * t543;
t540 = sin(qJ(1,3));
t532 = qJ(3,3) + pkin(5);
t522 = pkin(6) + t532;
t546 = cos(qJ(1,3));
t562 = pkin(1) * t540 - t546 * t522;
t577 = pkin(3) * (t525 - 0.1e1);
t530 = sin(pkin(7));
t593 = t530 * t539;
t623 = (t540 * t577 + t562 * t593) * pkin(3);
t542 = sin(qJ(1,2));
t533 = qJ(3,2) + pkin(5);
t523 = pkin(6) + t533;
t548 = cos(qJ(1,2));
t561 = pkin(1) * t542 - t548 * t523;
t592 = t530 * t541;
t622 = (t542 * t577 + t561 * t592) * pkin(3);
t544 = sin(qJ(1,1));
t534 = qJ(3,1) + pkin(5);
t524 = pkin(6) + t534;
t550 = cos(qJ(1,1));
t560 = pkin(1) * t544 - t550 * t524;
t591 = t530 * t543;
t621 = (t544 * t577 + t560 * t591) * pkin(3);
t620 = t530 * pkin(3);
t619 = t531 * pkin(3);
t535 = Ifges(3,2) - Ifges(3,1);
t618 = mrSges(3,2) * t530;
t617 = Ifges(3,4) * t530;
t616 = t530 * mrSges(3,1);
t615 = Ifges(2,5) + (-mrSges(2,1) - t627) * pkin(5);
t499 = pkin(2) + t619;
t580 = pkin(3) * t593;
t474 = 0.1e1 / (t499 * t545 - t580);
t505 = 0.1e1 / t522;
t614 = t474 * t505;
t579 = pkin(3) * t592;
t475 = 0.1e1 / (t499 * t547 - t579);
t506 = 0.1e1 / t523;
t613 = t475 * t506;
t578 = pkin(3) * t591;
t476 = 0.1e1 / (t499 * t549 - t578);
t507 = 0.1e1 / t524;
t612 = t476 * t507;
t566 = pkin(3) * cos(qJ(2,3) + pkin(7)) + t545 * pkin(2);
t480 = 0.1e1 / t566;
t536 = legFrame(3,2);
t508 = sin(t536);
t611 = t480 * t508;
t511 = cos(t536);
t610 = t480 * t511;
t565 = pkin(3) * cos(qJ(2,2) + pkin(7)) + t547 * pkin(2);
t481 = 0.1e1 / t565;
t537 = legFrame(2,2);
t509 = sin(t537);
t609 = t481 * t509;
t512 = cos(t537);
t608 = t481 * t512;
t564 = pkin(3) * cos(qJ(2,1) + pkin(7)) + t549 * pkin(2);
t482 = 0.1e1 / t564;
t538 = legFrame(1,2);
t510 = sin(t538);
t607 = t482 * t510;
t513 = cos(t538);
t606 = t482 * t513;
t605 = t499 * t511;
t604 = t499 * t512;
t603 = t499 * t513;
t602 = t508 * t499;
t601 = t508 * t540;
t600 = t509 * t499;
t599 = t509 * t542;
t598 = t510 * t499;
t597 = t510 * t544;
t596 = t511 * t540;
t595 = t512 * t542;
t594 = t513 * t544;
t496 = t535 * t525;
t554 = pkin(2) ^ 2;
t590 = m(3) * t554 - 0.2e1 * pkin(2) * t618;
t589 = pkin(2) * t634;
t588 = -0.2e1 * (mrSges(2,2) + t616) * pkin(1);
t587 = pkin(2) * t619;
t586 = t511 * t620;
t585 = t512 * t620;
t584 = t513 * t620;
t583 = t508 * t620;
t582 = t509 * t620;
t581 = t510 * t620;
t489 = t532 * mrSges(3,2) - Ifges(3,6);
t492 = -t532 * mrSges(3,1) + Ifges(3,5);
t567 = -pkin(5) * mrSges(2,2) + Ifges(2,6);
t449 = (t489 * t530 + t492 * t531 + (-m(3) * qJ(3,3) - mrSges(3,3)) * pkin(2) + t615) * t539 + (-t489 * t531 + t492 * t530 + t567) * t545;
t576 = t449 * t614;
t490 = t533 * mrSges(3,2) - Ifges(3,6);
t493 = -t533 * mrSges(3,1) + Ifges(3,5);
t450 = (t490 * t530 + t493 * t531 + (-m(3) * qJ(3,2) - mrSges(3,3)) * pkin(2) + t615) * t541 + (-t490 * t531 + t493 * t530 + t567) * t547;
t575 = t450 * t613;
t491 = t534 * mrSges(3,2) - Ifges(3,6);
t494 = -t534 * mrSges(3,1) + Ifges(3,5);
t451 = (t491 * t530 + t494 * t531 + (-m(3) * qJ(3,1) - mrSges(3,3)) * pkin(2) + t615) * t543 + (-t491 * t531 + t494 * t530 + t567) * t549;
t574 = t451 * t612;
t573 = t449 * t611;
t572 = t450 * t609;
t571 = t451 * t607;
t570 = t449 * t610;
t569 = t450 * t608;
t568 = t451 * t606;
t559 = -t618 + t627;
t563 = (t634 * t531 + (2 * mrSges(2,1)) + 0.2e1 * t559) * pkin(1);
t558 = mrSges(3,2) * t531 + t616;
t557 = 0.2e1 * Ifges(3,4) * t635 + 0.2e1 * (-pkin(2) * mrSges(3,2) - t530 * t535) * t531 - 0.2e1 * pkin(2) * t616 + (2 * Ifges(2,4)) - 0.2e1 * Ifges(3,4);
t556 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + 0.2e1 * (mrSges(2,3) + mrSges(3,3)) * pkin(5) + (m(2) + m(3)) * pkin(1) ^ 2 + m(2) * pkin(5) ^ 2 - t496;
t553 = pkin(3) ^ 2;
t555 = t553 * t635 - t553 + t554 + 0.2e1 * t587;
t501 = pkin(1) * t620;
t487 = -t620 + t624;
t486 = -t620 + t625;
t485 = -t620 + t626;
t479 = -mrSges(3,1) * t531 - t559;
t478 = t587 + t554 / 0.2e1 + (t525 - 0.1e1 / 0.2e1) * t553;
t477 = t531 * t589 + Ifges(2,3) + Ifges(3,3) + t590;
t473 = -0.2e1 * t544 * t578 + t560;
t472 = -0.2e1 * t542 * t579 + t561;
t471 = -0.2e1 * t540 * t580 + t562;
t470 = t544 * t524 + (pkin(1) + t564) * t550;
t469 = t542 * t523 + (pkin(1) + t565) * t548;
t468 = t540 * t522 + (pkin(1) + t566) * t546;
t467 = t555 * t543 + t501;
t466 = t555 * t541 + t501;
t465 = t555 * t539 + t501;
t461 = 0.2e1 * t496 + (t589 + 0.4e1 * t617) * t531 + Ifges(2,2) - Ifges(2,1) + t590 - t535;
t460 = t479 * t549 + t558 * t543 - t628;
t459 = t479 * t547 + t558 * t541 - t628;
t458 = t479 * t545 + t558 * t539 - t628;
t457 = (t499 * t594 + t581) * t549 + t543 * (-t544 * t584 + t598);
t456 = (t499 * t595 + t582) * t547 + t541 * (-t542 * t585 + t600);
t455 = (t499 * t596 + t583) * t545 + t539 * (-t540 * t586 + t602);
t454 = (-t499 * t597 + t584) * t549 + (t544 * t581 + t603) * t543;
t453 = (-t499 * t599 + t585) * t547 + (t542 * t582 + t604) * t541;
t452 = (-t499 * t601 + t586) * t545 + (t540 * t583 + t605) * t539;
t448 = (m(3) * t470 + t460 * t550) * t507;
t447 = (m(3) * t469 + t459 * t548) * t506;
t446 = (m(3) * t468 + t458 * t546) * t505;
t445 = t461 * t529 + (t557 * t543 + t563) * t549 + (-mrSges(3,2) * t624 - t617) * t629 + t543 * t588 + t534 ^ 2 * m(3) + qJ(3,1) * t633 + t556;
t444 = t461 * t528 + (t557 * t541 + t563) * t547 + (-mrSges(3,2) * t625 - t617) * t629 + t541 * t588 + t533 ^ 2 * m(3) + qJ(3,2) * t633 + t556;
t443 = t461 * t527 + (t557 * t539 + t563) * t545 + (-mrSges(3,2) * t626 - t617) * t629 + t539 * t588 + t532 ^ 2 * m(3) + qJ(3,3) * t633 + t556;
t442 = (t478 * t594 + t499 * t581) * t630 + (t467 * t510 + t473 * t603) * t549 - t513 * t621 + t487 * t598;
t441 = (-t478 * t597 + t499 * t584) * t630 + (t513 * t467 - t473 * t598) * t549 + t510 * t621 + t487 * t603;
t440 = (t478 * t595 + t499 * t582) * t631 + (t466 * t509 + t472 * t604) * t547 - t512 * t622 + t486 * t600;
t439 = (-t478 * t599 + t499 * t585) * t631 + (t512 * t466 - t472 * t600) * t547 + t509 * t622 + t486 * t604;
t438 = (t478 * t596 + t499 * t583) * t632 + (t465 * t508 + t471 * t605) * t545 - t511 * t623 + t485 * t602;
t437 = (-t478 * t601 + t499 * t586) * t632 + (t511 * t465 - t471 * t602) * t545 + t508 * t623 + t485 * t605;
t436 = t457 * t574 + t477 * t607;
t435 = t456 * t575 + t477 * t609;
t434 = t455 * t576 + t477 * t611;
t433 = t454 * t574 + t477 * t606;
t432 = t453 * t575 + t477 * t608;
t431 = t452 * t576 + t477 * t610;
t430 = (t445 * t550 + t460 * t470) * t507;
t429 = (t444 * t548 + t459 * t469) * t506;
t428 = (t443 * t546 + t458 * t468) * t505;
t427 = (m(3) * t442 + t457 * t460) * t612;
t426 = (m(3) * t440 + t456 * t459) * t613;
t425 = (m(3) * t438 + t455 * t458) * t614;
t424 = (m(3) * t441 + t454 * t460) * t612;
t423 = (m(3) * t439 + t453 * t459) * t613;
t422 = (m(3) * t437 + t452 * t458) * t614;
t421 = t571 + (t442 * t460 + t445 * t457) * t612;
t420 = t572 + (t440 * t459 + t444 * t456) * t613;
t419 = t573 + (t438 * t458 + t443 * t455) * t614;
t418 = t568 + (t441 * t460 + t445 * t454) * t612;
t417 = t569 + (t439 * t459 + t444 * t453) * t613;
t416 = t570 + (t437 * t458 + t443 * t452) * t614;
t1 = [t434 * t611 + t435 * t609 + t436 * t607 + m(4) + (t421 * t457 + t427 * t442) * t612 + (t420 * t456 + t426 * t440) * t613 + (t419 * t455 + t425 * t438) * t614, t434 * t610 + t435 * t608 + t436 * t606 + (t421 * t454 + t427 * t441) * t612 + (t420 * t453 + t426 * t439) * t613 + (t419 * t452 + t425 * t437) * t614, (t421 * t550 + t427 * t470) * t507 + (t420 * t548 + t426 * t469) * t506 + (t419 * t546 + t425 * t468) * t505; t431 * t611 + t432 * t609 + t433 * t607 + (t418 * t457 + t424 * t442) * t612 + (t417 * t456 + t423 * t440) * t613 + (t416 * t455 + t422 * t438) * t614, t431 * t610 + t432 * t608 + t433 * t606 + m(4) + (t418 * t454 + t424 * t441) * t612 + (t417 * t453 + t423 * t439) * t613 + (t416 * t452 + t422 * t437) * t614, (t418 * t550 + t424 * t470) * t507 + (t417 * t548 + t423 * t469) * t506 + (t416 * t546 + t422 * t468) * t505; (t550 * t571 + (t430 * t457 + t442 * t448) * t476) * t507 + (t548 * t572 + (t429 * t456 + t440 * t447) * t475) * t506 + (t546 * t573 + (t428 * t455 + t438 * t446) * t474) * t505, (t550 * t568 + (t430 * t454 + t441 * t448) * t476) * t507 + (t548 * t569 + (t429 * t453 + t439 * t447) * t475) * t506 + (t546 * t570 + (t428 * t452 + t437 * t446) * t474) * t505, m(4) + (t430 * t550 + t448 * t470) * t507 + (t429 * t548 + t447 * t469) * t506 + (t428 * t546 + t446 * t468) * t505;];
MX  = t1;
