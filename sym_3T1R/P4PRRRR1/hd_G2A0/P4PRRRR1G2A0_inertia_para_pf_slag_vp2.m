% Calculate inertia matrix for parallel robot
% P4PRRRR1G2A0
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
% koppelP [4x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-07 11:00
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-07 10:58:43
% EndTime: 2020-08-07 10:58:44
% DurationCPUTime: 1.54s
% Computational Cost: add. (2348->282), mult. (4404->493), div. (1884->13), fcn. (5440->26), ass. (0->234)
t645 = 2 * Ifges(3,4);
t644 = Ifges(3,1) + Ifges(2,3);
t524 = sin(qJ(3,4));
t526 = cos(qJ(3,4));
t493 = mrSges(3,1) * t524 + mrSges(3,2) * t526;
t512 = 0.1e1 / t526;
t643 = t493 * t512;
t525 = sin(qJ(2,4));
t642 = t493 * t525;
t534 = sin(qJ(3,3));
t540 = cos(qJ(3,3));
t497 = mrSges(3,1) * t534 + mrSges(3,2) * t540;
t518 = 0.1e1 / t540;
t641 = t497 * t518;
t535 = sin(qJ(2,3));
t640 = t497 * t535;
t536 = sin(qJ(3,2));
t542 = cos(qJ(3,2));
t498 = mrSges(3,1) * t536 + mrSges(3,2) * t542;
t520 = 0.1e1 / t542;
t639 = t498 * t520;
t537 = sin(qJ(2,2));
t638 = t498 * t537;
t538 = sin(qJ(3,1));
t544 = cos(qJ(3,1));
t499 = mrSges(3,1) * t538 + mrSges(3,2) * t544;
t522 = 0.1e1 / t544;
t637 = t499 * t522;
t539 = sin(qJ(2,1));
t636 = t499 * t539;
t530 = legFrame(4,2);
t501 = sin(t530);
t557 = 0.1e1 / pkin(2);
t635 = t501 * t557;
t531 = legFrame(3,2);
t502 = sin(t531);
t634 = t502 * t557;
t532 = legFrame(2,2);
t503 = sin(t532);
t633 = t503 * t557;
t533 = legFrame(1,2);
t504 = sin(t533);
t632 = t504 * t557;
t505 = cos(t530);
t631 = t505 * t557;
t506 = cos(t531);
t630 = t506 * t557;
t507 = cos(t532);
t629 = t507 * t557;
t508 = cos(t533);
t628 = t508 * t557;
t511 = 0.1e1 / t525;
t627 = t511 * t512;
t527 = cos(qJ(2,4));
t626 = t511 * t527;
t625 = t512 * t557;
t515 = 0.1e1 / t535;
t624 = t515 * t518;
t541 = cos(qJ(2,3));
t623 = t515 * t541;
t516 = 0.1e1 / t537;
t622 = t516 * t520;
t543 = cos(qJ(2,2));
t621 = t516 * t543;
t517 = 0.1e1 / t539;
t620 = t517 * t522;
t545 = cos(qJ(2,1));
t619 = t517 * t545;
t618 = t518 * t557;
t617 = t520 * t557;
t616 = t522 * t557;
t615 = t525 * t526;
t614 = t535 * t540;
t613 = t537 * t542;
t612 = t539 * t544;
t472 = t501 * t615 - t505 * t524;
t611 = t472 * t627;
t473 = t501 * t524 + t505 * t615;
t610 = t473 * t627;
t474 = t502 * t614 - t506 * t534;
t609 = t474 * t624;
t475 = t503 * t613 - t507 * t536;
t608 = t475 * t622;
t476 = t504 * t612 - t508 * t538;
t607 = t476 * t620;
t477 = t502 * t534 + t506 * t614;
t606 = t477 * t624;
t478 = t503 * t536 + t507 * t613;
t605 = t478 * t622;
t479 = t504 * t538 + t508 * t612;
t604 = t479 * t620;
t514 = m(1) + m(2) + m(3);
t603 = t514 * t627;
t602 = t514 * t624;
t601 = t514 * t622;
t600 = t514 * t620;
t546 = xP(4);
t509 = sin(t546);
t510 = cos(t546);
t547 = mrSges(4,2);
t548 = mrSges(4,1);
t599 = -t509 * t547 + t510 * t548;
t558 = t526 ^ 2;
t598 = 0.1e1 / t558 * t524 * t626;
t559 = t540 ^ 2;
t597 = 0.1e1 / t559 * t534 * t623;
t560 = t542 ^ 2;
t596 = 0.1e1 / t560 * t536 * t621;
t561 = t544 ^ 2;
t595 = 0.1e1 / t561 * t538 * t619;
t549 = koppelP(4,2);
t553 = koppelP(4,1);
t484 = -t509 * t553 - t510 * t549;
t488 = -t509 * t549 + t510 * t553;
t594 = t484 * t505 - t488 * t501;
t550 = koppelP(3,2);
t554 = koppelP(3,1);
t485 = -t509 * t554 - t510 * t550;
t489 = -t509 * t550 + t510 * t554;
t593 = t485 * t506 - t489 * t502;
t551 = koppelP(2,2);
t555 = koppelP(2,1);
t486 = -t509 * t555 - t510 * t551;
t490 = -t509 * t551 + t510 * t555;
t592 = t486 * t507 - t490 * t503;
t552 = koppelP(1,2);
t556 = koppelP(1,1);
t487 = -t509 * t556 - t510 * t552;
t491 = -t509 * t552 + t510 * t556;
t591 = t487 * t508 - t491 * t504;
t590 = -t509 * t548 - t510 * t547;
t492 = Ifges(3,5) * t524 + Ifges(3,6) * t526;
t589 = -Ifges(3,3) * t512 + t492 * t598;
t494 = Ifges(3,5) * t534 + Ifges(3,6) * t540;
t588 = -Ifges(3,3) * t518 + t494 * t597;
t495 = Ifges(3,5) * t536 + Ifges(3,6) * t542;
t587 = -Ifges(3,3) * t520 + t495 * t596;
t496 = Ifges(3,5) * t538 + Ifges(3,6) * t544;
t586 = -Ifges(3,3) * t522 + t496 * t595;
t436 = (t472 * t484 + t473 * t488) * t627;
t448 = t594 * t557 * t598;
t456 = t594 * t625;
t529 = mrSges(2,2) - mrSges(3,3);
t468 = (t526 * mrSges(3,1) - t524 * mrSges(3,2) + mrSges(2,1)) * t527 - t525 * t529;
t528 = -Ifges(3,1) + Ifges(3,2);
t480 = t524 * t526 * t645 + t528 * t558 + t644;
t408 = t436 * t468 + t448 * t480 - t456 * t492;
t416 = -t456 * Ifges(3,3) - t436 * t642 + t448 * t492;
t585 = t408 * t598 - t416 * t512;
t437 = (t474 * t485 + t477 * t489) * t624;
t449 = t593 * t557 * t597;
t457 = t593 * t618;
t469 = (t540 * mrSges(3,1) - t534 * mrSges(3,2) + mrSges(2,1)) * t541 - t535 * t529;
t481 = t534 * t540 * t645 + t528 * t559 + t644;
t409 = t437 * t469 + t449 * t481 - t457 * t494;
t417 = -t457 * Ifges(3,3) - t437 * t640 + t449 * t494;
t584 = t409 * t597 - t417 * t518;
t438 = (t475 * t486 + t478 * t490) * t622;
t450 = t592 * t557 * t596;
t458 = t592 * t617;
t470 = (t542 * mrSges(3,1) - t536 * mrSges(3,2) + mrSges(2,1)) * t543 - t537 * t529;
t482 = t536 * t542 * t645 + t528 * t560 + t644;
t410 = t438 * t470 + t450 * t482 - t458 * t495;
t418 = -t458 * Ifges(3,3) - t438 * t638 + t450 * t495;
t583 = t410 * t596 - t418 * t520;
t439 = (t476 * t487 + t479 * t491) * t620;
t451 = t591 * t557 * t595;
t459 = t591 * t616;
t471 = (t544 * mrSges(3,1) - t538 * mrSges(3,2) + mrSges(2,1)) * t545 - t539 * t529;
t483 = t538 * t544 * t645 + t528 * t561 + t644;
t411 = t439 * t471 + t451 * t483 - t459 * t496;
t419 = -t459 * Ifges(3,3) - t439 * t636 + t451 * t496;
t582 = t411 * t595 - t419 * t522;
t569 = t480 * t598 - t492 * t512;
t420 = t468 * t611 + t569 * t631;
t440 = -t472 * t643 + t589 * t631;
t581 = t420 * t598 - t440 * t512;
t421 = t468 * t610 - t569 * t635;
t441 = -t473 * t643 - t589 * t635;
t580 = t421 * t598 - t441 * t512;
t568 = t481 * t597 - t494 * t518;
t422 = t469 * t609 + t568 * t630;
t442 = -t474 * t641 + t588 * t630;
t579 = t422 * t597 - t442 * t518;
t567 = t482 * t596 - t495 * t520;
t423 = t470 * t608 + t567 * t629;
t443 = -t475 * t639 + t587 * t629;
t578 = t423 * t596 - t443 * t520;
t566 = t483 * t595 - t496 * t522;
t424 = t471 * t607 + t566 * t628;
t444 = -t476 * t637 + t586 * t628;
t577 = t424 * t595 - t444 * t522;
t425 = t469 * t606 - t568 * t634;
t445 = -t477 * t641 - t588 * t634;
t576 = t425 * t597 - t445 * t518;
t426 = t470 * t605 - t567 * t633;
t446 = -t478 * t639 - t587 * t633;
t575 = t426 * t596 - t446 * t520;
t427 = t471 * t604 - t566 * t632;
t447 = -t479 * t637 - t586 * t632;
t574 = t427 * t595 - t447 * t522;
t452 = (t468 * t527 - t480 * t625) * t511;
t464 = -t511 * t492 * t625 - t527 * t493;
t573 = t452 * t598 - t464 * t512;
t453 = (t469 * t541 - t481 * t618) * t515;
t465 = -t515 * t494 * t618 - t541 * t497;
t572 = t453 * t597 - t465 * t518;
t454 = (t470 * t543 - t482 * t617) * t516;
t466 = -t516 * t495 * t617 - t543 * t498;
t571 = t454 * t596 - t466 * t520;
t455 = (t471 * t545 - t483 * t616) * t517;
t467 = -t517 * t496 * t616 - t545 * t499;
t570 = t455 * t595 - t467 * t522;
t565 = t468 * t598 + t512 * t642;
t564 = t469 * t597 + t518 * t640;
t563 = t470 * t596 + t520 * t638;
t562 = t471 * t595 + t522 * t636;
t463 = (-t471 * t616 + t514 * t545) * t517;
t462 = (-t470 * t617 + t514 * t543) * t516;
t461 = (-t469 * t618 + t514 * t541) * t515;
t460 = (-t468 * t625 + t514 * t527) * t511;
t435 = t479 * t600 - t562 * t632;
t434 = t478 * t601 - t563 * t633;
t433 = t477 * t602 - t564 * t634;
t432 = t476 * t600 + t562 * t628;
t431 = t475 * t601 + t563 * t629;
t430 = t474 * t602 + t564 * t630;
t429 = t473 * t603 - t565 * t635;
t428 = t472 * t603 + t565 * t631;
t415 = t439 * t514 + t451 * t471 + t459 * t636;
t414 = t438 * t514 + t450 * t470 + t458 * t638;
t413 = t437 * t514 + t449 * t469 + t457 * t640;
t412 = t436 * t514 + t448 * t468 + t456 * t642;
t1 = [t428 * t611 + t430 * t609 + t431 * t608 + t432 * t607 + m(4) + (t581 * t505 + t579 * t506 + t578 * t507 + t577 * t508) * t557, t428 * t610 + t430 * t606 + t431 * t605 + t432 * t604 + (-t581 * t501 - t579 * t502 - t578 * t503 - t577 * t504) * t557, t428 * t626 + t430 * t623 + t431 * t621 + t432 * t619 + (-t420 * t627 - t422 * t624 - t423 * t622 - t424 * t620) * t557, t420 * t448 + t422 * t449 + t423 * t450 + t424 * t451 + t428 * t436 + t430 * t437 + t431 * t438 + t432 * t439 - t440 * t456 - t442 * t457 - t443 * t458 - t444 * t459 + t590; t429 * t611 + t433 * t609 + t434 * t608 + t435 * t607 + (t580 * t505 + t576 * t506 + t575 * t507 + t574 * t508) * t557, t429 * t610 + t433 * t606 + t434 * t605 + t435 * t604 + m(4) + (-t580 * t501 - t576 * t502 - t575 * t503 - t574 * t504) * t557, t429 * t626 + t433 * t623 + t434 * t621 + t435 * t619 + (-t421 * t627 - t425 * t624 - t426 * t622 - t427 * t620) * t557, t421 * t448 + t425 * t449 + t426 * t450 + t427 * t451 + t429 * t436 + t433 * t437 + t434 * t438 + t435 * t439 - t441 * t456 - t445 * t457 - t446 * t458 - t447 * t459 + t599; t460 * t611 + t461 * t609 + t462 * t608 + t463 * t607 + (t573 * t505 + t572 * t506 + t571 * t507 + t570 * t508) * t557, t460 * t610 + t461 * t606 + t462 * t605 + t463 * t604 + (-t573 * t501 - t572 * t502 - t571 * t503 - t570 * t504) * t557, t460 * t626 + t461 * t623 + t462 * t621 + t463 * t619 + m(4) + (-t452 * t627 - t453 * t624 - t454 * t622 - t455 * t620) * t557, t460 * t436 + t461 * t437 + t462 * t438 + t463 * t439 + t452 * t448 + t453 * t449 + t454 * t450 + t455 * t451 - t464 * t456 - t465 * t457 - t466 * t458 - t467 * t459; t412 * t611 + t413 * t609 + t414 * t608 + t415 * t607 + (t585 * t505 + t584 * t506 + t583 * t507 + t582 * t508) * t557 + t590, t412 * t610 + t413 * t606 + t414 * t605 + t415 * t604 + (-t585 * t501 - t584 * t502 - t583 * t503 - t582 * t504) * t557 + t599, t412 * t626 + t413 * t623 + t414 * t621 + t415 * t619 + (-t408 * t627 - t409 * t624 - t410 * t622 - t411 * t620) * t557, t408 * t448 + t409 * t449 + t410 * t450 + t411 * t451 + t412 * t436 + t413 * t437 + t414 * t438 + t415 * t439 - t416 * t456 - t417 * t457 - t418 * t458 - t419 * t459 + Ifges(4,3);];
MX  = t1;
