% Calculate inertia matrix for parallel robot
% P3RRPRR12V2G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4]';
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
% Datum: 2020-08-06 19:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V2G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:15:03
% EndTime: 2020-08-06 19:15:05
% DurationCPUTime: 1.66s
% Computational Cost: add. (3993->266), mult. (4713->454), div. (396->6), fcn. (3423->18), ass. (0->183)
t583 = pkin(1) ^ 2 + pkin(5) ^ 2;
t554 = pkin(2) + pkin(3);
t540 = sin(qJ(2,3));
t619 = pkin(1) * t540;
t542 = sin(qJ(2,2));
t618 = pkin(1) * t542;
t544 = sin(qJ(2,1));
t617 = pkin(1) * t544;
t616 = mrSges(3,3) - mrSges(2,2);
t615 = (-Ifges(2,1) - Ifges(3,1));
t614 = Ifges(2,6) - Ifges(3,6);
t613 = (mrSges(3,3) * qJ(3,1));
t612 = (mrSges(3,3) * qJ(3,2));
t611 = (mrSges(3,3) * qJ(3,3));
t537 = legFrame(3,3);
t516 = sin(t537);
t519 = cos(t537);
t541 = sin(qJ(1,3));
t547 = cos(qJ(1,3));
t470 = t516 * t541 - t519 * t547;
t589 = t540 * qJ(3,3);
t498 = pkin(1) + t589;
t553 = pkin(5) - pkin(6);
t509 = t553 * t547;
t477 = t541 * t498 - t509;
t506 = t541 * t553;
t575 = t498 * t547 + t506;
t546 = cos(qJ(2,3));
t586 = t554 * t546;
t449 = t470 * t586 + t477 * t516 - t575 * t519;
t610 = t449 * t546;
t471 = t516 * t547 + t519 * t541;
t450 = t471 * t586 + t477 * t519 + t516 * t575;
t609 = t450 * t546;
t538 = legFrame(2,3);
t517 = sin(t538);
t520 = cos(t538);
t543 = sin(qJ(1,2));
t549 = cos(qJ(1,2));
t472 = t517 * t543 - t520 * t549;
t588 = t542 * qJ(3,2);
t500 = pkin(1) + t588;
t510 = t553 * t549;
t479 = t543 * t500 - t510;
t507 = t543 * t553;
t573 = t500 * t549 + t507;
t548 = cos(qJ(2,2));
t585 = t554 * t548;
t451 = t472 * t585 + t479 * t517 - t573 * t520;
t608 = t451 * t548;
t473 = t517 * t549 + t520 * t543;
t452 = t473 * t585 + t479 * t520 + t517 * t573;
t607 = t452 * t548;
t539 = legFrame(1,3);
t518 = sin(t539);
t521 = cos(t539);
t545 = sin(qJ(1,1));
t551 = cos(qJ(1,1));
t474 = t518 * t545 - t521 * t551;
t587 = t544 * qJ(3,1);
t502 = pkin(1) + t587;
t511 = t553 * t551;
t481 = t545 * t502 - t511;
t508 = t545 * t553;
t571 = t502 * t551 + t508;
t550 = cos(qJ(2,1));
t584 = t554 * t550;
t453 = t474 * t584 + t481 * t518 - t571 * t521;
t606 = t453 * t550;
t475 = t518 * t551 + t521 * t545;
t454 = t475 * t584 + t481 * t521 + t518 * t571;
t605 = t454 * t550;
t530 = 2 * t611;
t557 = qJ(3,3) ^ 2;
t564 = pkin(2) ^ 2;
t533 = 2 * mrSges(3,1) * pkin(2);
t581 = Ifges(3,2) + Ifges(2,3) + t533;
t482 = m(3) * (t557 + t564) + t530 + t581;
t489 = -t546 * qJ(3,3) + t554 * t540;
t558 = 0.1e1 / qJ(3,3);
t580 = m(3) * pkin(2) + mrSges(3,1);
t458 = (t482 * t540 - t489 * t580) * t558;
t604 = t458 * t546;
t531 = 2 * t612;
t559 = qJ(3,2) ^ 2;
t483 = m(3) * (t559 + t564) + t531 + t581;
t490 = -t548 * qJ(3,2) + t554 * t542;
t560 = 0.1e1 / qJ(3,2);
t459 = (t483 * t542 - t490 * t580) * t560;
t603 = t459 * t548;
t532 = 2 * t613;
t561 = qJ(3,1) ^ 2;
t484 = m(3) * (t561 + t564) + t532 + t581;
t491 = -t550 * qJ(3,1) + t554 * t544;
t562 = 0.1e1 / qJ(3,1);
t460 = (t484 * t544 - t491 * t580) * t562;
t602 = t460 * t550;
t515 = mrSges(2,1) + t580;
t488 = (pkin(2) * mrSges(3,2)) + t515 * pkin(5) - Ifges(3,4) - Ifges(2,5);
t512 = m(3) * qJ(3,3) + t616;
t461 = t546 * (qJ(3,3) * mrSges(3,2) + t512 * pkin(5) + t614) - t488 * t540;
t601 = t461 * t546;
t513 = m(3) * qJ(3,2) + t616;
t462 = t548 * (qJ(3,2) * mrSges(3,2) + t513 * pkin(5) + t614) - t488 * t542;
t600 = t462 * t548;
t514 = qJ(3,1) * m(3) + t616;
t463 = t550 * (qJ(3,1) * mrSges(3,2) + t514 * pkin(5) + t614) - t488 * t544;
t599 = t463 * t550;
t598 = t482 * t546;
t597 = t483 * t548;
t596 = t484 * t550;
t522 = m(3) * pkin(5) + mrSges(3,2);
t595 = t522 * t540;
t594 = t522 * t542;
t593 = t522 * t544;
t592 = t580 * t546;
t591 = t580 * t548;
t590 = t580 * t550;
t582 = pkin(2) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t534 = t546 ^ 2;
t579 = (qJ(3,3) + t554) * (-qJ(3,3) + t554) * t534;
t535 = t548 ^ 2;
t578 = (qJ(3,2) + t554) * (-qJ(3,2) + t554) * t535;
t536 = t550 ^ 2;
t577 = (qJ(3,1) + t554) * (-qJ(3,1) + t554) * t536;
t497 = pkin(1) + 0.2e1 * t589;
t576 = t497 * t547 + t506;
t499 = pkin(1) + 0.2e1 * t588;
t574 = t499 * t549 + t507;
t501 = pkin(1) + 0.2e1 * t587;
t572 = t501 * t551 + t508;
t570 = Ifges(2,2) + Ifges(3,3) + t533 + t615;
t569 = Ifges(1,3) + 0.2e1 * (mrSges(3,2) + mrSges(2,3)) * pkin(5) - t615 + t583 * m(2);
t503 = qJ(3,3) + t619;
t568 = t503 * t547 + t540 * t506;
t504 = qJ(3,2) + t618;
t567 = t504 * t549 + t542 * t507;
t505 = qJ(3,1) + t617;
t566 = t505 * t551 + t544 * t508;
t495 = 0.2e1 * t515 * pkin(1);
t487 = 0.1e1 / (t502 + t584);
t486 = 0.1e1 / (t500 + t585);
t485 = 0.1e1 / (t498 + t586);
t480 = t545 * t501 - t511;
t478 = t543 * t499 - t510;
t476 = t541 * t497 - t509;
t469 = t545 * t505 - t544 * t511;
t468 = t543 * t504 - t542 * t510;
t467 = t541 * t503 - t540 * t509;
t466 = (m(3) * t491 - t544 * t580) * t562;
t465 = (m(3) * t490 - t542 * t580) * t560;
t464 = (m(3) * t489 - t540 * t580) * t558;
t457 = (t491 * t522 + t463) * t562 * t544;
t456 = (t490 * t522 + t462) * t560 * t542;
t455 = (t489 * t522 + t461) * t558 * t540;
t448 = ((-t561 + t564) * m(3) - (2 * t613) + t570) * t536 + (0.2e1 * (t580 * qJ(3,1) + t582) * t544 + t495) * t550 + 0.2e1 * t514 * t617 + (t561 + t583) * m(3) + t532 + t569;
t447 = ((-t559 + t564) * m(3) - (2 * t612) + t570) * t535 + (0.2e1 * (t580 * qJ(3,2) + t582) * t542 + t495) * t548 + 0.2e1 * t513 * t618 + (t559 + t583) * m(3) + t531 + t569;
t446 = ((-t557 + t564) * m(3) - (2 * t611) + t570) * t534 + (0.2e1 * (t580 * qJ(3,3) + t582) * t540 + t495) * t546 + 0.2e1 * t512 * t619 + (t557 + t583) * m(3) + t530 + t569;
t445 = t475 * t577 + (t480 * t521 + t572 * t518) * t584 + (t469 * t521 + t566 * t518) * qJ(3,1);
t444 = -t474 * t577 - (t518 * t480 - t572 * t521) * t584 - (t518 * t469 - t566 * t521) * qJ(3,1);
t443 = t473 * t578 + (t478 * t520 + t574 * t517) * t585 + (t468 * t520 + t567 * t517) * qJ(3,2);
t442 = -t472 * t578 - (t517 * t478 - t574 * t520) * t585 - (t517 * t468 - t567 * t520) * qJ(3,2);
t441 = t471 * t579 + (t476 * t519 + t576 * t516) * t586 + (t467 * t519 + t568 * t516) * qJ(3,3);
t440 = -t470 * t579 - (t516 * t476 - t576 * t519) * t586 - (t516 * t467 - t568 * t519) * qJ(3,3);
t439 = (-t475 * t593 + (m(3) * t444 + t453 * t590) * t562) * t487;
t438 = (-t474 * t593 + (m(3) * t445 - t454 * t590) * t562) * t487;
t437 = (-t473 * t594 + (m(3) * t442 + t451 * t591) * t560) * t486;
t436 = (-t472 * t594 + (m(3) * t443 - t452 * t591) * t560) * t486;
t435 = (-t471 * t595 + (m(3) * t440 + t449 * t592) * t558) * t485;
t434 = (-t470 * t595 + (m(3) * t441 - t450 * t592) * t558) * t485;
t433 = (-t463 * t475 + (-t444 * t580 - t453 * t596) * t562) * t487;
t432 = (-t463 * t474 + (-t445 * t580 + t454 * t596) * t562) * t487;
t431 = (-t462 * t473 + (-t442 * t580 - t451 * t597) * t560) * t486;
t430 = (-t462 * t472 + (-t443 * t580 + t452 * t597) * t560) * t486;
t429 = (-t461 * t471 + (-t440 * t580 - t449 * t598) * t558) * t485;
t428 = (-t461 * t470 + (-t441 * t580 + t450 * t598) * t558) * t485;
t427 = (-t448 * t475 + (t444 * t593 - t453 * t599) * t562) * t487;
t426 = (-t448 * t474 + (t445 * t593 + t454 * t599) * t562) * t487;
t425 = (-t447 * t473 + (t442 * t594 - t451 * t600) * t560) * t486;
t424 = (-t447 * t472 + (t443 * t594 + t452 * t600) * t560) * t486;
t423 = (-t446 * t471 + (t440 * t595 - t449 * t601) * t558) * t485;
t422 = (-t446 * t470 + (t441 * t595 + t450 * t601) * t558) * t485;
t1 = [m(4) + (-t427 * t475 + (-t433 * t606 + t439 * t444) * t562) * t487 + (-t425 * t473 + (-t431 * t608 + t437 * t442) * t560) * t486 + (-t423 * t471 + (-t429 * t610 + t435 * t440) * t558) * t485, (-t427 * t474 + (t433 * t605 + t439 * t445) * t562) * t487 + (-t425 * t472 + (t431 * t607 + t437 * t443) * t560) * t486 + (-t423 * t470 + (t429 * t609 + t435 * t441) * t558) * t485, (t433 * t544 + t439 * t491) * t562 + (t431 * t542 + t437 * t490) * t560 + (t429 * t540 + t435 * t489) * t558; (-t426 * t475 + (-t432 * t606 + t438 * t444) * t562) * t487 + (-t424 * t473 + (-t430 * t608 + t436 * t442) * t560) * t486 + (-t422 * t471 + (-t428 * t610 + t434 * t440) * t558) * t485, m(4) + (-t426 * t474 + (t432 * t605 + t438 * t445) * t562) * t487 + (-t424 * t472 + (t430 * t607 + t436 * t443) * t560) * t486 + (-t422 * t470 + (t428 * t609 + t434 * t441) * t558) * t485, (t432 * t544 + t438 * t491) * t562 + (t430 * t542 + t436 * t490) * t560 + (t428 * t540 + t434 * t489) * t558; (-t457 * t475 + (t444 * t466 - t453 * t602) * t562) * t487 + (-t456 * t473 + (t442 * t465 - t451 * t603) * t560) * t486 + (-t455 * t471 + (t440 * t464 - t449 * t604) * t558) * t485, (-t457 * t474 + (t445 * t466 + t454 * t602) * t562) * t487 + (-t456 * t472 + (t443 * t465 + t452 * t603) * t560) * t486 + (-t455 * t470 + (t441 * t464 + t450 * t604) * t558) * t485, m(4) + (t460 * t544 + t466 * t491) * t562 + (t459 * t542 + t465 * t490) * t560 + (t458 * t540 + t464 * t489) * t558;];
MX  = t1;
