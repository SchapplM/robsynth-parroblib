% Calculate inertia matrix for parallel robot
% P3RRPRR8V2G1A0
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
% Datum: 2022-11-07 13:12
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V2G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-07 13:11:15
% EndTime: 2022-11-07 13:11:15
% DurationCPUTime: 0.94s
% Computational Cost: add. (3183->225), mult. (4212->341), div. (249->6), fcn. (2325->35), ass. (0->142)
t568 = 2 * pkin(1);
t594 = pkin(3) ^ 2 / 0.2e1;
t560 = pkin(2) ^ 2;
t593 = t560 / 0.2e1;
t592 = t568 / 0.2e1;
t591 = pkin(2) * pkin(3);
t590 = 2 * mrSges(3,1);
t589 = 2 * mrSges(3,3);
t534 = cos(pkin(7));
t588 = 0.2e1 * t534;
t587 = m(3) * pkin(1);
t586 = m(3) * pkin(2);
t585 = mrSges(3,2) * pkin(1);
t529 = qJ(2,3) + pkin(7);
t542 = sin(qJ(2,3));
t487 = t542 * pkin(2) + pkin(3) * sin(t529);
t556 = 0.2e1 * qJ(2,3);
t584 = sin(t556 + pkin(7)) * t591 + sin(0.2e1 * t529) * t594 + sin(t556) * t593 + t487 * t592;
t530 = qJ(2,2) + pkin(7);
t544 = sin(qJ(2,2));
t488 = t544 * pkin(2) + pkin(3) * sin(t530);
t557 = 0.2e1 * qJ(2,2);
t583 = sin(t557 + pkin(7)) * t591 + sin(0.2e1 * t530) * t594 + sin(t557) * t593 + t488 * t592;
t531 = qJ(2,1) + pkin(7);
t546 = sin(qJ(2,1));
t489 = t546 * pkin(2) + pkin(3) * sin(t531);
t558 = 0.2e1 * qJ(2,1);
t582 = sin(t558 + pkin(7)) * t591 + sin(0.2e1 * t531) * t594 + sin(t558) * t593 + t489 * t592;
t581 = pkin(3) * cos(t529);
t580 = pkin(3) * cos(t530);
t579 = pkin(3) * cos(t531);
t541 = Ifges(3,2) - Ifges(3,1);
t540 = (qJ(3,1) + pkin(5));
t539 = (qJ(3,2) + pkin(5));
t538 = (qJ(3,3) + pkin(5));
t533 = sin(pkin(7));
t578 = mrSges(3,2) * t533;
t577 = Ifges(3,4) * t533;
t576 = t533 * mrSges(3,1);
t575 = Ifges(2,5) + (-mrSges(2,1) - t586) * pkin(5);
t548 = cos(qJ(2,3));
t521 = t548 * pkin(2);
t483 = 0.1e1 / (t521 + t581);
t525 = -pkin(6) - t538;
t518 = 0.1e1 / t525;
t574 = t483 * t518;
t550 = cos(qJ(2,2));
t522 = t550 * pkin(2);
t484 = 0.1e1 / (t522 + t580);
t526 = -pkin(6) - t539;
t519 = 0.1e1 / t526;
t573 = t484 * t519;
t552 = cos(qJ(2,1));
t523 = t552 * pkin(2);
t485 = 0.1e1 / (t523 + t579);
t527 = -pkin(6) - t540;
t520 = 0.1e1 / t527;
t572 = t485 * t520;
t528 = t534 ^ 2;
t499 = t541 * t528;
t570 = m(3) * t560 - 0.2e1 * pkin(2) * t578;
t569 = pkin(2) * t590;
t567 = -0.2e1 * pkin(1) * (mrSges(2,2) + t576);
t566 = -pkin(5) * mrSges(2,2) + Ifges(2,6);
t564 = -t578 + t586;
t565 = pkin(1) * t590 * t534 + (mrSges(2,1) + t564) * t568;
t563 = mrSges(3,2) * t534 + t576;
t562 = 0.4e1 * Ifges(3,4) * t528 + 0.2e1 * (-pkin(2) * mrSges(3,2) - t533 * t541) * t534 - 0.2e1 * pkin(2) * t576 + (2 * Ifges(2,4)) - 0.2e1 * Ifges(3,4);
t561 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + 0.2e1 * (mrSges(2,3) + mrSges(3,3)) * pkin(5) + ((m(2) + m(3)) * pkin(1) ^ 2) + m(2) * pkin(5) ^ 2 - t499;
t553 = cos(qJ(1,1));
t551 = cos(qJ(1,2));
t549 = cos(qJ(1,3));
t547 = sin(qJ(1,1));
t545 = sin(qJ(1,2));
t543 = sin(qJ(1,3));
t537 = legFrame(1,3);
t536 = legFrame(2,3);
t535 = legFrame(3,3);
t517 = cos(t537);
t516 = cos(t536);
t515 = cos(t535);
t514 = sin(t537);
t513 = sin(t536);
t512 = sin(t535);
t504 = t523 + pkin(1);
t503 = t522 + pkin(1);
t502 = t521 + pkin(1);
t497 = -t540 * mrSges(3,1) + Ifges(3,5);
t496 = -t539 * mrSges(3,1) + Ifges(3,5);
t495 = -t538 * mrSges(3,1) + Ifges(3,5);
t494 = t540 * mrSges(3,2) - Ifges(3,6);
t493 = t539 * mrSges(3,2) - Ifges(3,6);
t492 = t538 * mrSges(3,2) - Ifges(3,6);
t482 = -mrSges(3,1) * t534 - t564;
t481 = -t514 * t547 + t517 * t553;
t480 = t514 * t553 + t517 * t547;
t479 = -t513 * t545 + t516 * t551;
t478 = t513 * t551 + t516 * t545;
t477 = -t512 * t543 + t515 * t549;
t476 = t512 * t549 + t515 * t543;
t475 = t534 * t569 + Ifges(2,3) + Ifges(3,3) + t570;
t474 = t504 * t553 - t547 * t527;
t473 = t503 * t551 - t545 * t526;
t472 = t502 * t549 - t543 * t525;
t471 = t547 * t504 + t553 * t527;
t470 = t545 * t503 + t551 * t526;
t469 = t543 * t502 + t549 * t525;
t468 = 0.2e1 * t499 + (t569 + 0.4e1 * t577) * t534 + Ifges(2,2) - Ifges(2,1) + t570 - t541;
t467 = t482 * t552 + t563 * t546 - t587;
t466 = t482 * t550 + t563 * t544 - t587;
t465 = t482 * t548 + t563 * t542 - t587;
t461 = -t514 * t471 + t474 * t517 + t481 * t579;
t460 = t471 * t517 + t474 * t514 + t480 * t579;
t459 = -t513 * t470 + t473 * t516 + t479 * t580;
t458 = t470 * t516 + t473 * t513 + t478 * t580;
t457 = -t512 * t469 + t472 * t515 + t477 * t581;
t456 = t469 * t515 + t472 * t512 + t476 * t581;
t455 = (t494 * t533 + t497 * t534 + (-m(3) * qJ(3,1) - mrSges(3,3)) * pkin(2) + t575) * t546 + t552 * (-t494 * t534 + t497 * t533 + t566);
t454 = (t493 * t533 + t496 * t534 + (-m(3) * qJ(3,2) - mrSges(3,3)) * pkin(2) + t575) * t544 + t550 * (-t493 * t534 + t496 * t533 + t566);
t453 = (t492 * t533 + t495 * t534 + (-m(3) * qJ(3,3) - mrSges(3,3)) * pkin(2) + t575) * t542 + t548 * (-t492 * t534 + t495 * t533 + t566);
t452 = (m(3) * t582 + t489 * t467) * t572;
t451 = (m(3) * t583 + t488 * t466) * t573;
t450 = (m(3) * t584 + t487 * t465) * t574;
t449 = (m(3) * t461 + t467 * t481) * t520;
t448 = (m(3) * t460 + t467 * t480) * t520;
t447 = (m(3) * t459 + t466 * t479) * t519;
t446 = (m(3) * t458 + t466 * t478) * t519;
t445 = (m(3) * t457 + t465 * t477) * t518;
t444 = (m(3) * t456 + t465 * t476) * t518;
t443 = (-t546 * t585 - t577) * t588 + t546 * t567 + (t540 ^ 2 * m(3)) + (qJ(3,1) * t589) + t561 + (t468 * t552 + t562 * t546 + t565) * t552;
t442 = (-t544 * t585 - t577) * t588 + t544 * t567 + (t539 ^ 2 * m(3)) + (qJ(3,2) * t589) + t561 + (t468 * t550 + t562 * t544 + t565) * t550;
t441 = (-t542 * t585 - t577) * t588 + t542 * t567 + (t538 ^ 2 * m(3)) + (qJ(3,3) * t589) + t561 + (t468 * t548 + t562 * t542 + t565) * t548;
t440 = (t443 * t481 + t461 * t467) * t520;
t439 = (t443 * t480 + t460 * t467) * t520;
t438 = (t442 * t479 + t459 * t466) * t519;
t437 = (t442 * t478 + t458 * t466) * t519;
t436 = (t441 * t477 + t457 * t465) * t518;
t435 = (t441 * t476 + t456 * t465) * t518;
t434 = (t455 - (t489 * t443 + t467 * t582) * t520) * t485;
t433 = (t454 - (t488 * t442 + t466 * t583) * t519) * t484;
t432 = (t453 - (t487 * t441 + t465 * t584) * t518) * t483;
t1 = [m(4) - (-t440 * t481 - t449 * t461) * t520 - (-t438 * t479 - t447 * t459) * t519 - (-t436 * t477 - t445 * t457) * t518, -(-t440 * t480 - t449 * t460) * t520 - (-t438 * t478 - t447 * t458) * t519 - (-t436 * t476 - t445 * t456) * t518, -(-t440 * t489 - t449 * t582 + t481 * t455) * t572 - (-t438 * t488 - t447 * t583 + t479 * t454) * t573 - (-t436 * t487 - t445 * t584 + t477 * t453) * t574; -(-t439 * t481 - t448 * t461) * t520 - (-t437 * t479 - t446 * t459) * t519 - (-t435 * t477 - t444 * t457) * t518, m(4) - (-t439 * t480 - t448 * t460) * t520 - (-t437 * t478 - t446 * t458) * t519 - (-t435 * t476 - t444 * t456) * t518, -(-t439 * t489 - t448 * t582 + t480 * t455) * t572 - (-t437 * t488 - t446 * t583 + t478 * t454) * t573 - (-t435 * t487 - t444 * t584 + t476 * t453) * t574; -(t434 * t481 - t452 * t461) * t520 - (t433 * t479 - t451 * t459) * t519 - (t432 * t477 - t450 * t457) * t518, -(t434 * t480 - t452 * t460) * t520 - (t433 * t478 - t451 * t458) * t519 - (t432 * t476 - t450 * t456) * t518, m(4) + (t485 * t475 - (-t452 * t582 + (t485 * t455 + t434) * t489) * t520) * t485 + (t484 * t475 - (-t451 * t583 + (t484 * t454 + t433) * t488) * t519) * t484 + (t483 * t475 - (-t450 * t584 + (t483 * t453 + t432) * t487) * t518) * t483;];
MX  = t1;
