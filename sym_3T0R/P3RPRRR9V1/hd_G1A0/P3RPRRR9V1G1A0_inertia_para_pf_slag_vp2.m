% Calculate inertia matrix for parallel robot
% P3RPRRR9V1G1A0
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
% Datum: 2020-08-06 18:48
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR9V1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:47:06
% EndTime: 2020-08-06 18:47:07
% DurationCPUTime: 1.07s
% Computational Cost: add. (3099->229), mult. (3567->365), div. (270->7), fcn. (2130->32), ass. (0->152)
t497 = pkin(7) + qJ(3,3);
t474 = sin(t497);
t517 = sin(qJ(3,3));
t532 = 2 * pkin(7);
t574 = 2 * pkin(1);
t449 = pkin(3) * sin(0.2e1 * t497) + t474 * t574 + (sin(t532 + qJ(3,3)) + t517) * pkin(2);
t570 = t449 / 0.2e1;
t498 = pkin(7) + qJ(3,2);
t475 = sin(t498);
t519 = sin(qJ(3,2));
t450 = pkin(3) * sin((2 * t498)) + t475 * t574 + (sin((t532 + qJ(3,2))) + t519) * pkin(2);
t569 = t450 / 0.2e1;
t499 = pkin(7) + qJ(3,1);
t476 = sin(t499);
t521 = sin(qJ(3,1));
t451 = pkin(3) * sin((2 * t499)) + t476 * t574 + (sin((t532 + qJ(3,1))) + t521) * pkin(2);
t568 = t451 / 0.2e1;
t573 = 4 * Ifges(3,4);
t572 = 2 * mrSges(2,3) + 2 * mrSges(3,3);
t571 = pkin(2) * mrSges(3,1);
t529 = m(2) + m(3);
t567 = t529 / 0.2e1;
t507 = sin(pkin(7));
t566 = pkin(1) * t507;
t565 = pkin(1) * t529;
t477 = cos(t497);
t564 = pkin(3) * t477;
t478 = cos(t498);
t563 = pkin(3) * t478;
t479 = cos(t499);
t562 = pkin(3) * t479;
t515 = Ifges(3,1) - Ifges(3,2);
t561 = 2 * Ifges(2,4) - 2 * Ifges(3,4);
t514 = pkin(5) + qJ(2,1);
t513 = pkin(5) + qJ(2,2);
t512 = pkin(5) + qJ(2,3);
t560 = Ifges(3,4) * t517;
t559 = Ifges(3,4) * t519;
t558 = Ifges(3,4) * t521;
t533 = 0.1e1 / pkin(3);
t557 = Ifges(3,3) * t533;
t556 = m(3) * pkin(2) + mrSges(2,1);
t555 = t517 * mrSges(3,1) + mrSges(2,2);
t554 = t519 * mrSges(3,1) + mrSges(2,2);
t553 = t521 * mrSges(3,1) + mrSges(2,2);
t471 = 0.1e1 / t477;
t493 = -pkin(6) - t512;
t489 = 0.1e1 / t493;
t552 = t471 * t489;
t472 = 0.1e1 / t478;
t494 = -pkin(6) - t513;
t490 = 0.1e1 / t494;
t551 = t472 * t490;
t473 = 0.1e1 / t479;
t495 = -pkin(6) - t514;
t491 = 0.1e1 / t495;
t550 = t473 * t491;
t549 = t474 * t489;
t548 = t475 * t490;
t547 = t476 * t491;
t523 = cos(qJ(3,3));
t504 = t523 ^ 2;
t546 = t515 * t504;
t525 = cos(qJ(3,2));
t505 = t525 ^ 2;
t545 = t515 * t505;
t527 = cos(qJ(3,1));
t506 = t527 ^ 2;
t544 = t515 * t506;
t464 = t512 * mrSges(3,2) - Ifges(3,6);
t467 = t512 * mrSges(3,1) - Ifges(3,5);
t508 = cos(pkin(7));
t443 = (-t464 * t523 - t467 * t517) * t508 + t507 * (t464 * t517 - t467 * t523);
t543 = t533 * t443;
t465 = t513 * mrSges(3,2) - Ifges(3,6);
t468 = t513 * mrSges(3,1) - Ifges(3,5);
t444 = (-t465 * t525 - t468 * t519) * t508 + t507 * (t465 * t519 - t468 * t525);
t542 = t533 * t444;
t466 = t514 * mrSges(3,2) - Ifges(3,6);
t469 = t514 * mrSges(3,1) - Ifges(3,5);
t445 = (-t466 * t527 - t469 * t521) * t508 + t507 * (t466 * t521 - t469 * t527);
t541 = t533 * t445;
t500 = -0.2e1 * pkin(2) * mrSges(3,2);
t540 = -0.2e1 * t566;
t539 = mrSges(3,2) * t566;
t538 = -mrSges(3,2) * t517 + t556;
t537 = -mrSges(3,2) * t519 + t556;
t536 = -mrSges(3,2) * t521 + t556;
t535 = pkin(2) ^ 2 * m(3) - Ifges(2,1) + Ifges(2,2) + t515;
t534 = t529 * (pkin(1) ^ 2) + (2 * mrSges(3,3) * pkin(5)) + Ifges(2,1) + Ifges(3,2) + Ifges(1,3);
t528 = cos(qJ(1,1));
t526 = cos(qJ(1,2));
t524 = cos(qJ(1,3));
t522 = sin(qJ(1,1));
t520 = sin(qJ(1,2));
t518 = sin(qJ(1,3));
t511 = legFrame(1,3);
t510 = legFrame(2,3);
t509 = legFrame(3,3);
t503 = mrSges(3,1) * t574;
t501 = 0.2e1 * t571;
t496 = t508 ^ 2;
t485 = cos(t511);
t484 = cos(t510);
t483 = cos(t509);
t482 = sin(t511);
t481 = sin(t510);
t480 = sin(t509);
t470 = t508 * pkin(2) + pkin(1);
t463 = t482 * t528 + t485 * t522;
t462 = t481 * t526 + t484 * t520;
t461 = t480 * t524 + t483 * t518;
t460 = -t482 * t522 + t485 * t528;
t459 = -t481 * t520 + t484 * t526;
t458 = -t480 * t518 + t483 * t524;
t457 = t470 * t528 - t522 * t495;
t456 = t470 * t526 - t520 * t494;
t455 = t470 * t524 - t518 * t493;
t454 = t522 * t470 + t528 * t495;
t453 = t520 * t470 + t526 * t494;
t452 = t518 * t470 + t524 * t493;
t448 = (-mrSges(3,1) * t527 - t536) * t508 + (mrSges(3,2) * t527 + t553) * t507 - t565;
t447 = (-mrSges(3,1) * t525 - t537) * t508 + (mrSges(3,2) * t525 + t554) * t507 - t565;
t446 = (-mrSges(3,1) * t523 - t538) * t508 + (mrSges(3,2) * t523 + t555) * t507 - t565;
t442 = t454 * t485 + t457 * t482 + t463 * t562;
t441 = t453 * t484 + t456 * t481 + t462 * t563;
t440 = t452 * t483 + t455 * t480 + t461 * t564;
t439 = -t454 * t482 + t457 * t485 + t460 * t562;
t438 = -t453 * t481 + t456 * t484 + t459 * t563;
t437 = -t452 * t480 + t455 * t483 + t458 * t564;
t436 = (t476 * t448 + t451 * t567) * t550;
t435 = (t475 * t447 + t450 * t567) * t551;
t434 = (t474 * t446 + t449 * t567) * t552;
t433 = (t442 * t529 + t448 * t463) * t491;
t432 = (t441 * t529 + t447 * t462) * t490;
t431 = (t440 * t529 + t446 * t461) * t489;
t430 = (t439 * t529 + t448 * t460) * t491;
t429 = (t438 * t529 + t447 * t459) * t490;
t428 = (t437 * t529 + t446 * t458) * t489;
t427 = (-0.2e1 * t544 + (t501 + 0.4e1 * t558) * t527 + t521 * t500 + t535) * t496 + (t503 * t527 + t536 * t574 + (t506 * t573 + t500 * t527 + 0.2e1 * (t515 * t527 - t571) * t521 + t561) * t507) * t508 + t544 + 0.2e1 * (-t539 - t558) * t527 + t553 * t540 + t514 ^ 2 * m(3) + t534 + ((m(2) * qJ(2,1) + t572) * qJ(2,1));
t426 = (-0.2e1 * t545 + (t501 + 0.4e1 * t559) * t525 + t519 * t500 + t535) * t496 + (t503 * t525 + t537 * t574 + (t505 * t573 + t500 * t525 + 0.2e1 * (t515 * t525 - t571) * t519 + t561) * t507) * t508 + t545 + 0.2e1 * (-t539 - t559) * t525 + t554 * t540 + t513 ^ 2 * m(3) + t534 + ((m(2) * qJ(2,2) + t572) * qJ(2,2));
t425 = (-0.2e1 * t546 + (t501 + 0.4e1 * t560) * t523 + t517 * t500 + t535) * t496 + (t503 * t523 + t538 * t574 + (t504 * t573 + t500 * t523 + 0.2e1 * (t515 * t523 - t571) * t517 + t561) * t507) * t508 + t546 + 0.2e1 * (-t539 - t560) * t523 + t555 * t540 + t512 ^ 2 * m(3) + t534 + ((m(2) * qJ(2,3) + t572) * qJ(2,3));
t424 = (t427 * t463 + t442 * t448) * t491;
t423 = (t426 * t462 + t441 * t447) * t490;
t422 = (t425 * t461 + t440 * t446) * t489;
t421 = (t427 * t460 + t439 * t448) * t491;
t420 = (t426 * t459 + t438 * t447) * t490;
t419 = (t425 * t458 + t437 * t446) * t489;
t418 = (t541 - (t476 * t427 + t448 * t568) * t491) * t473;
t417 = (t542 - (t475 * t426 + t447 * t569) * t490) * t472;
t416 = (t543 - (t474 * t425 + t446 * t570) * t489) * t471;
t1 = [m(4) - (-t421 * t460 - t430 * t439) * t491 - (-t420 * t459 - t429 * t438) * t490 - (-t419 * t458 - t428 * t437) * t489, -(-t421 * t463 - t430 * t442) * t491 - (-t420 * t462 - t429 * t441) * t490 - (-t419 * t461 - t428 * t440) * t489, -(-t421 * t476 - t430 * t568 + t460 * t541) * t550 - (-t420 * t475 - t429 * t569 + t459 * t542) * t551 - (-t419 * t474 - t428 * t570 + t458 * t543) * t552; -(-t424 * t460 - t433 * t439) * t491 - (-t423 * t459 - t432 * t438) * t490 - (-t422 * t458 - t431 * t437) * t489, m(4) - (-t424 * t463 - t433 * t442) * t491 - (-t423 * t462 - t432 * t441) * t490 - (-t422 * t461 - t431 * t440) * t489, -(-t424 * t476 - t433 * t568 + t463 * t541) * t550 - (-t423 * t475 - t432 * t569 + t462 * t542) * t551 - (-t422 * t474 - t431 * t570 + t461 * t543) * t552; -(t418 * t460 - t436 * t439) * t491 - (t417 * t459 - t435 * t438) * t490 - (t416 * t458 - t434 * t437) * t489, -(t418 * t463 - t436 * t442) * t491 - (t417 * t462 - t435 * t441) * t490 - (t416 * t461 - t434 * t440) * t489, m(4) + (-t418 * t547 + t436 * t491 * t568 + (-t445 * t547 + t557) * t533 * t473) * t473 + (-t417 * t548 + t435 * t490 * t569 + (-t444 * t548 + t557) * t533 * t472) * t472 + (-t416 * t549 + t434 * t489 * t570 + (-t443 * t549 + t557) * t533 * t471) * t471;];
MX  = t1;
