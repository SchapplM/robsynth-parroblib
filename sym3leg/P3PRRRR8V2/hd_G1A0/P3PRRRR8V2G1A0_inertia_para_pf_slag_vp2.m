% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G1A0
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2020-08-06 17:36
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:35:31
% EndTime: 2020-08-06 17:35:32
% DurationCPUTime: 1.16s
% Computational Cost: add. (3573->207), mult. (6735->393), div. (240->7), fcn. (6828->22), ass. (0->186)
t515 = cos(qJ(2,1));
t509 = sin(qJ(2,1));
t517 = pkin(7) + pkin(6);
t548 = t509 * t517;
t472 = pkin(2) * t515 + t548;
t496 = sin(pkin(8));
t498 = cos(pkin(8));
t478 = t517 * t515;
t469 = pkin(2) * t509 - t478;
t499 = cos(pkin(4));
t497 = sin(pkin(4));
t508 = sin(qJ(3,1));
t558 = t497 * t508;
t522 = pkin(3) * t558 - t469 * t499;
t590 = t472 * t498 + t522 * t496;
t513 = cos(qJ(2,2));
t507 = sin(qJ(2,2));
t550 = t507 * t517;
t471 = pkin(2) * t513 + t550;
t477 = t517 * t513;
t468 = pkin(2) * t507 - t477;
t506 = sin(qJ(3,2));
t560 = t497 * t506;
t523 = pkin(3) * t560 - t468 * t499;
t589 = t471 * t498 + t523 * t496;
t511 = cos(qJ(2,3));
t505 = sin(qJ(2,3));
t552 = t505 * t517;
t470 = pkin(2) * t511 + t552;
t476 = t517 * t511;
t467 = pkin(2) * t505 - t476;
t504 = sin(qJ(3,3));
t562 = t497 * t504;
t524 = pkin(3) * t562 - t467 * t499;
t588 = t470 * t498 + t524 * t496;
t510 = cos(qJ(3,3));
t493 = t510 ^ 2;
t587 = pkin(3) * t493;
t512 = cos(qJ(3,2));
t494 = t512 ^ 2;
t586 = pkin(3) * t494;
t514 = cos(qJ(3,1));
t495 = t514 ^ 2;
t585 = pkin(3) * t495;
t584 = m(3) * pkin(2) + mrSges(2,1);
t500 = legFrame(3,3);
t480 = sin(t500);
t483 = cos(t500);
t446 = -t496 * t480 + t483 * t498;
t449 = t498 * t480 + t483 * t496;
t473 = t510 * pkin(3) + pkin(2);
t464 = t505 * t473 - t476;
t568 = (t473 * t511 + t552) * t499;
t422 = -t464 * t446 - t449 * t568;
t547 = t510 * t497;
t553 = t504 * t499;
t434 = 0.1e1 / (t464 * t547 + t473 * t553);
t583 = t422 * t434;
t501 = legFrame(2,3);
t481 = sin(t501);
t484 = cos(t501);
t447 = -t496 * t481 + t484 * t498;
t450 = t498 * t481 + t484 * t496;
t474 = t512 * pkin(3) + pkin(2);
t465 = t507 * t474 - t477;
t567 = (t474 * t513 + t550) * t499;
t423 = -t465 * t447 - t450 * t567;
t546 = t512 * t497;
t551 = t506 * t499;
t435 = 0.1e1 / (t465 * t546 + t474 * t551);
t582 = t423 * t435;
t502 = legFrame(1,3);
t482 = sin(t502);
t485 = cos(t502);
t448 = -t496 * t482 + t485 * t498;
t451 = t498 * t482 + t485 * t496;
t475 = t514 * pkin(3) + pkin(2);
t466 = t509 * t475 - t478;
t566 = (t475 * t515 + t548) * t499;
t424 = -t466 * t448 - t451 * t566;
t545 = t514 * t497;
t549 = t508 * t499;
t436 = 0.1e1 / (t466 * t545 + t475 * t549);
t581 = t424 * t436;
t425 = -t446 * t568 + t464 * t449;
t580 = t425 * t434;
t426 = -t447 * t567 + t465 * t450;
t579 = t426 * t435;
t427 = -t448 * t566 + t466 * t451;
t578 = t427 * t436;
t561 = t497 * t505;
t428 = 0.1e1 / (t561 * t587 + (pkin(3) * t553 + t467 * t497) * t510 + pkin(2) * t553);
t491 = m(1) + m(2) + m(3);
t577 = t428 * t491;
t559 = t497 * t507;
t429 = 0.1e1 / (t559 * t586 + (pkin(3) * t551 + t468 * t497) * t512 + pkin(2) * t551);
t576 = t429 * t491;
t557 = t497 * t509;
t430 = 0.1e1 / (t557 * t585 + (pkin(3) * t549 + t469 * t497) * t514 + pkin(2) * t549);
t575 = t430 * t491;
t521 = 0.1e1 / pkin(3);
t574 = t434 * t521;
t573 = t435 * t521;
t572 = t436 * t521;
t479 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t528 = t510 * mrSges(3,1) - mrSges(3,2) * t504;
t571 = ((t528 + t584) * t511 + t505 * t479) * t497;
t527 = t512 * mrSges(3,1) - mrSges(3,2) * t506;
t570 = ((t527 + t584) * t513 + t507 * t479) * t497;
t526 = t514 * mrSges(3,1) - mrSges(3,2) * t508;
t569 = ((t526 + t584) * t515 + t509 * t479) * t497;
t556 = t499 * t505;
t555 = t499 * t507;
t554 = t499 * t509;
t544 = -0.2e1 * mrSges(3,2) * pkin(2);
t543 = pkin(2) * t562;
t542 = pkin(2) * t560;
t541 = pkin(2) * t558;
t540 = Ifges(3,3) * t574;
t539 = Ifges(3,3) * t573;
t538 = Ifges(3,3) * t572;
t537 = t428 * t571;
t536 = t429 * t570;
t535 = t430 * t569;
t440 = t528 * t499 - (t504 * mrSges(3,1) + t510 * mrSges(3,2)) * t561;
t534 = t440 * t574;
t486 = -pkin(6) * mrSges(3,2) + Ifges(3,6);
t487 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t458 = t486 * t510 - t504 * t487;
t533 = t458 * t574;
t441 = t527 * t499 - (t506 * mrSges(3,1) + t512 * mrSges(3,2)) * t559;
t532 = t441 * t573;
t459 = t486 * t512 - t506 * t487;
t531 = t459 * t573;
t442 = t526 * t499 - (t508 * mrSges(3,1) + t514 * mrSges(3,2)) * t557;
t530 = t442 * t572;
t460 = t486 * t514 - t508 * t487;
t529 = t460 * t572;
t431 = t496 * t470 - t524 * t498;
t452 = t496 * t556 - t498 * t511;
t455 = t496 * t511 + t498 * t556;
t398 = -(t452 * t483 + t480 * t455) * t587 + (-t431 * t480 + t588 * t483) * t510 + t449 * t543;
t416 = -t446 * t547 - (t446 * t556 + t511 * t449) * t504;
t380 = t398 * t577 + t416 * t537 + t425 * t534;
t432 = t496 * t471 - t523 * t498;
t453 = t496 * t555 - t498 * t513;
t456 = t496 * t513 + t498 * t555;
t399 = -(t453 * t484 + t481 * t456) * t586 + (-t432 * t481 + t589 * t484) * t512 + t450 * t542;
t417 = -t447 * t546 - (t447 * t555 + t513 * t450) * t506;
t381 = t399 * t576 + t417 * t536 + t426 * t532;
t433 = t496 * t472 - t522 * t498;
t454 = t496 * t554 - t498 * t515;
t457 = t496 * t515 + t498 * t554;
t400 = -(t454 * t485 + t482 * t457) * t585 + (-t433 * t482 + t590 * t485) * t514 + t451 * t541;
t418 = -t448 * t545 - (t448 * t554 + t515 * t451) * t508;
t382 = t400 * t575 + t418 * t535 + t427 * t530;
t401 = (-t480 * t452 + t455 * t483) * t587 + (t431 * t483 + t588 * t480) * t510 - t446 * t543;
t419 = -t449 * t547 - (-t511 * t446 + t449 * t556) * t504;
t383 = t401 * t577 + t419 * t537 + t422 * t534;
t402 = (-t481 * t453 + t456 * t484) * t586 + (t432 * t484 + t589 * t481) * t512 - t447 * t542;
t420 = -t450 * t546 - (-t513 * t447 + t450 * t555) * t506;
t384 = t402 * t576 + t420 * t536 + t423 * t532;
t403 = (-t482 * t454 + t457 * t485) * t585 + (t433 * t485 + t590 * t482) * t514 - t448 * t541;
t421 = -t451 * t545 - (-t515 * t448 + t451 * t554) * t508;
t385 = t403 * t575 + t421 * t535 + t424 * t530;
t525 = 0.2e1 * pkin(6) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3);
t516 = mrSges(3,1) * pkin(2);
t503 = -Ifges(3,1) + Ifges(3,2);
t439 = t503 * t495 + 0.2e1 * (Ifges(3,4) * t508 + t516) * t514 + t508 * t544 + t525;
t438 = t503 * t494 + 0.2e1 * (Ifges(3,4) * t506 + t516) * t512 + t506 * t544 + t525;
t437 = t503 * t493 + 0.2e1 * (Ifges(3,4) * t504 + t516) * t510 + t504 * t544 + t525;
t391 = t424 * t538 + (t403 * t442 + t421 * t460) * t430;
t390 = t423 * t539 + (t402 * t441 + t420 * t459) * t429;
t389 = t422 * t540 + (t401 * t440 + t419 * t458) * t428;
t388 = t427 * t538 + (t400 * t442 + t418 * t460) * t430;
t387 = t426 * t539 + (t399 * t441 + t417 * t459) * t429;
t386 = t425 * t540 + (t398 * t440 + t416 * t458) * t428;
t379 = t424 * t529 + (t403 * t569 + t421 * t439) * t430;
t378 = t423 * t531 + (t402 * t570 + t420 * t438) * t429;
t377 = t422 * t533 + (t401 * t571 + t419 * t437) * t428;
t376 = t427 * t529 + (t400 * t569 + t418 * t439) * t430;
t375 = t426 * t531 + (t399 * t570 + t417 * t438) * t429;
t374 = t425 * t533 + (t398 * t571 + t416 * t437) * t428;
t373 = t385 + t384 + t383;
t372 = t382 + t381 + t380;
t1 = [m(4) + (t376 * t418 + t382 * t400) * t430 + (t375 * t417 + t381 * t399) * t429 + (t374 * t416 + t380 * t398) * t428 + (t386 * t580 + t387 * t579 + t388 * t578) * t521, (t376 * t421 + t382 * t403) * t430 + (t375 * t420 + t381 * t402) * t429 + (t374 * t419 + t380 * t401) * t428 + (t386 * t583 + t387 * t582 + t388 * t581) * t521, t372; (t379 * t418 + t385 * t400) * t430 + (t378 * t417 + t384 * t399) * t429 + (t377 * t416 + t383 * t398) * t428 + (t389 * t580 + t390 * t579 + t391 * t578) * t521, m(4) + (t379 * t421 + t385 * t403) * t430 + (t378 * t420 + t384 * t402) * t429 + (t377 * t419 + t383 * t401) * t428 + (t389 * t583 + t390 * t582 + t391 * t581) * t521, t373; t372, t373, 0.3e1 * m(1) + 0.3e1 * m(2) + 0.3e1 * m(3) + m(4);];
MX  = t1;
