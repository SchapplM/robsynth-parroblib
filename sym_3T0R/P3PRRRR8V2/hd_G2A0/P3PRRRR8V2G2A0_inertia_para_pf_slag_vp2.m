% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G2A0
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
% Datum: 2020-08-06 17:50
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:49:28
% EndTime: 2020-08-06 17:49:29
% DurationCPUTime: 1.53s
% Computational Cost: add. (5178->243), mult. (10728->483), div. (432->4), fcn. (10116->22), ass. (0->188)
t489 = sin(qJ(3,3));
t579 = pkin(2) * t489;
t491 = sin(qJ(3,2));
t578 = pkin(2) * t491;
t493 = sin(qJ(3,1));
t577 = pkin(2) * t493;
t495 = cos(qJ(3,3));
t478 = t495 ^ 2;
t576 = pkin(3) * t478;
t497 = cos(qJ(3,2));
t479 = t497 ^ 2;
t575 = pkin(3) * t479;
t499 = cos(qJ(3,1));
t480 = t499 ^ 2;
t574 = pkin(3) * t480;
t573 = pkin(3) * t495;
t572 = pkin(3) * t497;
t571 = pkin(3) * t499;
t570 = m(3) * pkin(2) + mrSges(2,1);
t490 = sin(qJ(2,3));
t496 = cos(qJ(2,3));
t502 = pkin(7) + pkin(6);
t458 = pkin(2) * t490 - t502 * t496;
t461 = pkin(2) * t496 + t490 * t502;
t481 = sin(pkin(8));
t483 = cos(pkin(8));
t484 = cos(pkin(4));
t541 = t484 * t496;
t548 = t483 * t484;
t419 = (t481 * t490 - t483 * t541) * t573 - t461 * t548 + t481 * t458;
t506 = 0.1e1 / pkin(3);
t569 = t419 * t506;
t492 = sin(qJ(2,2));
t498 = cos(qJ(2,2));
t459 = pkin(2) * t492 - t502 * t498;
t462 = pkin(2) * t498 + t492 * t502;
t540 = t484 * t498;
t420 = (t481 * t492 - t483 * t540) * t572 - t462 * t548 + t481 * t459;
t568 = t420 * t506;
t494 = sin(qJ(2,1));
t500 = cos(qJ(2,1));
t460 = pkin(2) * t494 - t502 * t500;
t463 = pkin(2) * t500 + t494 * t502;
t539 = t484 * t500;
t421 = (t481 * t494 - t483 * t539) * t571 - t463 * t548 + t481 * t460;
t567 = t421 * t506;
t553 = t481 * t484;
t422 = (t481 * t541 + t483 * t490) * t573 + t461 * t553 + t458 * t483;
t566 = t422 * t506;
t423 = (t481 * t540 + t483 * t492) * t572 + t462 * t553 + t459 * t483;
t565 = t423 * t506;
t424 = (t481 * t539 + t483 * t494) * t571 + t463 * t553 + t460 * t483;
t564 = t424 * t506;
t531 = t495 * mrSges(3,1) - mrSges(3,2) * t489;
t482 = sin(pkin(4));
t537 = t490 * t482;
t440 = t531 * t484 - (t489 * mrSges(3,1) + t495 * mrSges(3,2)) * t537;
t563 = t440 * t506;
t530 = t497 * mrSges(3,1) - mrSges(3,2) * t491;
t535 = t492 * t482;
t441 = t530 * t484 - (t491 * mrSges(3,1) + t497 * mrSges(3,2)) * t535;
t562 = t441 * t506;
t529 = t499 * mrSges(3,1) - mrSges(3,2) * t493;
t533 = t494 * t482;
t442 = t529 * t484 - (t493 * mrSges(3,1) + t499 * mrSges(3,2)) * t533;
t561 = t442 * t506;
t464 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t560 = ((t531 + t570) * t496 + t490 * t464) * t482;
t559 = ((t530 + t570) * t498 + t492 * t464) * t482;
t558 = ((t529 + t570) * t500 + t494 * t464) * t482;
t465 = -mrSges(3,2) * pkin(6) + Ifges(3,6);
t466 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t455 = t465 * t495 - t489 * t466;
t557 = t455 * t506;
t456 = t465 * t497 - t491 * t466;
t556 = t456 * t506;
t457 = t465 * t499 - t493 * t466;
t555 = t457 * t506;
t554 = t481 * t482;
t552 = t482 * t483;
t551 = t482 * t489;
t550 = t482 * t491;
t549 = t482 * t493;
t547 = t483 * t495;
t546 = t483 * t497;
t545 = t483 * t499;
t544 = t484 * t490;
t543 = t484 * t492;
t542 = t484 * t494;
t538 = t489 * t484;
t536 = t491 * t484;
t534 = t493 * t484;
t532 = -0.2e1 * pkin(2) * mrSges(3,2);
t528 = 0.2e1 * pkin(6) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3);
t527 = pkin(3) * t551 - t458 * t484;
t526 = pkin(3) * t550 - t459 * t484;
t525 = pkin(3) * t549 - t460 * t484;
t449 = t481 * t544 - t483 * t496;
t434 = t489 * t449 + t495 * t554;
t524 = Ifges(3,3) * t566 + t434 * t455;
t450 = t481 * t543 - t483 * t498;
t435 = t491 * t450 + t497 * t554;
t523 = Ifges(3,3) * t565 + t435 * t456;
t451 = t481 * t542 - t483 * t500;
t436 = t493 * t451 + t499 * t554;
t522 = Ifges(3,3) * t564 + t436 * t457;
t428 = -t481 * t461 + t527 * t483;
t446 = pkin(3) * t538 + t458 * t482;
t452 = t481 * t496 + t483 * t544;
t486 = legFrame(3,2);
t470 = sin(t486);
t473 = cos(t486);
t410 = -(t452 * t470 - t473 * t537) * t576 + (t428 * t470 + t473 * t446) * t495 + (t470 * t552 + t473 * t484) * t579;
t425 = 0.1e1 / (pkin(2) * t538 + t446 * t495 + t537 * t576);
t485 = -Ifges(3,1) + Ifges(3,2);
t501 = mrSges(3,1) * pkin(2);
t431 = t485 * t478 + 0.2e1 * (Ifges(3,4) * t489 + t501) * t495 + t489 * t532 + t528;
t512 = t422 * t557 + t431 * t434;
t383 = (t410 * t560 + t512 * t470) * t425;
t395 = (t410 * t440 + t524 * t470) * t425;
t521 = t383 * t434 + t395 * t566;
t411 = (t452 * t473 + t470 * t537) * t576 + (-t428 * t473 + t446 * t470) * t495 + (t484 * t470 - t473 * t552) * t579;
t384 = (t411 * t560 - t512 * t473) * t425;
t396 = (t411 * t440 - t524 * t473) * t425;
t520 = t384 * t434 + t396 * t566;
t429 = -t481 * t462 + t526 * t483;
t447 = pkin(3) * t536 + t459 * t482;
t453 = t481 * t498 + t483 * t543;
t487 = legFrame(2,2);
t471 = sin(t487);
t474 = cos(t487);
t412 = -(t453 * t471 - t474 * t535) * t575 + (t429 * t471 + t474 * t447) * t497 + (t471 * t552 + t474 * t484) * t578;
t426 = 0.1e1 / (pkin(2) * t536 + t447 * t497 + t535 * t575);
t432 = t485 * t479 + 0.2e1 * (Ifges(3,4) * t491 + t501) * t497 + t491 * t532 + t528;
t511 = t423 * t556 + t432 * t435;
t385 = (t412 * t559 + t511 * t471) * t426;
t397 = (t412 * t441 + t523 * t471) * t426;
t519 = t385 * t435 + t397 * t565;
t413 = (t453 * t474 + t471 * t535) * t575 + (-t429 * t474 + t447 * t471) * t497 + (t484 * t471 - t474 * t552) * t578;
t386 = (t413 * t559 - t511 * t474) * t426;
t398 = (t413 * t441 - t523 * t474) * t426;
t518 = t386 * t435 + t398 * t565;
t430 = -t481 * t463 + t525 * t483;
t448 = pkin(3) * t534 + t460 * t482;
t454 = t481 * t500 + t483 * t542;
t488 = legFrame(1,2);
t472 = sin(t488);
t475 = cos(t488);
t414 = -(t454 * t472 - t475 * t533) * t574 + (t430 * t472 + t475 * t448) * t499 + (t472 * t552 + t475 * t484) * t577;
t427 = 0.1e1 / (pkin(2) * t534 + t448 * t499 + t533 * t574);
t433 = t485 * t480 + 0.2e1 * (Ifges(3,4) * t493 + t501) * t499 + t493 * t532 + t528;
t510 = t424 * t555 + t433 * t436;
t387 = (t414 * t558 + t510 * t472) * t427;
t399 = (t414 * t442 + t522 * t472) * t427;
t517 = t387 * t436 + t399 * t564;
t415 = (t454 * t475 + t472 * t533) * t574 + (-t430 * t475 + t448 * t472) * t499 + (t484 * t472 - t475 * t552) * t577;
t388 = (t415 * t558 - t510 * t475) * t427;
t400 = (t415 * t442 - t522 * t475) * t427;
t516 = t388 * t436 + t400 * t564;
t416 = -t449 * t576 + t461 * t547 + (pkin(2) * t551 + t527 * t495) * t481;
t437 = -t489 * t452 - t482 * t547;
t401 = (t416 * t560 + t419 * t557 + t431 * t437) * t425;
t407 = (Ifges(3,3) * t569 + t416 * t440 + t437 * t455) * t425;
t515 = t401 * t434 + t407 * t566;
t417 = -t450 * t575 + t462 * t546 + (pkin(2) * t550 + t526 * t497) * t481;
t438 = -t491 * t453 - t482 * t546;
t402 = (t417 * t559 + t420 * t556 + t432 * t438) * t426;
t408 = (Ifges(3,3) * t568 + t417 * t441 + t438 * t456) * t426;
t514 = t402 * t435 + t408 * t565;
t418 = -t451 * t574 + t463 * t545 + (pkin(2) * t549 + t525 * t499) * t481;
t439 = -t493 * t454 - t482 * t545;
t403 = (t418 * t558 + t421 * t555 + t433 * t439) * t427;
t409 = (Ifges(3,3) * t567 + t418 * t442 + t439 * t457) * t427;
t513 = t403 * t436 + t409 * t564;
t509 = t422 * t563 + t434 * t560;
t508 = t423 * t562 + t435 * t559;
t507 = t424 * t561 + t436 * t558;
t476 = m(1) + m(2) + m(3);
t406 = (t418 * t476 + t421 * t561 + t439 * t558) * t427;
t405 = (t417 * t476 + t420 * t562 + t438 * t559) * t426;
t404 = (t416 * t476 + t419 * t563 + t437 * t560) * t425;
t394 = (t415 * t476 - t507 * t475) * t427;
t393 = (t414 * t476 + t507 * t472) * t427;
t392 = (t413 * t476 - t508 * t474) * t426;
t391 = (t412 * t476 + t508 * t471) * t426;
t390 = (t411 * t476 - t509 * t473) * t425;
t389 = (t410 * t476 + t509 * t470) * t425;
t1 = [m(4) + (t394 * t415 - t516 * t475) * t427 + (t392 * t413 - t518 * t474) * t426 + (t390 * t411 - t520 * t473) * t425, (t394 * t414 + t516 * t472) * t427 + (t392 * t412 + t518 * t471) * t426 + (t390 * t410 + t520 * t470) * t425, (t388 * t439 + t394 * t418 + t400 * t567) * t427 + (t386 * t438 + t392 * t417 + t398 * t568) * t426 + (t384 * t437 + t390 * t416 + t396 * t569) * t425; (t393 * t415 - t517 * t475) * t427 + (t391 * t413 - t519 * t474) * t426 + (t389 * t411 - t521 * t473) * t425, m(4) + (t393 * t414 + t517 * t472) * t427 + (t391 * t412 + t519 * t471) * t426 + (t389 * t410 + t521 * t470) * t425, (t387 * t439 + t393 * t418 + t399 * t567) * t427 + (t385 * t438 + t391 * t417 + t397 * t568) * t426 + (t383 * t437 + t389 * t416 + t395 * t569) * t425; (t406 * t415 - t513 * t475) * t427 + (t405 * t413 - t514 * t474) * t426 + (t404 * t411 - t515 * t473) * t425, (t406 * t414 + t513 * t472) * t427 + (t405 * t412 + t514 * t471) * t426 + (t404 * t410 + t515 * t470) * t425, m(4) + (t403 * t439 + t406 * t418 + t409 * t567) * t427 + (t402 * t438 + t405 * t417 + t408 * t568) * t426 + (t401 * t437 + t404 * t416 + t407 * t569) * t425;];
MX  = t1;
