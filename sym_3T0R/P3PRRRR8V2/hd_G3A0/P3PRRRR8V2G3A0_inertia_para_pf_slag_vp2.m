% Calculate inertia matrix for parallel robot
% P3PRRRR8V2G3A0
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
% Datum: 2020-08-06 18:05
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(8,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V2G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:04:34
% EndTime: 2020-08-06 18:04:36
% DurationCPUTime: 1.50s
% Computational Cost: add. (5178->243), mult. (10728->485), div. (432->4), fcn. (10116->22), ass. (0->187)
t491 = sin(qJ(3,3));
t580 = pkin(2) * t491;
t493 = sin(qJ(3,2));
t579 = pkin(2) * t493;
t495 = sin(qJ(3,1));
t578 = pkin(2) * t495;
t497 = cos(qJ(3,3));
t480 = t497 ^ 2;
t577 = pkin(3) * t480;
t499 = cos(qJ(3,2));
t481 = t499 ^ 2;
t576 = pkin(3) * t481;
t501 = cos(qJ(3,1));
t482 = t501 ^ 2;
t575 = pkin(3) * t482;
t574 = pkin(3) * t497;
t573 = pkin(3) * t499;
t572 = pkin(3) * t501;
t571 = m(3) * pkin(2) + mrSges(2,1);
t492 = sin(qJ(2,3));
t498 = cos(qJ(2,3));
t504 = pkin(7) + pkin(6);
t460 = pkin(2) * t492 - t504 * t498;
t463 = pkin(2) * t498 + t492 * t504;
t483 = sin(pkin(8));
t485 = cos(pkin(8));
t486 = cos(pkin(4));
t537 = t486 * t498;
t544 = t485 * t486;
t421 = (t483 * t492 - t485 * t537) * t574 - t463 * t544 + t460 * t483;
t508 = 0.1e1 / pkin(3);
t570 = t421 * t508;
t494 = sin(qJ(2,2));
t500 = cos(qJ(2,2));
t461 = pkin(2) * t494 - t504 * t500;
t464 = pkin(2) * t500 + t494 * t504;
t536 = t486 * t500;
t422 = (t483 * t494 - t485 * t536) * t573 - t464 * t544 + t461 * t483;
t569 = t422 * t508;
t496 = sin(qJ(2,1));
t502 = cos(qJ(2,1));
t462 = pkin(2) * t496 - t504 * t502;
t465 = pkin(2) * t502 + t496 * t504;
t535 = t486 * t502;
t423 = (t483 * t496 - t485 * t535) * t572 - t465 * t544 + t462 * t483;
t568 = t423 * t508;
t554 = t483 * t486;
t424 = (t483 * t537 + t485 * t492) * t574 + t463 * t554 + t460 * t485;
t567 = t424 * t508;
t425 = (t483 * t536 + t485 * t494) * t573 + t464 * t554 + t461 * t485;
t566 = t425 * t508;
t426 = (t483 * t535 + t485 * t496) * t572 + t465 * t554 + t462 * t485;
t565 = t426 * t508;
t533 = t497 * mrSges(3,1) - mrSges(3,2) * t491;
t484 = sin(pkin(4));
t552 = t484 * t492;
t442 = t533 * t486 - (t491 * mrSges(3,1) + t497 * mrSges(3,2)) * t552;
t564 = t442 * t508;
t532 = t499 * mrSges(3,1) - mrSges(3,2) * t493;
t550 = t484 * t494;
t443 = t532 * t486 - (t493 * mrSges(3,1) + t499 * mrSges(3,2)) * t550;
t563 = t443 * t508;
t531 = t501 * mrSges(3,1) - mrSges(3,2) * t495;
t548 = t484 * t496;
t444 = t531 * t486 - (t495 * mrSges(3,1) + t501 * mrSges(3,2)) * t548;
t562 = t444 * t508;
t466 = m(3) * pkin(6) - mrSges(2,2) + mrSges(3,3);
t561 = ((t533 + t571) * t498 + t466 * t492) * t484;
t560 = ((t532 + t571) * t500 + t466 * t494) * t484;
t559 = ((t531 + t571) * t502 + t466 * t496) * t484;
t467 = -pkin(6) * mrSges(3,2) + Ifges(3,6);
t468 = mrSges(3,1) * pkin(6) - Ifges(3,5);
t457 = t467 * t497 - t491 * t468;
t558 = t457 * t508;
t458 = t467 * t499 - t493 * t468;
t557 = t458 * t508;
t459 = t467 * t501 - t495 * t468;
t556 = t459 * t508;
t555 = t483 * t484;
t553 = t484 * t491;
t551 = t484 * t493;
t549 = t484 * t495;
t547 = t484 * t497;
t546 = t484 * t499;
t545 = t484 * t501;
t543 = t486 * t491;
t542 = t486 * t492;
t541 = t486 * t493;
t540 = t486 * t494;
t539 = t486 * t495;
t538 = t486 * t496;
t534 = -0.2e1 * mrSges(3,2) * pkin(2);
t530 = 0.2e1 * pkin(6) * mrSges(3,3) + Ifges(3,1) + Ifges(2,3) + (pkin(2) ^ 2 + pkin(6) ^ 2) * m(3);
t529 = pkin(3) * t553 - t460 * t486;
t528 = pkin(3) * t551 - t461 * t486;
t527 = pkin(3) * t549 - t462 * t486;
t454 = t483 * t498 + t485 * t542;
t439 = t491 * t454 + t485 * t547;
t526 = Ifges(3,3) * t570 - t439 * t457;
t455 = t483 * t500 + t485 * t540;
t440 = t493 * t455 + t485 * t546;
t525 = Ifges(3,3) * t569 - t440 * t458;
t456 = t483 * t502 + t485 * t538;
t441 = t495 * t456 + t485 * t545;
t524 = Ifges(3,3) * t568 - t441 * t459;
t430 = t485 * t463 + t529 * t483;
t448 = pkin(3) * t543 + t484 * t460;
t451 = t483 * t542 - t485 * t498;
t488 = legFrame(3,2);
t472 = sin(t488);
t475 = cos(t488);
t412 = -(t451 * t475 - t472 * t552) * t577 + (t430 * t475 + t472 * t448) * t497 + (t486 * t472 + t475 * t555) * t580;
t427 = 0.1e1 / (pkin(2) * t543 + t448 * t497 + t552 * t577);
t487 = -Ifges(3,1) + Ifges(3,2);
t503 = mrSges(3,1) * pkin(2);
t433 = t487 * t480 + 0.2e1 * (Ifges(3,4) * t491 + t503) * t497 + t491 * t534 + t530;
t514 = t421 * t558 - t433 * t439;
t385 = (t412 * t561 + t514 * t475) * t427;
t397 = (t412 * t442 + t526 * t475) * t427;
t523 = -t385 * t439 + t397 * t570;
t431 = t485 * t464 + t528 * t483;
t449 = pkin(3) * t541 + t484 * t461;
t452 = t483 * t540 - t485 * t500;
t489 = legFrame(2,2);
t473 = sin(t489);
t476 = cos(t489);
t413 = -(t452 * t476 - t473 * t550) * t576 + (t431 * t476 + t473 * t449) * t499 + (t486 * t473 + t476 * t555) * t579;
t428 = 0.1e1 / (pkin(2) * t541 + t449 * t499 + t550 * t576);
t434 = t487 * t481 + 0.2e1 * (Ifges(3,4) * t493 + t503) * t499 + t493 * t534 + t530;
t513 = t422 * t557 - t434 * t440;
t386 = (t413 * t560 + t513 * t476) * t428;
t398 = (t413 * t443 + t525 * t476) * t428;
t522 = -t386 * t440 + t398 * t569;
t432 = t485 * t465 + t527 * t483;
t450 = pkin(3) * t539 + t484 * t462;
t453 = t483 * t538 - t485 * t502;
t490 = legFrame(1,2);
t474 = sin(t490);
t477 = cos(t490);
t414 = -(t453 * t477 - t474 * t548) * t575 + (t432 * t477 + t474 * t450) * t501 + (t486 * t474 + t477 * t555) * t578;
t429 = 0.1e1 / (pkin(2) * t539 + t450 * t501 + t548 * t575);
t435 = t487 * t482 + 0.2e1 * (Ifges(3,4) * t495 + t503) * t501 + t495 * t534 + t530;
t512 = t423 * t556 - t435 * t441;
t387 = (t414 * t559 + t512 * t477) * t429;
t399 = (t414 * t444 + t524 * t477) * t429;
t521 = -t387 * t441 + t399 * t568;
t415 = (t451 * t472 + t475 * t552) * t577 + (-t430 * t472 + t475 * t448) * t497 + (-t472 * t555 + t475 * t486) * t580;
t388 = (t415 * t561 - t514 * t472) * t427;
t400 = (t415 * t442 - t526 * t472) * t427;
t520 = -t388 * t439 + t400 * t570;
t416 = (t452 * t473 + t476 * t550) * t576 + (-t431 * t473 + t476 * t449) * t499 + (-t473 * t555 + t476 * t486) * t579;
t389 = (t416 * t560 - t513 * t473) * t428;
t401 = (t416 * t443 - t525 * t473) * t428;
t519 = -t389 * t440 + t401 * t569;
t417 = (t453 * t474 + t477 * t548) * t575 + (-t432 * t474 + t477 * t450) * t501 + (-t474 * t555 + t477 * t486) * t578;
t390 = (t417 * t559 - t512 * t474) * t429;
t402 = (t417 * t444 - t524 * t474) * t429;
t518 = -t390 * t441 + t402 * t568;
t418 = -t454 * t577 - t463 * t483 * t497 + (pkin(2) * t553 + t529 * t497) * t485;
t436 = t491 * t451 + t483 * t547;
t403 = (t418 * t561 + t424 * t558 + t433 * t436) * t427;
t409 = (Ifges(3,3) * t567 + t418 * t442 + t436 * t457) * t427;
t517 = -t403 * t439 + t409 * t570;
t419 = -t455 * t576 - t464 * t483 * t499 + (pkin(2) * t551 + t528 * t499) * t485;
t437 = t493 * t452 + t483 * t546;
t404 = (t419 * t560 + t425 * t557 + t434 * t437) * t428;
t410 = (Ifges(3,3) * t566 + t419 * t443 + t437 * t458) * t428;
t516 = -t404 * t440 + t410 * t569;
t420 = -t456 * t575 - t465 * t483 * t501 + (pkin(2) * t549 + t527 * t501) * t485;
t438 = t495 * t453 + t483 * t545;
t405 = (t420 * t559 + t426 * t556 + t435 * t438) * t429;
t411 = (Ifges(3,3) * t565 + t420 * t444 + t438 * t459) * t429;
t515 = -t405 * t441 + t411 * t568;
t511 = t421 * t564 - t439 * t561;
t510 = t422 * t563 - t440 * t560;
t509 = t423 * t562 - t441 * t559;
t478 = m(1) + m(2) + m(3);
t408 = (t420 * t478 + t426 * t562 + t438 * t559) * t429;
t407 = (t419 * t478 + t425 * t563 + t437 * t560) * t428;
t406 = (t418 * t478 + t424 * t564 + t436 * t561) * t427;
t396 = (t417 * t478 - t509 * t474) * t429;
t395 = (t416 * t478 - t510 * t473) * t428;
t394 = (t415 * t478 - t511 * t472) * t427;
t393 = (t414 * t478 + t509 * t477) * t429;
t392 = (t413 * t478 + t510 * t476) * t428;
t391 = (t412 * t478 + t511 * t475) * t427;
t1 = [m(4) + (t393 * t414 + t521 * t477) * t429 + (t392 * t413 + t522 * t476) * t428 + (t391 * t412 + t523 * t475) * t427, (t393 * t417 - t521 * t474) * t429 + (t392 * t416 - t522 * t473) * t428 + (t391 * t415 - t523 * t472) * t427, (t387 * t438 + t393 * t420 + t399 * t565) * t429 + (t386 * t437 + t392 * t419 + t398 * t566) * t428 + (t385 * t436 + t391 * t418 + t397 * t567) * t427; (t396 * t414 + t518 * t477) * t429 + (t395 * t413 + t519 * t476) * t428 + (t394 * t412 + t520 * t475) * t427, m(4) + (t396 * t417 - t518 * t474) * t429 + (t395 * t416 - t519 * t473) * t428 + (t394 * t415 - t520 * t472) * t427, (t390 * t438 + t396 * t420 + t402 * t565) * t429 + (t389 * t437 + t395 * t419 + t401 * t566) * t428 + (t388 * t436 + t394 * t418 + t400 * t567) * t427; (t408 * t414 + t515 * t477) * t429 + (t407 * t413 + t516 * t476) * t428 + (t406 * t412 + t517 * t475) * t427, (t408 * t417 - t515 * t474) * t429 + (t407 * t416 - t516 * t473) * t428 + (t406 * t415 - t517 * t472) * t427, m(4) + (t405 * t438 + t408 * t420 + t411 * t565) * t429 + (t404 * t437 + t407 * t419 + t410 * t566) * t428 + (t403 * t436 + t406 * t418 + t409 * t567) * t427;];
MX  = t1;
