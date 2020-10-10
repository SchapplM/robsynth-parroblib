% Calculate inertia matrix for parallel robot
% P3PRRRR8V1G3A0
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
%   pkin=[a2,a4,alpha2,d2,d4,theta1]';
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
% Datum: 2020-08-06 17:17
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR8V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 17:16:10
% EndTime: 2020-08-06 17:16:11
% DurationCPUTime: 1.14s
% Computational Cost: add. (2667->199), mult. (7119->392), div. (648->7), fcn. (7497->22), ass. (0->192)
t561 = 2 * Ifges(3,4);
t457 = sin(qJ(2,3));
t463 = cos(qJ(2,3));
t462 = cos(qJ(3,3));
t554 = pkin(2) * t462;
t422 = -pkin(5) * t463 + t457 * t554;
t448 = sin(pkin(3));
t450 = cos(pkin(3));
t456 = sin(qJ(3,3));
t557 = pkin(2) * t456;
t413 = t422 * t448 + t450 * t557;
t560 = 0.1e1 / t413;
t459 = sin(qJ(2,2));
t465 = cos(qJ(2,2));
t464 = cos(qJ(3,2));
t553 = pkin(2) * t464;
t423 = -pkin(5) * t465 + t459 * t553;
t458 = sin(qJ(3,2));
t556 = pkin(2) * t458;
t414 = t423 * t448 + t450 * t556;
t559 = 0.1e1 / t414;
t461 = sin(qJ(2,1));
t467 = cos(qJ(2,1));
t466 = cos(qJ(3,1));
t552 = pkin(2) * t466;
t424 = -pkin(5) * t467 + t461 * t552;
t460 = sin(qJ(3,1));
t555 = pkin(2) * t460;
t415 = t424 * t448 + t450 * t555;
t558 = 0.1e1 / t415;
t551 = Ifges(3,1) + Ifges(2,3);
t468 = 0.1e1 / pkin(2);
t550 = Ifges(3,3) * t468;
t425 = pkin(5) * t457 + t463 * t554;
t447 = sin(pkin(6));
t449 = cos(pkin(6));
t483 = -t422 * t450 + t448 * t557;
t386 = t425 * t449 + t483 * t447;
t453 = legFrame(3,2);
t437 = sin(t453);
t440 = cos(t453);
t380 = t386 * t440 + t413 * t437;
t549 = t380 * t560;
t381 = -t386 * t437 + t413 * t440;
t548 = t381 * t560;
t426 = pkin(5) * t459 + t465 * t553;
t482 = -t423 * t450 + t448 * t556;
t387 = t426 * t449 + t482 * t447;
t454 = legFrame(2,2);
t438 = sin(t454);
t441 = cos(t454);
t382 = t387 * t441 + t414 * t438;
t547 = t382 * t559;
t383 = -t387 * t438 + t414 * t441;
t546 = t383 * t559;
t427 = pkin(5) * t461 + t467 * t552;
t481 = -t424 * t450 + t448 * t555;
t388 = t427 * t449 + t481 * t447;
t455 = legFrame(1,2);
t439 = sin(t455);
t442 = cos(t455);
t384 = t388 * t442 + t415 * t439;
t545 = t384 * t558;
t385 = -t388 * t439 + t415 * t442;
t544 = t385 * t558;
t389 = -t425 * t447 + t483 * t449;
t543 = t389 * t560;
t390 = -t426 * t447 + t482 * t449;
t542 = t390 * t559;
t391 = -t427 * t447 + t481 * t449;
t541 = t391 * t558;
t498 = t462 * mrSges(3,1) - mrSges(3,2) * t456;
t404 = t498 * t450 - t448 * (mrSges(3,1) * t456 + mrSges(3,2) * t462) * t457;
t540 = t404 * t560;
t539 = t404 * t468;
t497 = t464 * mrSges(3,1) - mrSges(3,2) * t458;
t405 = t497 * t450 - t448 * (mrSges(3,1) * t458 + mrSges(3,2) * t464) * t459;
t538 = t405 * t559;
t537 = t405 * t468;
t496 = t466 * mrSges(3,1) - mrSges(3,2) * t460;
t406 = t496 * t450 - t448 * (mrSges(3,1) * t460 + mrSges(3,2) * t466) * t461;
t536 = t406 * t558;
t535 = t406 * t468;
t534 = t560 / t462;
t533 = t559 / t464;
t532 = t558 / t466;
t443 = m(1) + m(2) + m(3);
t531 = t560 * t443;
t530 = t559 * t443;
t529 = t558 * t443;
t452 = mrSges(2,2) - mrSges(3,3);
t528 = ((mrSges(2,1) + t498) * t463 - t457 * t452) * t448;
t527 = ((mrSges(2,1) + t497) * t465 - t459 * t452) * t448;
t526 = ((mrSges(2,1) + t496) * t467 - t461 * t452) * t448;
t428 = Ifges(3,5) * t456 + Ifges(3,6) * t462;
t525 = t428 * t468;
t429 = Ifges(3,5) * t458 + Ifges(3,6) * t464;
t524 = t429 * t468;
t430 = Ifges(3,5) * t460 + Ifges(3,6) * t466;
t523 = t430 * t468;
t522 = t450 * t457;
t521 = t450 * t459;
t520 = t450 * t461;
t519 = t450 * t463;
t518 = t450 * t465;
t517 = t450 * t467;
t516 = t456 * t463;
t515 = t458 * t465;
t514 = t460 * t467;
t395 = (t447 * t519 + t449 * t457) * t554 + (t447 * t522 - t449 * t463) * pkin(5);
t513 = t395 * t534;
t396 = (t447 * t518 + t449 * t459) * t553 + (t447 * t521 - t449 * t465) * pkin(5);
t512 = t396 * t533;
t397 = (t447 * t517 + t449 * t461) * t552 + (t447 * t520 - t449 * t467) * pkin(5);
t511 = t397 * t532;
t474 = t448 * t462 + t456 * t522;
t398 = t474 * t447 - t449 * t516;
t510 = t398 * t534;
t473 = t448 * t464 + t458 * t521;
t399 = t473 * t447 - t449 * t515;
t509 = t399 * t533;
t472 = t448 * t466 + t460 * t520;
t400 = t472 * t447 - t449 * t514;
t508 = t400 * t532;
t507 = t437 * t534;
t506 = t440 * t534;
t505 = t438 * t533;
t504 = t441 * t533;
t503 = t439 * t532;
t502 = t442 * t532;
t501 = t560 * t528;
t500 = t559 * t527;
t499 = t558 * t526;
t392 = (-t447 * t457 + t449 * t519) * t554 + pkin(5) * (t447 * t463 + t449 * t522);
t495 = t392 * t507;
t494 = t392 * t506;
t393 = (-t447 * t459 + t449 * t518) * t553 + pkin(5) * (t447 * t465 + t449 * t521);
t493 = t393 * t505;
t492 = t393 * t504;
t394 = (-t447 * t461 + t449 * t517) * t552 + pkin(5) * (t447 * t467 + t449 * t520);
t491 = t394 * t503;
t490 = t394 * t502;
t401 = t447 * t516 + t474 * t449;
t489 = t401 * t507;
t488 = t401 * t506;
t402 = t447 * t515 + t473 * t449;
t487 = t402 * t505;
t486 = t402 * t504;
t403 = t447 * t514 + t472 * t449;
t485 = t403 * t503;
t484 = t403 * t502;
t480 = t392 * t550 + t401 * t428;
t479 = t393 * t550 + t402 * t429;
t478 = t394 * t550 + t403 * t430;
t451 = -Ifges(3,1) + Ifges(3,2);
t419 = (t451 * t462 + t456 * t561) * t462 + t551;
t477 = t392 * t525 + t401 * t419;
t420 = (t451 * t464 + t458 * t561) * t464 + t551;
t476 = t393 * t524 + t402 * t420;
t421 = (t451 * t466 + t460 * t561) * t466 + t551;
t475 = t394 * t523 + t403 * t421;
t471 = t392 * t539 + t401 * t528;
t470 = t393 * t537 + t402 * t527;
t469 = t394 * t535 + t403 * t526;
t379 = t391 * t536 + (t397 * t550 + t400 * t430) * t532;
t378 = t390 * t538 + (t396 * t550 + t399 * t429) * t533;
t377 = t389 * t540 + (t395 * t550 + t398 * t428) * t534;
t376 = t391 * t499 + (t397 * t523 + t400 * t421) * t532;
t375 = t390 * t500 + (t396 * t524 + t399 * t420) * t533;
t374 = t389 * t501 + (t395 * t525 + t398 * t419) * t534;
t373 = t391 * t529 + (t397 * t535 + t400 * t526) * t532;
t372 = t390 * t530 + (t396 * t537 + t399 * t527) * t533;
t371 = t389 * t531 + (t395 * t539 + t398 * t528) * t534;
t370 = t385 * t536 + t478 * t503;
t369 = t384 * t536 - t478 * t502;
t368 = t383 * t538 + t479 * t505;
t367 = t382 * t538 - t479 * t504;
t366 = t381 * t540 + t480 * t507;
t365 = t380 * t540 - t480 * t506;
t364 = t385 * t499 + t475 * t503;
t363 = t384 * t499 - t475 * t502;
t362 = t383 * t500 + t476 * t505;
t361 = t382 * t500 - t476 * t504;
t360 = t381 * t501 + t477 * t507;
t359 = t380 * t501 - t477 * t506;
t358 = t385 * t529 + t469 * t503;
t357 = t384 * t529 - t469 * t502;
t356 = t383 * t530 + t470 * t505;
t355 = t382 * t530 - t470 * t504;
t354 = t381 * t531 + t471 * t507;
t353 = t380 * t531 - t471 * t506;
t1 = [-t359 * t488 - t361 * t486 - t363 * t484 + t353 * t549 + t355 * t547 + t357 * t545 + m(4) + (-t365 * t494 - t367 * t492 - t369 * t490) * t468, t359 * t489 + t361 * t487 + t363 * t485 + t353 * t548 + t355 * t546 + t357 * t544 + (t365 * t495 + t367 * t493 + t369 * t491) * t468, t359 * t510 + t361 * t509 + t363 * t508 + t353 * t543 + t355 * t542 + t357 * t541 + (t365 * t513 + t367 * t512 + t369 * t511) * t468; -t360 * t488 - t362 * t486 - t364 * t484 + t354 * t549 + t356 * t547 + t358 * t545 + (-t366 * t494 - t368 * t492 - t370 * t490) * t468, t360 * t489 + t362 * t487 + t364 * t485 + t354 * t548 + t356 * t546 + t358 * t544 + m(4) + (t366 * t495 + t368 * t493 + t370 * t491) * t468, t360 * t510 + t362 * t509 + t364 * t508 + t354 * t543 + t356 * t542 + t358 * t541 + (t366 * t513 + t368 * t512 + t370 * t511) * t468; -t374 * t488 - t375 * t486 - t376 * t484 + t371 * t549 + t372 * t547 + t373 * t545 + (-t377 * t494 - t378 * t492 - t379 * t490) * t468, t374 * t489 + t375 * t487 + t376 * t485 + t371 * t548 + t372 * t546 + t373 * t544 + (t377 * t495 + t378 * t493 + t379 * t491) * t468, t374 * t510 + t375 * t509 + t376 * t508 + t371 * t543 + t372 * t542 + t373 * t541 + m(4) + (t377 * t513 + t378 * t512 + t379 * t511) * t468;];
MX  = t1;
