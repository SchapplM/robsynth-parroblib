% Calculate inertia matrix for parallel robot
% P3RRPRR12V1G2A0
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4]';
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
% Datum: 2020-08-06 19:07
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:05:44
% EndTime: 2020-08-06 19:05:45
% DurationCPUTime: 1.07s
% Computational Cost: add. (2739->224), mult. (3825->390), div. (432->6), fcn. (2727->18), ass. (0->177)
t447 = 0.1e1 / qJ(3,1);
t431 = sin(qJ(2,1));
t437 = cos(qJ(2,1));
t440 = pkin(1) + pkin(2);
t408 = t431 * qJ(3,1) + t440 * t437;
t432 = sin(qJ(1,1));
t438 = cos(qJ(1,1));
t450 = -t432 * pkin(4) + t408 * t438;
t495 = t450 * t447;
t445 = 0.1e1 / qJ(3,2);
t429 = sin(qJ(2,2));
t435 = cos(qJ(2,2));
t407 = t429 * qJ(3,2) + t440 * t435;
t430 = sin(qJ(1,2));
t436 = cos(qJ(1,2));
t451 = -t430 * pkin(4) + t407 * t436;
t496 = t451 * t445;
t443 = 0.1e1 / qJ(3,3);
t427 = sin(qJ(2,3));
t433 = cos(qJ(2,3));
t406 = t427 * qJ(3,3) + t440 * t433;
t428 = sin(qJ(1,3));
t434 = cos(qJ(1,3));
t452 = -t428 * pkin(4) + t406 * t434;
t497 = t452 * t443;
t442 = qJ(3,3) ^ 2;
t523 = m(3) * t442;
t444 = qJ(3,2) ^ 2;
t522 = m(3) * t444;
t446 = qJ(3,1) ^ 2;
t521 = m(3) * t446;
t520 = t434 * pkin(4);
t519 = t436 * pkin(4);
t518 = t438 * pkin(4);
t517 = Ifges(2,1) + Ifges(3,1);
t516 = Ifges(2,6) - Ifges(3,6);
t515 = mrSges(3,2) * t427;
t514 = mrSges(3,2) * t429;
t513 = mrSges(3,2) * t431;
t512 = mrSges(3,3) * qJ(3,1);
t511 = mrSges(3,3) * qJ(3,2);
t510 = mrSges(3,3) * qJ(3,3);
t473 = t427 * t428;
t397 = qJ(3,3) * t473 + t520;
t424 = legFrame(3,2);
t410 = sin(t424);
t413 = cos(t424);
t421 = t433 ^ 2;
t467 = t440 * t427;
t472 = t428 * t440;
t488 = t410 * qJ(3,3);
t367 = (t413 * t472 - t488) * t421 + (t397 * t413 + t410 * t467) * t433 + t488;
t509 = t367 * t443;
t471 = t429 * t430;
t398 = qJ(3,2) * t471 + t519;
t425 = legFrame(2,2);
t411 = sin(t425);
t414 = cos(t425);
t422 = t435 ^ 2;
t466 = t440 * t429;
t470 = t430 * t440;
t486 = t411 * qJ(3,2);
t368 = (t414 * t470 - t486) * t422 + (t398 * t414 + t411 * t466) * t435 + t486;
t508 = t368 * t445;
t469 = t431 * t432;
t399 = qJ(3,1) * t469 + t518;
t426 = legFrame(1,2);
t412 = sin(t426);
t415 = cos(t426);
t423 = t437 ^ 2;
t465 = t440 * t431;
t468 = t432 * t440;
t484 = t412 * qJ(3,1);
t369 = (t415 * t468 - t484) * t423 + (t399 * t415 + t412 * t465) * t437 + t484;
t507 = t369 * t447;
t482 = t413 * qJ(3,3);
t370 = (-t410 * t472 - t482) * t421 + (-t410 * t397 + t413 * t467) * t433 + t482;
t506 = t370 * t443;
t480 = t414 * qJ(3,2);
t371 = (-t411 * t470 - t480) * t422 + (-t411 * t398 + t414 * t466) * t435 + t480;
t505 = t371 * t445;
t478 = t415 * qJ(3,1);
t372 = (-t412 * t468 - t478) * t423 + (-t412 * t399 + t415 * t465) * t437 + t478;
t504 = t372 * t447;
t382 = t406 * t428 + t520;
t403 = -t433 * qJ(3,3) + t467;
t376 = t382 * t413 + t410 * t403;
t503 = t376 * t443;
t377 = -t382 * t410 + t403 * t413;
t502 = t377 * t443;
t383 = t407 * t430 + t519;
t404 = -t435 * qJ(3,2) + t466;
t378 = t383 * t414 + t411 * t404;
t501 = t378 * t445;
t379 = -t383 * t411 + t404 * t414;
t500 = t379 * t445;
t384 = t408 * t432 + t518;
t405 = -t437 * qJ(3,1) + t465;
t380 = t384 * t415 + t412 * t405;
t499 = t380 * t447;
t381 = -t384 * t412 + t405 * t415;
t498 = t381 * t447;
t409 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t391 = (mrSges(3,2) * qJ(3,3) + t516) * t433 - t427 * t409;
t494 = t391 * t443;
t392 = (mrSges(3,2) * qJ(3,2) + t516) * t435 - t429 * t409;
t493 = t392 * t445;
t393 = (qJ(3,1) * mrSges(3,2) + t516) * t437 - t431 * t409;
t492 = t393 * t447;
t417 = 0.2e1 * t510;
t448 = pkin(1) ^ 2;
t420 = 0.2e1 * mrSges(3,1) * pkin(1);
t462 = Ifges(3,2) + Ifges(2,3) + t420;
t394 = (t442 + t448) * m(3) + t417 + t462;
t491 = t394 * t443;
t418 = 0.2e1 * t511;
t395 = (t444 + t448) * m(3) + t418 + t462;
t490 = t395 * t445;
t419 = 0.2e1 * t512;
t396 = (t446 + t448) * m(3) + t419 + t462;
t489 = t396 * t447;
t487 = t410 * t434;
t485 = t411 * t436;
t483 = t412 * t438;
t481 = t413 * t434;
t479 = t414 * t436;
t477 = t415 * t438;
t416 = m(3) * pkin(1) + mrSges(3,1);
t476 = t416 * t443;
t475 = t416 * t445;
t474 = t416 * t447;
t464 = Ifges(1,3) + t517;
t463 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t461 = t434 * t515;
t460 = t443 * t515;
t459 = t436 * t514;
t458 = t445 * t514;
t457 = t438 * t513;
t456 = t447 * t513;
t455 = t433 * t497;
t454 = t435 * t496;
t453 = t437 * t495;
t449 = m(3) * t448 + Ifges(2,2) + Ifges(3,3) + t420 - t517;
t402 = 0.1e1 / t408;
t401 = 0.1e1 / t407;
t400 = 0.1e1 / t406;
t375 = (t449 - 0.2e1 * t512 - t521) * t423 + 0.2e1 * (t416 * qJ(3,1) + t463) * t431 * t437 + t521 + t419 + t464;
t374 = (t449 - 0.2e1 * t511 - t522) * t422 + 0.2e1 * (t416 * qJ(3,2) + t463) * t429 * t435 + t522 + t418 + t464;
t373 = (t449 - 0.2e1 * t510 - t523) * t421 + 0.2e1 * (t416 * qJ(3,3) + t463) * t427 * t433 + t523 + t417 + t464;
t366 = m(3) * t495 + (-mrSges(3,2) * t469 - t416 * t453) * t402;
t365 = m(3) * t496 + (-mrSges(3,2) * t471 - t416 * t454) * t401;
t364 = m(3) * t497 + (-mrSges(3,2) * t473 - t416 * t455) * t400;
t363 = -t450 * t474 + (-t393 * t432 + t396 * t453) * t402;
t362 = -t451 * t475 + (-t392 * t430 + t395 * t454) * t401;
t361 = -t452 * t476 + (-t391 * t428 + t394 * t455) * t400;
t360 = t450 * t456 + (-t375 * t432 + t393 * t453) * t402;
t359 = t451 * t458 + (-t374 * t430 + t392 * t454) * t401;
t358 = t452 * t460 + (-t373 * t428 + t391 * t455) * t400;
t357 = m(3) * t499 + (-t369 * t474 + t415 * t457) * t402;
t356 = m(3) * t501 + (-t368 * t475 + t414 * t459) * t401;
t355 = m(3) * t503 + (-t367 * t476 + t413 * t461) * t400;
t354 = m(3) * t498 + (-t372 * t474 - t412 * t457) * t402;
t353 = m(3) * t500 + (-t371 * t475 - t411 * t459) * t401;
t352 = m(3) * t502 + (-t370 * t476 - t410 * t461) * t400;
t351 = -t380 * t474 + (t369 * t489 + t393 * t477) * t402;
t350 = -t378 * t475 + (t368 * t490 + t392 * t479) * t401;
t349 = -t376 * t476 + (t367 * t491 + t391 * t481) * t400;
t348 = -t381 * t474 + (t372 * t489 - t393 * t483) * t402;
t347 = -t379 * t475 + (t371 * t490 - t392 * t485) * t401;
t346 = -t377 * t476 + (t370 * t491 - t391 * t487) * t400;
t345 = t380 * t456 + (t369 * t492 + t375 * t477) * t402;
t344 = t378 * t458 + (t368 * t493 + t374 * t479) * t401;
t343 = t376 * t460 + (t367 * t494 + t373 * t481) * t400;
t342 = t381 * t456 + (t372 * t492 - t375 * t483) * t402;
t341 = t379 * t458 + (t371 * t493 - t374 * t485) * t401;
t340 = t377 * t460 + (t370 * t494 - t373 * t487) * t400;
t1 = [t355 * t503 + t356 * t501 + t357 * t499 + m(4) + (t345 * t477 + t351 * t507) * t402 + (t344 * t479 + t350 * t508) * t401 + (t343 * t481 + t349 * t509) * t400, t355 * t502 + t356 * t500 + t357 * t498 + (-t345 * t483 + t351 * t504) * t402 + (-t344 * t485 + t350 * t505) * t401 + (-t343 * t487 + t349 * t506) * t400, t355 * t497 + t356 * t496 + t357 * t495 + (-t345 * t432 + t351 * t453) * t402 + (-t344 * t430 + t350 * t454) * t401 + (-t343 * t428 + t349 * t455) * t400; t352 * t503 + t353 * t501 + t354 * t499 + (t342 * t477 + t348 * t507) * t402 + (t341 * t479 + t347 * t508) * t401 + (t340 * t481 + t346 * t509) * t400, t352 * t502 + t353 * t500 + t354 * t498 + m(4) + (-t342 * t483 + t348 * t504) * t402 + (-t341 * t485 + t347 * t505) * t401 + (-t340 * t487 + t346 * t506) * t400, t352 * t497 + t353 * t496 + t354 * t495 + (-t342 * t432 + t348 * t453) * t402 + (-t341 * t430 + t347 * t454) * t401 + (-t340 * t428 + t346 * t455) * t400; t364 * t503 + t365 * t501 + t366 * t499 + (t360 * t477 + t363 * t507) * t402 + (t359 * t479 + t362 * t508) * t401 + (t358 * t481 + t361 * t509) * t400, t364 * t502 + t365 * t500 + t366 * t498 + (-t360 * t483 + t363 * t504) * t402 + (-t359 * t485 + t362 * t505) * t401 + (-t358 * t487 + t361 * t506) * t400, t364 * t497 + t365 * t496 + t366 * t495 + m(4) + (-t360 * t432 + t363 * t453) * t402 + (-t359 * t430 + t362 * t454) * t401 + (-t358 * t428 + t361 * t455) * t400;];
MX  = t1;
