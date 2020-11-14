% Calculate inertia matrix for parallel robot
% P3RRPRR12V1G3A0
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
% Datum: 2020-08-06 19:11
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR12V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 19:10:03
% EndTime: 2020-08-06 19:10:04
% DurationCPUTime: 1.63s
% Computational Cost: add. (2739->224), mult. (3753->390), div. (432->6), fcn. (2655->18), ass. (0->177)
t449 = 0.1e1 / qJ(3,1);
t433 = sin(qJ(2,1));
t439 = cos(qJ(2,1));
t442 = pkin(1) + pkin(2);
t407 = t433 * qJ(3,1) + t442 * t439;
t434 = sin(qJ(1,1));
t440 = cos(qJ(1,1));
t527 = t440 * pkin(4) + t407 * t434;
t494 = t527 * t449;
t447 = 0.1e1 / qJ(3,2);
t431 = sin(qJ(2,2));
t437 = cos(qJ(2,2));
t406 = t431 * qJ(3,2) + t442 * t437;
t432 = sin(qJ(1,2));
t438 = cos(qJ(1,2));
t529 = t438 * pkin(4) + t406 * t432;
t495 = t529 * t447;
t445 = 0.1e1 / qJ(3,3);
t429 = sin(qJ(2,3));
t435 = cos(qJ(2,3));
t405 = t429 * qJ(3,3) + t442 * t435;
t430 = sin(qJ(1,3));
t436 = cos(qJ(1,3));
t531 = t436 * pkin(4) + t405 * t430;
t496 = t531 * t445;
t522 = t430 * pkin(4);
t530 = t405 * t436 - t522;
t521 = t432 * pkin(4);
t528 = t406 * t438 - t521;
t520 = t434 * pkin(4);
t526 = t407 * t440 - t520;
t444 = qJ(3,3) ^ 2;
t525 = m(3) * t444;
t446 = qJ(3,2) ^ 2;
t524 = m(3) * t446;
t448 = qJ(3,1) ^ 2;
t523 = m(3) * t448;
t516 = Ifges(2,1) + Ifges(3,1);
t515 = Ifges(2,6) - Ifges(3,6);
t514 = mrSges(3,2) * t429;
t513 = mrSges(3,2) * t431;
t512 = mrSges(3,2) * t433;
t511 = mrSges(3,3) * qJ(3,1);
t510 = mrSges(3,3) * qJ(3,2);
t509 = mrSges(3,3) * qJ(3,3);
t472 = t429 * t436;
t396 = qJ(3,3) * t472 - t522;
t426 = legFrame(3,2);
t412 = sin(t426);
t415 = cos(t426);
t423 = t435 ^ 2;
t466 = t442 * t429;
t469 = t436 * t442;
t487 = t412 * qJ(3,3);
t369 = (t415 * t469 - t487) * t423 + (t396 * t415 + t412 * t466) * t435 + t487;
t508 = t369 * t445;
t471 = t431 * t438;
t397 = qJ(3,2) * t471 - t521;
t427 = legFrame(2,2);
t413 = sin(t427);
t416 = cos(t427);
t424 = t437 ^ 2;
t465 = t442 * t431;
t468 = t438 * t442;
t485 = t413 * qJ(3,2);
t370 = (t416 * t468 - t485) * t424 + (t397 * t416 + t413 * t465) * t437 + t485;
t507 = t370 * t447;
t470 = t433 * t440;
t398 = qJ(3,1) * t470 - t520;
t428 = legFrame(1,2);
t414 = sin(t428);
t417 = cos(t428);
t425 = t439 ^ 2;
t464 = t442 * t433;
t467 = t440 * t442;
t483 = t414 * qJ(3,1);
t371 = (t417 * t467 - t483) * t425 + (t398 * t417 + t414 * t464) * t439 + t483;
t506 = t371 * t449;
t481 = t415 * qJ(3,3);
t372 = (-t412 * t469 - t481) * t423 + (-t412 * t396 + t415 * t466) * t435 + t481;
t505 = t372 * t445;
t479 = t416 * qJ(3,2);
t373 = (-t413 * t468 - t479) * t424 + (-t413 * t397 + t416 * t465) * t437 + t479;
t504 = t373 * t447;
t477 = t417 * qJ(3,1);
t374 = (-t414 * t467 - t477) * t425 + (-t414 * t398 + t417 * t464) * t439 + t477;
t503 = t374 * t449;
t402 = -t435 * qJ(3,3) + t466;
t378 = t402 * t412 + t530 * t415;
t502 = t378 * t445;
t403 = -t437 * qJ(3,2) + t465;
t379 = t403 * t413 + t528 * t416;
t501 = t379 * t447;
t404 = -t439 * qJ(3,1) + t464;
t380 = t404 * t414 + t526 * t417;
t500 = t380 * t449;
t381 = t402 * t415 - t530 * t412;
t499 = t381 * t445;
t382 = t403 * t416 - t528 * t413;
t498 = t382 * t447;
t383 = t404 * t417 - t526 * t414;
t497 = t383 * t449;
t411 = pkin(1) * mrSges(3,2) - Ifges(3,4) - Ifges(2,5);
t390 = (mrSges(3,2) * qJ(3,3) + t515) * t435 - t429 * t411;
t493 = t390 * t445;
t391 = (mrSges(3,2) * qJ(3,2) + t515) * t437 - t431 * t411;
t492 = t391 * t447;
t392 = (qJ(3,1) * mrSges(3,2) + t515) * t439 - t433 * t411;
t491 = t392 * t449;
t419 = 0.2e1 * t509;
t450 = pkin(1) ^ 2;
t422 = 0.2e1 * mrSges(3,1) * pkin(1);
t461 = Ifges(3,2) + Ifges(2,3) + t422;
t393 = (t444 + t450) * m(3) + t419 + t461;
t490 = t393 * t445;
t420 = 0.2e1 * t510;
t394 = (t446 + t450) * m(3) + t420 + t461;
t489 = t394 * t447;
t421 = 0.2e1 * t511;
t395 = (t448 + t450) * m(3) + t421 + t461;
t488 = t395 * t449;
t486 = t412 * t430;
t484 = t413 * t432;
t482 = t414 * t434;
t480 = t415 * t430;
t478 = t416 * t432;
t476 = t417 * t434;
t418 = m(3) * pkin(1) + mrSges(3,1);
t475 = t418 * t445;
t474 = t418 * t447;
t473 = t418 * t449;
t463 = Ifges(1,3) + t516;
t462 = pkin(1) * mrSges(3,3) + Ifges(2,4) - Ifges(3,5);
t460 = t430 * t514;
t459 = t445 * t514;
t458 = t432 * t513;
t457 = t447 * t513;
t456 = t434 * t512;
t455 = t449 * t512;
t454 = t435 * t496;
t453 = t437 * t495;
t452 = t439 * t494;
t451 = m(3) * t450 + Ifges(2,2) + Ifges(3,3) + t422 - t516;
t401 = 0.1e1 / t407;
t400 = 0.1e1 / t406;
t399 = 0.1e1 / t405;
t377 = (t451 - 0.2e1 * t511 - t523) * t425 + 0.2e1 * (t418 * qJ(3,1) + t462) * t433 * t439 + t523 + t421 + t463;
t376 = (t451 - 0.2e1 * t510 - t524) * t424 + 0.2e1 * (t418 * qJ(3,2) + t462) * t431 * t437 + t524 + t420 + t463;
t375 = (t451 - 0.2e1 * t509 - t525) * t423 + 0.2e1 * (t418 * qJ(3,3) + t462) * t429 * t435 + t525 + t419 + t463;
t368 = -m(3) * t494 + (-mrSges(3,2) * t470 + t418 * t452) * t401;
t367 = -m(3) * t495 + (-mrSges(3,2) * t471 + t418 * t453) * t400;
t366 = -m(3) * t496 + (-mrSges(3,2) * t472 + t418 * t454) * t399;
t365 = t527 * t473 + (-t392 * t440 - t395 * t452) * t401;
t364 = t529 * t474 + (-t391 * t438 - t394 * t453) * t400;
t363 = t531 * t475 + (-t390 * t436 - t393 * t454) * t399;
t362 = -t527 * t455 + (-t377 * t440 - t392 * t452) * t401;
t361 = -t529 * t457 + (-t376 * t438 - t391 * t453) * t400;
t360 = -t531 * t459 + (-t375 * t436 - t390 * t454) * t399;
t359 = m(3) * t500 + (-t371 * t473 - t417 * t456) * t401;
t358 = m(3) * t501 + (-t370 * t474 - t416 * t458) * t400;
t357 = m(3) * t502 + (-t369 * t475 - t415 * t460) * t399;
t356 = m(3) * t497 + (-t374 * t473 + t414 * t456) * t401;
t355 = m(3) * t498 + (-t373 * t474 + t413 * t458) * t400;
t354 = m(3) * t499 + (-t372 * t475 + t412 * t460) * t399;
t353 = -t380 * t473 + (t371 * t488 - t392 * t476) * t401;
t352 = -t379 * t474 + (t370 * t489 - t391 * t478) * t400;
t351 = -t378 * t475 + (t369 * t490 - t390 * t480) * t399;
t350 = -t383 * t473 + (t374 * t488 + t392 * t482) * t401;
t349 = -t382 * t474 + (t373 * t489 + t391 * t484) * t400;
t348 = -t381 * t475 + (t372 * t490 + t390 * t486) * t399;
t347 = t380 * t455 + (t371 * t491 - t377 * t476) * t401;
t346 = t379 * t457 + (t370 * t492 - t376 * t478) * t400;
t345 = t378 * t459 + (t369 * t493 - t375 * t480) * t399;
t344 = t383 * t455 + (t374 * t491 + t377 * t482) * t401;
t343 = t382 * t457 + (t373 * t492 + t376 * t484) * t400;
t342 = t381 * t459 + (t372 * t493 + t375 * t486) * t399;
t1 = [t357 * t502 + t358 * t501 + t359 * t500 + m(4) + (-t347 * t476 + t353 * t506) * t401 + (-t346 * t478 + t352 * t507) * t400 + (-t345 * t480 + t351 * t508) * t399, t357 * t499 + t358 * t498 + t359 * t497 + (t347 * t482 + t353 * t503) * t401 + (t346 * t484 + t352 * t504) * t400 + (t345 * t486 + t351 * t505) * t399, -t357 * t496 - t358 * t495 - t359 * t494 + (-t347 * t440 - t353 * t452) * t401 + (-t346 * t438 - t352 * t453) * t400 + (-t345 * t436 - t351 * t454) * t399; t354 * t502 + t355 * t501 + t356 * t500 + (-t344 * t476 + t350 * t506) * t401 + (-t343 * t478 + t349 * t507) * t400 + (-t342 * t480 + t348 * t508) * t399, t354 * t499 + t355 * t498 + t356 * t497 + m(4) + (t344 * t482 + t350 * t503) * t401 + (t343 * t484 + t349 * t504) * t400 + (t342 * t486 + t348 * t505) * t399, -t354 * t496 - t355 * t495 - t356 * t494 + (-t344 * t440 - t350 * t452) * t401 + (-t343 * t438 - t349 * t453) * t400 + (-t342 * t436 - t348 * t454) * t399; t366 * t502 + t367 * t501 + t368 * t500 + (-t362 * t476 + t365 * t506) * t401 + (-t361 * t478 + t364 * t507) * t400 + (-t360 * t480 + t363 * t508) * t399, t366 * t499 + t367 * t498 + t368 * t497 + (t362 * t482 + t365 * t503) * t401 + (t361 * t484 + t364 * t504) * t400 + (t360 * t486 + t363 * t505) * t399, -t366 * t496 - t367 * t495 - t368 * t494 + m(4) + (-t362 * t440 - t365 * t452) * t401 + (-t361 * t438 - t364 * t453) * t400 + (-t360 * t436 - t363 * t454) * t399;];
MX  = t1;
