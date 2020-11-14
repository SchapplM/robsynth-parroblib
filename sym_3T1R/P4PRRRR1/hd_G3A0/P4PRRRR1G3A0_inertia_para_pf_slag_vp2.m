% Calculate inertia matrix for parallel robot
% P4PRRRR1G3A0
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
% Datum: 2020-03-02 19:06
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRRR1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRRR1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 19:02:57
% EndTime: 2020-03-02 19:02:59
% DurationCPUTime: 2.44s
% Computational Cost: add. (1800->254), mult. (2548->464), div. (1184->13), fcn. (3344->26), ass. (0->202)
t533 = 2 * Ifges(3,4);
t426 = cos(qJ(3,4));
t412 = 0.1e1 / t426;
t532 = t412 ^ 2;
t440 = cos(qJ(3,3));
t418 = 0.1e1 / t440;
t531 = t418 ^ 2;
t442 = cos(qJ(3,2));
t420 = 0.1e1 / t442;
t530 = t420 ^ 2;
t444 = cos(qJ(3,1));
t422 = 0.1e1 / t444;
t529 = t422 ^ 2;
t528 = Ifges(3,1) + Ifges(2,3);
t430 = legFrame(4,2);
t401 = sin(t430);
t405 = cos(t430);
t425 = sin(qJ(2,4));
t427 = cos(qJ(2,4));
t373 = -t401 * t427 + t405 * t425;
t411 = 0.1e1 / t425;
t527 = t373 * t411;
t374 = t425 * t401 + t405 * t427;
t526 = t374 * t411;
t431 = legFrame(3,2);
t402 = sin(t431);
t406 = cos(t431);
t435 = sin(qJ(2,3));
t441 = cos(qJ(2,3));
t378 = -t402 * t441 + t406 * t435;
t415 = 0.1e1 / t435;
t525 = t378 * t415;
t379 = t435 * t402 + t406 * t441;
t524 = t379 * t415;
t432 = legFrame(2,2);
t403 = sin(t432);
t407 = cos(t432);
t437 = sin(qJ(2,2));
t443 = cos(qJ(2,2));
t380 = -t403 * t443 + t407 * t437;
t416 = 0.1e1 / t437;
t523 = t380 * t416;
t381 = t437 * t403 + t407 * t443;
t522 = t381 * t416;
t433 = legFrame(1,2);
t404 = sin(t433);
t408 = cos(t433);
t439 = sin(qJ(2,1));
t445 = cos(qJ(2,1));
t382 = -t404 * t445 + t408 * t439;
t417 = 0.1e1 / t439;
t521 = t382 * t417;
t383 = t439 * t404 + t408 * t445;
t520 = t383 * t417;
t424 = sin(qJ(3,4));
t393 = mrSges(3,1) * t424 + mrSges(3,2) * t426;
t519 = t393 * t412;
t518 = t393 * t425;
t434 = sin(qJ(3,3));
t397 = mrSges(3,1) * t434 + mrSges(3,2) * t440;
t517 = t397 * t418;
t516 = t397 * t435;
t436 = sin(qJ(3,2));
t398 = mrSges(3,1) * t436 + mrSges(3,2) * t442;
t515 = t398 * t420;
t514 = t398 * t437;
t438 = sin(qJ(3,1));
t399 = mrSges(3,1) * t438 + mrSges(3,2) * t444;
t513 = t399 * t422;
t512 = t399 * t439;
t511 = t411 * t412;
t510 = t411 * t424;
t457 = 0.1e1 / pkin(2);
t509 = t412 * t457;
t458 = t426 ^ 2;
t508 = 0.1e1 / t458 * t427;
t507 = t415 * t418;
t506 = t415 * t434;
t505 = t416 * t420;
t504 = t416 * t436;
t503 = t417 * t422;
t502 = t417 * t438;
t501 = t418 * t457;
t459 = t440 ^ 2;
t500 = 0.1e1 / t459 * t441;
t499 = t420 * t457;
t460 = t442 ^ 2;
t498 = 0.1e1 / t460 * t443;
t497 = t422 * t457;
t461 = t444 ^ 2;
t496 = 0.1e1 / t461 * t445;
t392 = Ifges(3,5) * t424 + Ifges(3,6) * t426;
t495 = t392 * t411 * t532;
t394 = Ifges(3,5) * t434 + Ifges(3,6) * t440;
t494 = t394 * t415 * t531;
t395 = Ifges(3,5) * t436 + Ifges(3,6) * t442;
t493 = t395 * t416 * t530;
t396 = Ifges(3,5) * t438 + Ifges(3,6) * t444;
t492 = t396 * t417 * t529;
t491 = t401 * t511;
t490 = t401 * t509;
t489 = t402 * t507;
t488 = t402 * t501;
t487 = t403 * t505;
t486 = t403 * t499;
t485 = t404 * t503;
t484 = t404 * t497;
t483 = t405 * t511;
t482 = t405 * t509;
t481 = t406 * t507;
t480 = t406 * t501;
t479 = t407 * t505;
t478 = t407 * t499;
t477 = t408 * t503;
t476 = t408 * t497;
t475 = t412 * t510;
t474 = t457 * t508;
t473 = t418 * t506;
t472 = t420 * t504;
t471 = t422 * t502;
t470 = t457 * t500;
t469 = t457 * t498;
t468 = t457 * t496;
t446 = xP(4);
t409 = sin(t446);
t410 = cos(t446);
t447 = mrSges(4,2);
t448 = mrSges(4,1);
t467 = -t409 * t447 + t410 * t448;
t466 = t508 * t510;
t465 = t500 * t506;
t464 = t498 * t504;
t463 = t496 * t502;
t462 = -t409 * t448 - t410 * t447;
t456 = koppelP(1,1);
t455 = koppelP(2,1);
t454 = koppelP(3,1);
t453 = koppelP(4,1);
t452 = koppelP(1,2);
t451 = koppelP(2,2);
t450 = koppelP(3,2);
t449 = koppelP(4,2);
t429 = mrSges(2,2) - mrSges(3,3);
t428 = -Ifges(3,1) + Ifges(3,2);
t414 = m(1) + m(2) + m(3);
t391 = -t409 * t452 + t410 * t456;
t390 = -t409 * t451 + t410 * t455;
t389 = -t409 * t450 + t410 * t454;
t388 = -t409 * t449 + t410 * t453;
t387 = -t409 * t456 - t410 * t452;
t386 = -t409 * t455 - t410 * t451;
t385 = -t409 * t454 - t410 * t450;
t384 = -t409 * t453 - t410 * t449;
t377 = t438 * t444 * t533 + t428 * t461 + t528;
t376 = t436 * t442 * t533 + t428 * t460 + t528;
t375 = t434 * t440 * t533 + t428 * t459 + t528;
t372 = t424 * t426 * t533 + t428 * t458 + t528;
t371 = (t444 * mrSges(3,1) - t438 * mrSges(3,2) + mrSges(2,1)) * t445 - t439 * t429;
t370 = (t442 * mrSges(3,1) - t436 * mrSges(3,2) + mrSges(2,1)) * t443 - t437 * t429;
t369 = (t440 * mrSges(3,1) - t434 * mrSges(3,2) + mrSges(2,1)) * t441 - t435 * t429;
t368 = (t426 * mrSges(3,1) - t424 * mrSges(3,2) + mrSges(2,1)) * t427 - t425 * t429;
t367 = (-t387 * t408 + t391 * t404) * t417 * t497;
t366 = (-t386 * t407 + t390 * t403) * t416 * t499;
t365 = (-t385 * t406 + t389 * t402) * t415 * t501;
t364 = (-t384 * t405 + t388 * t401) * t411 * t509;
t363 = (-t371 * t476 + t383 * t414) * t417;
t362 = (t371 * t484 + t382 * t414) * t417;
t361 = (-t370 * t478 + t381 * t414) * t416;
t360 = (t370 * t486 + t380 * t414) * t416;
t359 = (-t369 * t480 + t379 * t414) * t415;
t358 = (t369 * t488 + t378 * t414) * t415;
t357 = (-t368 * t482 + t374 * t414) * t411;
t356 = (t368 * t490 + t373 * t414) * t411;
t355 = -t497 * t512 + (-t371 * t468 + t414 * t422) * t502;
t354 = -t499 * t514 + (-t370 * t469 + t414 * t420) * t504;
t353 = -t501 * t516 + (-t369 * t470 + t414 * t418) * t506;
t352 = -t509 * t518 + (-t368 * t474 + t412 * t414) * t510;
t351 = (t382 * t391 + t383 * t387) * t417;
t350 = (t380 * t390 + t381 * t386) * t416;
t349 = (t378 * t389 + t379 * t385) * t415;
t348 = (t373 * t388 + t374 * t384) * t411;
t347 = (t371 * t383 - t377 * t476) * t417;
t346 = (t371 * t382 + t377 * t484) * t417;
t345 = (t370 * t381 - t376 * t478) * t416;
t344 = (t370 * t380 + t376 * t486) * t416;
t343 = (t369 * t379 - t375 * t480) * t415;
t342 = (t369 * t378 + t375 * t488) * t415;
t341 = (t368 * t374 - t372 * t482) * t411;
t340 = (t368 * t373 + t372 * t490) * t411;
t339 = t396 * t497 + (t371 * t422 - t377 * t468) * t502;
t338 = t395 * t499 + (t370 * t420 - t376 * t469) * t504;
t337 = t394 * t501 + (t369 * t418 - t375 * t470) * t506;
t336 = t392 * t509 + (t368 * t412 - t372 * t474) * t510;
t335 = t351 * t414 + t367 * t371;
t334 = t350 * t414 + t366 * t370;
t333 = t349 * t414 + t365 * t369;
t332 = t348 * t414 + t364 * t368;
t331 = t351 * t371 + t367 * t377;
t330 = t350 * t370 + t366 * t376;
t329 = t349 * t369 + t365 * t375;
t328 = t348 * t368 + t364 * t372;
t1 = [t357 * t526 + t359 * t524 + t361 * t522 + t363 * t520 + m(4) + (-t341 * t483 - t343 * t481 - t345 * t479 - t347 * t477) * t457, t357 * t527 + t359 * t525 + t361 * t523 + t363 * t521 + (t341 * t491 + t343 * t489 + t345 * t487 + t347 * t485) * t457, t357 * t475 + t359 * t473 + t361 * t472 + t363 * t471 + (-t347 * t463 - t383 * t513 - t345 * t464 - t381 * t515 - t343 * t465 - t379 * t517 - t341 * t466 - t374 * t519 + (-t405 * t495 - t406 * t494 - t407 * t493 - t408 * t492) * t457) * t457, t341 * t364 + t343 * t365 + t345 * t366 + t347 * t367 + t357 * t348 + t359 * t349 + t361 * t350 + t363 * t351 + t462; t356 * t526 + t358 * t524 + t360 * t522 + t362 * t520 + (-t340 * t483 - t342 * t481 - t344 * t479 - t346 * t477) * t457, t356 * t527 + t358 * t525 + t360 * t523 + t362 * t521 + m(4) + (t340 * t491 + t342 * t489 + t344 * t487 + t346 * t485) * t457, t356 * t475 + t358 * t473 + t360 * t472 + t362 * t471 + (-t346 * t463 - t382 * t513 - t344 * t464 - t380 * t515 - t342 * t465 - t378 * t517 - t340 * t466 - t373 * t519 + (t401 * t495 + t402 * t494 + t403 * t493 + t404 * t492) * t457) * t457, t340 * t364 + t342 * t365 + t344 * t366 + t346 * t367 + t356 * t348 + t358 * t349 + t360 * t350 + t362 * t351 + t467; t352 * t526 + t353 * t524 + t354 * t522 + t355 * t520 + (-t336 * t483 - t337 * t481 - t338 * t479 - t339 * t477) * t457, t352 * t527 + t353 * t525 + t354 * t523 + t355 * t521 + (t336 * t491 + t337 * t489 + t338 * t487 + t339 * t485) * t457, t352 * t475 + t353 * t473 + t354 * t472 + t355 * t471 + m(4) + (-t336 * t466 - t337 * t465 - t338 * t464 - t339 * t463 - t424 * t532 * t393 - t434 * t531 * t397 - t436 * t530 * t398 - t438 * t529 * t399 + ((Ifges(3,3) * t422 - t396 * t463) * t422 + (Ifges(3,3) * t420 - t395 * t464) * t420 + (Ifges(3,3) * t418 - t394 * t465) * t418 + (Ifges(3,3) * t412 - t392 * t466) * t412) * t457) * t457, t336 * t364 + t337 * t365 + t338 * t366 + t339 * t367 + t352 * t348 + t353 * t349 + t354 * t350 + t355 * t351; t332 * t526 + t333 * t524 + t334 * t522 + t335 * t520 + (-t328 * t483 - t329 * t481 - t330 * t479 - t331 * t477) * t457 + t462, t332 * t527 + t333 * t525 + t334 * t523 + t335 * t521 + (t328 * t491 + t329 * t489 + t330 * t487 + t331 * t485) * t457 + t467, t332 * t475 + t333 * t473 + t334 * t472 + t335 * t471 + (-t331 * t463 + (-t351 * t512 + t367 * t396) * t422 - t330 * t464 + (-t350 * t514 + t366 * t395) * t420 - t329 * t465 + (-t349 * t516 + t365 * t394) * t418 - t328 * t466 + (-t348 * t518 + t364 * t392) * t412) * t457, t328 * t364 + t329 * t365 + t330 * t366 + t331 * t367 + t332 * t348 + t333 * t349 + t334 * t350 + t335 * t351 + Ifges(4,3);];
MX  = t1;
