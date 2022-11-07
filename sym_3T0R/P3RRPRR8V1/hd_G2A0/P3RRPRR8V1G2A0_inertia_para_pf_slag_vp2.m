% Calculate inertia matrix for parallel robot
% P3RRPRR8V1G2A0
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,d1,d4,theta3]';
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
% Datum: 2022-11-04 17:05
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(5,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR8V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:04:25
% EndTime: 2022-11-04 17:04:26
% DurationCPUTime: 0.95s
% Computational Cost: add. (2682->207), mult. (4131->345), div. (294->9), fcn. (3012->23), ass. (0->159)
t493 = m(3) * pkin(1);
t492 = 2 * mrSges(3,3);
t418 = cos(pkin(5));
t413 = t418 ^ 2;
t417 = sin(pkin(5));
t422 = Ifges(3,2) - Ifges(3,1);
t484 = t417 * mrSges(3,1);
t491 = 0.4e1 * Ifges(3,4) * t413 - 0.2e1 * (pkin(1) * mrSges(3,2) + t422 * t417) * t418 - 0.2e1 * pkin(1) * t484 + (2 * Ifges(2,4)) - 0.2e1 * Ifges(3,4);
t490 = (m(3) * qJ(3,1));
t489 = (m(3) * qJ(3,2));
t488 = (m(3) * qJ(3,3));
t487 = pkin(2) * t417;
t486 = mrSges(3,2) * t417;
t485 = Ifges(3,4) * t417;
t426 = sin(qJ(2,3));
t427 = sin(qJ(1,3));
t432 = cos(qJ(2,3));
t423 = legFrame(3,2);
t407 = sin(t423);
t457 = t407 * t487;
t410 = cos(t423);
t460 = t410 * t487;
t397 = t418 * pkin(2) + pkin(1);
t465 = t407 * t397;
t468 = t397 * t410;
t362 = (-t427 * t465 + t460) * t432 + (t427 * t457 + t468) * t426;
t441 = t397 * t432 - t426 * t487;
t382 = 0.1e1 / t441;
t483 = t362 * t382;
t428 = sin(qJ(2,2));
t429 = sin(qJ(1,2));
t434 = cos(qJ(2,2));
t424 = legFrame(2,2);
t408 = sin(t424);
t456 = t408 * t487;
t411 = cos(t424);
t459 = t411 * t487;
t464 = t408 * t397;
t467 = t397 * t411;
t363 = (-t429 * t464 + t459) * t434 + (t429 * t456 + t467) * t428;
t440 = t397 * t434 - t428 * t487;
t383 = 0.1e1 / t440;
t482 = t363 * t383;
t430 = sin(qJ(2,1));
t431 = sin(qJ(1,1));
t436 = cos(qJ(2,1));
t425 = legFrame(1,2);
t409 = sin(t425);
t455 = t409 * t487;
t412 = cos(t425);
t458 = t412 * t487;
t463 = t409 * t397;
t466 = t397 * t412;
t364 = (-t431 * t463 + t458) * t436 + (t431 * t455 + t466) * t430;
t439 = t397 * t436 - t430 * t487;
t384 = 0.1e1 / t439;
t481 = t364 * t384;
t365 = (t427 * t468 + t457) * t432 + (-t427 * t460 + t465) * t426;
t480 = t365 * t382;
t366 = (t429 * t467 + t456) * t434 + (-t429 * t459 + t464) * t428;
t479 = t366 * t383;
t367 = (t431 * t466 + t455) * t436 + (-t431 * t458 + t463) * t430;
t478 = t367 * t384;
t389 = -mrSges(3,1) * t418 + t486 - t493;
t393 = t418 * mrSges(3,2) + t484;
t372 = t389 * t432 + t426 * t393;
t477 = t372 * t382;
t373 = t389 * t434 + t428 * t393;
t476 = t373 * t383;
t374 = t389 * t436 + t430 * t393;
t475 = t374 * t384;
t444 = t432 * pkin(1) + pkin(2) * cos(qJ(2,3) + pkin(5));
t390 = 0.1e1 / t444;
t474 = t390 * t407;
t473 = t390 * t410;
t443 = t434 * pkin(1) + pkin(2) * cos(qJ(2,2) + pkin(5));
t391 = 0.1e1 / t443;
t472 = t391 * t408;
t471 = t391 * t411;
t442 = t436 * pkin(1) + pkin(2) * cos(qJ(2,1) + pkin(5));
t392 = 0.1e1 / t442;
t470 = t392 * t409;
t469 = t392 * t412;
t395 = t422 * t413;
t462 = (-0.2e1 * t486 + t493) * pkin(1);
t461 = 0.2e1 * pkin(1) * mrSges(3,1);
t401 = mrSges(3,2) * qJ(3,3) - Ifges(3,6);
t404 = mrSges(3,1) * qJ(3,3) - Ifges(3,5);
t368 = (t401 * t417 - t404 * t418 + Ifges(2,5) + (-mrSges(3,3) - t488) * pkin(1)) * t426 - (t401 * t418 + t404 * t417 - Ifges(2,6)) * t432;
t419 = pkin(4) + qJ(3,3);
t414 = 0.1e1 / t419;
t454 = t368 * t382 * t414;
t402 = mrSges(3,2) * qJ(3,2) - Ifges(3,6);
t405 = mrSges(3,1) * qJ(3,2) - Ifges(3,5);
t369 = (t402 * t417 - t405 * t418 + Ifges(2,5) + (-mrSges(3,3) - t489) * pkin(1)) * t428 - (t402 * t418 + t405 * t417 - Ifges(2,6)) * t434;
t420 = pkin(4) + qJ(3,2);
t415 = 0.1e1 / t420;
t453 = t369 * t383 * t415;
t403 = mrSges(3,2) * qJ(3,1) - Ifges(3,6);
t406 = mrSges(3,1) * qJ(3,1) - Ifges(3,5);
t370 = (t403 * t417 - t406 * t418 + Ifges(2,5) + (-mrSges(3,3) - t490) * pkin(1)) * t430 - (t403 * t418 + t406 * t417 - Ifges(2,6)) * t436;
t421 = pkin(4) + qJ(3,1);
t416 = 0.1e1 / t421;
t452 = t370 * t384 * t416;
t451 = t368 * t474;
t450 = t369 * t472;
t449 = t370 * t470;
t448 = t368 * t473;
t447 = t369 * t471;
t446 = t370 * t469;
t445 = -0.2e1 * t418 * t485 + Ifges(2,1) + Ifges(3,2) + Ifges(1,3) - t395;
t437 = cos(qJ(1,1));
t435 = cos(qJ(1,2));
t433 = cos(qJ(1,3));
t388 = t418 * t461 + Ifges(2,3) + Ifges(3,3) + t462;
t387 = t430 * t397 + t436 * t487;
t386 = t428 * t397 + t434 * t487;
t385 = t426 * t397 + t432 * t487;
t381 = t431 * t421 + t442 * t437;
t380 = t429 * t420 + t443 * t435;
t379 = t427 * t419 + t444 * t433;
t377 = -t437 * t421 + t439 * t431;
t376 = -t435 * t420 + t440 * t429;
t375 = -t433 * t419 + t441 * t427;
t371 = 0.2e1 * t395 + (t461 + 0.4e1 * t485) * t418 - Ifges(2,1) + Ifges(2,2) + t462 - t422;
t361 = -t377 * t409 + t412 * t387;
t360 = t377 * t412 + t409 * t387;
t359 = -t376 * t408 + t411 * t386;
t358 = t376 * t411 + t408 * t386;
t357 = -t375 * t407 + t410 * t385;
t356 = t375 * t410 + t407 * t385;
t355 = (m(3) * t381 + t374 * t437) * t416;
t354 = (m(3) * t380 + t373 * t435) * t415;
t353 = (m(3) * t379 + t372 * t433) * t414;
t352 = (t371 * t436 + t430 * t491) * t436 + (t492 + t490) * qJ(3,1) + t445;
t351 = (t371 * t434 + t428 * t491) * t434 + (t492 + t489) * qJ(3,2) + t445;
t350 = (t371 * t432 + t426 * t491) * t432 + (t492 + t488) * qJ(3,3) + t445;
t349 = (t352 * t437 + t374 * t381) * t416;
t348 = (t351 * t435 + t373 * t380) * t415;
t347 = (t350 * t433 + t372 * t379) * t414;
t346 = t367 * t452 + t388 * t470;
t345 = t366 * t453 + t388 * t472;
t344 = t365 * t454 + t388 * t474;
t343 = t364 * t452 + t388 * t469;
t342 = t363 * t453 + t388 * t471;
t341 = t362 * t454 + t388 * t473;
t340 = (m(3) * t360 + t367 * t475) * t416;
t339 = (m(3) * t358 + t366 * t476) * t415;
t338 = (m(3) * t356 + t365 * t477) * t414;
t337 = (m(3) * t361 + t364 * t475) * t416;
t336 = (m(3) * t359 + t363 * t476) * t415;
t335 = (m(3) * t357 + t362 * t477) * t414;
t334 = t449 + (t352 * t478 + t360 * t374) * t416;
t333 = t450 + (t351 * t479 + t358 * t373) * t415;
t332 = t451 + (t350 * t480 + t356 * t372) * t414;
t331 = t446 + (t352 * t481 + t361 * t374) * t416;
t330 = t447 + (t351 * t482 + t359 * t373) * t415;
t329 = t448 + (t350 * t483 + t357 * t372) * t414;
t1 = [t344 * t474 + t345 * t472 + t346 * t470 + m(4) + (t334 * t478 + t340 * t360) * t416 + (t333 * t479 + t339 * t358) * t415 + (t332 * t480 + t338 * t356) * t414, t344 * t473 + t345 * t471 + t346 * t469 + (t334 * t481 + t340 * t361) * t416 + (t333 * t482 + t339 * t359) * t415 + (t332 * t483 + t338 * t357) * t414, (t334 * t437 + t340 * t381) * t416 + (t333 * t435 + t339 * t380) * t415 + (t332 * t433 + t338 * t379) * t414; t341 * t474 + t342 * t472 + t343 * t470 + (t331 * t478 + t337 * t360) * t416 + (t330 * t479 + t336 * t358) * t415 + (t329 * t480 + t335 * t356) * t414, t341 * t473 + t342 * t471 + t343 * t469 + m(4) + (t331 * t481 + t337 * t361) * t416 + (t330 * t482 + t336 * t359) * t415 + (t329 * t483 + t335 * t357) * t414, (t331 * t437 + t337 * t381) * t416 + (t330 * t435 + t336 * t380) * t415 + (t329 * t433 + t335 * t379) * t414; (t349 * t478 + t355 * t360 + t437 * t449) * t416 + (t348 * t479 + t354 * t358 + t435 * t450) * t415 + (t347 * t480 + t353 * t356 + t433 * t451) * t414, (t349 * t481 + t355 * t361 + t437 * t446) * t416 + (t348 * t482 + t354 * t359 + t435 * t447) * t415 + (t347 * t483 + t353 * t357 + t433 * t448) * t414, m(4) + (t349 * t437 + t355 * t381) * t416 + (t348 * t435 + t354 * t380) * t415 + (t347 * t433 + t353 * t379) * t414;];
MX  = t1;
