% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3RPRRR12V1G1A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x14]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-08-06 18:21
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3RPRRR12V1G1A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_regmin: pkin has to be [6x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G1A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:21:34
% EndTime: 2020-08-06 18:21:35
% DurationCPUTime: 0.94s
% Computational Cost: add. (1668->167), mult. (2133->330), div. (441->14), fcn. (2277->18), ass. (0->166)
t444 = pkin(1) + pkin(5);
t432 = sin(qJ(3,3));
t405 = t432 * pkin(3) + qJ(2,3);
t419 = pkin(6) + t444;
t433 = sin(qJ(1,3));
t439 = cos(qJ(1,3));
t381 = t433 * t405 + t419 * t439;
t384 = -t405 * t439 + t419 * t433;
t429 = legFrame(3,3);
t409 = sin(t429);
t412 = cos(t429);
t375 = t381 * t412 - t409 * t384;
t399 = 0.1e1 / t405;
t369 = t375 * t399;
t434 = sin(qJ(3,2));
t406 = t434 * pkin(3) + qJ(2,2);
t435 = sin(qJ(1,2));
t441 = cos(qJ(1,2));
t382 = t435 * t406 + t419 * t441;
t385 = -t406 * t441 + t419 * t435;
t430 = legFrame(2,3);
t410 = sin(t430);
t413 = cos(t430);
t376 = t382 * t413 - t410 * t385;
t401 = 0.1e1 / t406;
t370 = t376 * t401;
t436 = sin(qJ(3,1));
t407 = -t436 * pkin(3) - qJ(2,1);
t437 = sin(qJ(1,1));
t443 = cos(qJ(1,1));
t383 = -t437 * t407 + t419 * t443;
t386 = t407 * t443 + t419 * t437;
t431 = legFrame(1,3);
t411 = sin(t431);
t414 = cos(t431);
t377 = t383 * t414 - t411 * t386;
t403 = 0.1e1 / t407;
t505 = t377 * t403;
t378 = t409 * t381 + t384 * t412;
t372 = t378 * t399;
t400 = 0.1e1 / t405 ^ 2;
t504 = t378 * t400;
t379 = t410 * t382 + t385 * t413;
t373 = t379 * t401;
t402 = 0.1e1 / t406 ^ 2;
t503 = t379 * t402;
t380 = t411 * t383 + t386 * t414;
t502 = t380 * t403;
t404 = 0.1e1 / t407 ^ 2;
t501 = t380 * t404;
t393 = -t409 * t433 + t412 * t439;
t387 = t393 ^ 2;
t500 = t387 * t400;
t394 = -t410 * t435 + t413 * t441;
t388 = t394 ^ 2;
t499 = t388 * t402;
t395 = -t411 * t437 + t414 * t443;
t389 = t395 ^ 2;
t498 = t389 * t404;
t396 = t409 * t439 + t412 * t433;
t390 = t396 ^ 2;
t497 = t390 * t400;
t397 = t410 * t441 + t413 * t435;
t391 = t397 ^ 2;
t496 = t391 * t402;
t398 = t411 * t443 + t414 * t437;
t392 = t398 ^ 2;
t495 = t392 * t404;
t494 = t393 * t399;
t493 = t394 * t401;
t492 = t395 * t403;
t491 = t396 * t399;
t490 = t396 * t400;
t489 = t397 * t401;
t488 = t397 * t402;
t487 = t398 * t403;
t486 = t398 * t404;
t485 = t399 * t444;
t438 = cos(qJ(3,3));
t426 = t438 ^ 2;
t484 = t400 * t426;
t483 = t401 * t444;
t440 = cos(qJ(3,2));
t427 = t440 ^ 2;
t482 = t402 * t427;
t481 = t403 * t444;
t442 = cos(qJ(3,1));
t428 = t442 ^ 2;
t480 = t404 * t428;
t479 = 0.1e1 / t432 * t438;
t478 = 0.1e1 / t434 * t440;
t477 = 0.1e1 / t436 * t442;
t476 = pkin(1) * t494;
t475 = pkin(1) * t493;
t474 = pkin(1) * t492;
t473 = pkin(1) * t491;
t472 = pkin(1) * t489;
t471 = pkin(1) * t487;
t470 = qJ(2,3) * t500;
t469 = qJ(2,2) * t499;
t468 = qJ(2,1) * t498;
t467 = qJ(2,3) * t497;
t466 = qJ(2,2) * t496;
t465 = qJ(2,1) * t495;
t464 = t393 * t490;
t463 = t394 * t488;
t462 = t395 * t486;
t461 = t399 * t479;
t460 = t400 * t432 * t438;
t459 = t401 * t478;
t458 = t402 * t434 * t440;
t457 = t403 * t477;
t456 = t404 * t436 * t442;
t455 = qJ(2,3) * t464;
t454 = qJ(2,2) * t463;
t453 = qJ(2,1) * t462;
t452 = -t492 + t493 + t494;
t451 = -t487 + t489 + t491;
t421 = 0.1e1 / t432 ^ 2;
t423 = 0.1e1 / t434 ^ 2;
t425 = 0.1e1 / t436 ^ 2;
t450 = t426 * t421 + t427 * t423 + t428 * t425;
t449 = -t375 * t461 - t376 * t459 + t377 * t457;
t448 = -t378 * t461 - t379 * t459 + t380 * t457;
t341 = t393 * t461 + t394 * t459 - t395 * t457;
t342 = t396 * t461 + t397 * t459 - t398 * t457;
t447 = pkin(1) ^ 2;
t445 = 0.1e1 / pkin(3);
t417 = qJ(2,1) ^ 2 + t447;
t416 = qJ(2,2) ^ 2 + t447;
t415 = qJ(2,3) ^ 2 + t447;
t368 = t451 * t445;
t367 = t452 * t445;
t366 = -t502 + 0.2e1 * t471;
t365 = -t502 + t471;
t364 = t373 - 0.2e1 * t472;
t363 = t373 - t472;
t362 = t372 - 0.2e1 * t473;
t361 = t372 - t473;
t360 = -t505 + 0.2e1 * t474;
t359 = -t505 + t474;
t358 = t370 - 0.2e1 * t475;
t357 = t370 - t475;
t356 = t369 - 0.2e1 * t476;
t355 = t369 - t476;
t354 = -t398 * t481 + t502;
t353 = t397 * t483 - t373;
t352 = t396 * t485 - t372;
t351 = -t395 * t481 + t505;
t350 = t394 * t483 - t370;
t349 = t393 * t485 - t369;
t348 = (-pkin(1) * t380 + t398 * t417) * t403;
t347 = (-pkin(1) * t377 + t395 * t417) * t403;
t346 = (-pkin(1) * t379 + t397 * t416) * t401;
t345 = (-pkin(1) * t376 + t394 * t416) * t401;
t344 = (-pkin(1) * t378 + t396 * t415) * t399;
t343 = (-pkin(1) * t375 + t393 * t415) * t399;
t340 = t342 * t445;
t339 = t341 * t445;
t338 = t462 + t463 + t464;
t337 = 0.2e1 * t453 + 0.2e1 * t454 + 0.2e1 * t455;
t336 = t426 * t464 + t427 * t463 + t428 * t462;
t335 = 0.2e1 * t438 * t455 + 0.2e1 * t440 * t454 + 0.2e1 * t442 * t453;
t334 = 0.2e1 * t432 * t455 + 0.2e1 * t434 * t454 + 0.2e1 * t436 * t453;
t333 = -0.2e1 * t393 * t396 * t460 - 0.2e1 * t394 * t397 * t458 - 0.2e1 * t395 * t398 * t456;
t1 = [t498 + t499 + t500, 0, 0, (-t360 * t403 + t377 * t404) * t395 + (t358 * t401 + t376 * t402) * t394 + (t356 * t399 + t375 * t400) * t393, 0.2e1 * t468 + 0.2e1 * t469 + 0.2e1 * t470, -(-t347 * t395 + t359 * t377) * t403 + (t345 * t394 + t357 * t376) * t401 + (t343 * t393 + t355 * t375) * t399, t387 * t484 + t388 * t482 + t389 * t480, -0.2e1 * t387 * t460 - 0.2e1 * t388 * t458 - 0.2e1 * t389 * t456, 0, 0, 0, 0.2e1 * t432 * t470 + 0.2e1 * t434 * t469 + 0.2e1 * t436 * t468, 0.2e1 * t438 * t470 + 0.2e1 * t440 * t469 + 0.2e1 * t442 * t468, 1; t338, 0, 0, t362 * t494 + t364 * t493 - t366 * t492 + t375 * t490 + t376 * t488 + t377 * t486, t337, -(-t348 * t395 + t365 * t377) * t403 + (t346 * t394 + t363 * t376) * t401 + (t344 * t393 + t361 * t375) * t399, t336, t333, 0, 0, 0, t334, t335, 0; 0, 0, 0, t341, 0, -t341 * pkin(1) - t449, 0, 0, -t339, t367, 0, (t341 * t444 + t449) * t445, (-t452 * t444 + t369 + t370 - t505) * t445, 0; t338, 0, 0, t356 * t491 + t358 * t489 - t360 * t487 + t393 * t504 + t394 * t503 + t395 * t501, t337, -(-t347 * t398 + t359 * t380) * t403 + (t345 * t397 + t357 * t379) * t401 + (t343 * t396 + t355 * t378) * t399, t336, t333, 0, 0, 0, t334, t335, 0; t495 + t496 + t497, 0, 0, (-t366 * t403 + t501) * t398 + (t364 * t401 + t503) * t397 + (t362 * t399 + t504) * t396, 0.2e1 * t465 + 0.2e1 * t466 + 0.2e1 * t467, -(-t348 * t398 + t365 * t380) * t403 + (t346 * t397 + t363 * t379) * t401 + (t344 * t396 + t361 * t378) * t399, t390 * t484 + t391 * t482 + t392 * t480, -0.2e1 * t390 * t460 - 0.2e1 * t391 * t458 - 0.2e1 * t392 * t456, 0, 0, 0, 0.2e1 * t432 * t467 + 0.2e1 * t434 * t466 + 0.2e1 * t436 * t465, 0.2e1 * t438 * t467 + 0.2e1 * t440 * t466 + 0.2e1 * t442 * t465, 1; 0, 0, 0, t342, 0, -t342 * pkin(1) - t448, 0, 0, -t340, t368, 0, (t342 * t444 + t448) * t445, (-t451 * t444 + t372 + t373 - t502) * t445, 0; 0, 0, 0, t341, 0, t355 * t479 + t357 * t478 + t359 * t477, 0, 0, -t339, t367, 0, (t349 * t479 + t350 * t478 + t351 * t477) * t445, (-t349 - t350 - t351) * t445, 0; 0, 0, 0, t342, 0, t361 * t479 + t363 * t478 + t365 * t477, 0, 0, -t340, t368, 0, (t352 * t479 + t353 * t478 + t354 * t477) * t445, (-t352 - t353 - t354) * t445, 0; 0, 0, 0, 0, 0, t450, 0, 0, 0, 0, (t421 + t423 + t425) / pkin(3) ^ 2, -0.2e1 * t450 * t445, 0.2e1 * (t477 + t478 + t479) * t445, 1;];
tau_reg  = t1;
