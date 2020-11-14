% Calculate inertia matrix for parallel robot
% P4PRRR1G1A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRR1G1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:12:30
% EndTime: 2020-03-02 20:12:32
% DurationCPUTime: 1.67s
% Computational Cost: add. (4685->165), mult. (3107->314), div. (720->6), fcn. (3392->34), ass. (0->164)
t431 = pkin(7) + qJ(2,4);
t411 = qJ(3,4) + t431;
t401 = sin(t411);
t402 = cos(t411);
t409 = sin(t431);
t410 = cos(t431);
t508 = 0.1e1 / (t401 * t410 - t409 * t402);
t432 = pkin(7) + qJ(2,3);
t418 = qJ(3,3) + t432;
t403 = sin(t418);
t406 = cos(t418);
t412 = sin(t432);
t415 = cos(t432);
t507 = 0.1e1 / (t403 * t415 - t412 * t406);
t433 = pkin(7) + qJ(2,2);
t419 = qJ(3,2) + t433;
t404 = sin(t419);
t407 = cos(t419);
t413 = sin(t433);
t416 = cos(t433);
t506 = 0.1e1 / (t404 * t416 - t413 * t407);
t434 = pkin(7) + qJ(2,1);
t420 = qJ(3,1) + t434;
t405 = sin(t420);
t408 = cos(t420);
t414 = sin(t434);
t417 = cos(t434);
t505 = 0.1e1 / (t405 * t417 - t414 * t408);
t435 = legFrame(4,3);
t421 = sin(t435);
t425 = cos(t435);
t376 = t425 * t401 + t421 * t402;
t356 = pkin(2) * (t425 * t409 + t421 * t410) + t376 * pkin(3);
t504 = t356 * t508;
t377 = -t401 * t421 + t425 * t402;
t357 = -pkin(2) * (t409 * t421 - t425 * t410) + t377 * pkin(3);
t503 = t357 * t508;
t436 = legFrame(3,3);
t422 = sin(t436);
t426 = cos(t436);
t378 = t426 * t403 + t422 * t406;
t358 = pkin(2) * (t426 * t412 + t422 * t415) + t378 * pkin(3);
t502 = t358 * t507;
t437 = legFrame(2,3);
t423 = sin(t437);
t427 = cos(t437);
t380 = t427 * t404 + t423 * t407;
t359 = pkin(2) * (t427 * t413 + t423 * t416) + t380 * pkin(3);
t501 = t359 * t506;
t438 = legFrame(1,3);
t424 = sin(t438);
t428 = cos(t438);
t382 = t428 * t405 + t424 * t408;
t360 = pkin(2) * (t428 * t414 + t424 * t417) + t382 * pkin(3);
t500 = t360 * t505;
t379 = -t403 * t422 + t426 * t406;
t361 = -pkin(2) * (t412 * t422 - t426 * t415) + t379 * pkin(3);
t499 = t361 * t507;
t381 = -t404 * t423 + t427 * t407;
t362 = -pkin(2) * (t413 * t423 - t427 * t416) + t381 * pkin(3);
t498 = t362 * t506;
t383 = -t405 * t424 + t428 * t408;
t363 = -pkin(2) * (t414 * t424 - t428 * t417) + t383 * pkin(3);
t497 = t363 * t505;
t496 = t508 * t376;
t495 = t508 * t377;
t464 = (mrSges(3,1) * cos(qJ(3,4)) - mrSges(3,2) * sin(qJ(3,4))) * pkin(2);
t475 = m(3) * pkin(2) ^ 2 + Ifges(2,3) + Ifges(3,3);
t384 = 0.2e1 * t464 + t475;
t494 = t508 * t384;
t396 = Ifges(3,3) + t464;
t493 = t508 * t396;
t459 = 0.1e1 / pkin(3);
t492 = t508 * t459;
t491 = t507 * t378;
t490 = t507 * t379;
t463 = (mrSges(3,1) * cos(qJ(3,3)) - mrSges(3,2) * sin(qJ(3,3))) * pkin(2);
t385 = 0.2e1 * t463 + t475;
t489 = t507 * t385;
t397 = Ifges(3,3) + t463;
t488 = t507 * t397;
t487 = t506 * t380;
t486 = t506 * t381;
t462 = (mrSges(3,1) * cos(qJ(3,2)) - mrSges(3,2) * sin(qJ(3,2))) * pkin(2);
t386 = 0.2e1 * t462 + t475;
t485 = t506 * t386;
t398 = Ifges(3,3) + t462;
t484 = t506 * t398;
t483 = t505 * t382;
t482 = t505 * t383;
t461 = (mrSges(3,1) * cos(qJ(3,1)) - mrSges(3,2) * sin(qJ(3,1))) * pkin(2);
t387 = 0.2e1 * t461 + t475;
t481 = t505 * t387;
t399 = Ifges(3,3) + t461;
t480 = t505 * t399;
t479 = t507 * t459;
t478 = t506 * t459;
t477 = t505 * t459;
t460 = 0.1e1 / pkin(2);
t476 = t459 * t460;
t474 = Ifges(3,3) * t492;
t473 = Ifges(3,3) * t479;
t472 = Ifges(3,3) * t478;
t471 = Ifges(3,3) * t477;
t470 = t396 * t492;
t469 = t397 * t479;
t468 = t398 * t478;
t467 = t399 * t477;
t448 = xP(4);
t429 = sin(t448);
t430 = cos(t448);
t449 = mrSges(4,2);
t450 = mrSges(4,1);
t466 = -t429 * t449 + t430 * t450;
t465 = -t429 * t450 - t430 * t449;
t458 = koppelP(1,1);
t457 = koppelP(2,1);
t456 = koppelP(3,1);
t455 = koppelP(4,1);
t454 = koppelP(1,2);
t453 = koppelP(2,2);
t452 = koppelP(3,2);
t451 = koppelP(4,2);
t395 = -t429 * t454 + t430 * t458;
t394 = -t429 * t453 + t430 * t457;
t393 = -t429 * t452 + t430 * t456;
t392 = -t429 * t451 + t430 * t455;
t391 = -t429 * t458 - t430 * t454;
t390 = -t429 * t457 - t430 * t453;
t389 = -t429 * t456 - t430 * t452;
t388 = -t429 * t455 - t430 * t451;
t355 = (t382 * t395 + t383 * t391) * t460 * t505;
t354 = (t380 * t394 + t381 * t390) * t460 * t506;
t353 = (t378 * t393 + t379 * t389) * t460 * t507;
t352 = (t376 * t392 + t377 * t388) * t460 * t508;
t351 = (-t363 * t471 + t383 * t480) * t460;
t350 = (-t360 * t471 + t382 * t480) * t460;
t349 = (-t362 * t472 + t381 * t484) * t460;
t348 = (-t359 * t472 + t380 * t484) * t460;
t347 = (-t361 * t473 + t379 * t488) * t460;
t346 = (-t358 * t473 + t378 * t488) * t460;
t345 = (-t357 * t474 + t377 * t493) * t460;
t344 = (-t356 * t474 + t376 * t493) * t460;
t343 = (-t363 * t467 + t383 * t481) * t460;
t342 = (-t360 * t467 + t382 * t481) * t460;
t341 = (-t362 * t468 + t381 * t485) * t460;
t340 = (-t359 * t468 + t380 * t485) * t460;
t339 = (-t361 * t469 + t379 * t489) * t460;
t338 = (-t358 * t469 + t378 * t489) * t460;
t337 = (-t357 * t470 + t377 * t494) * t460;
t336 = (-t356 * t470 + t376 * t494) * t460;
t335 = (t360 * t395 + t363 * t391) * t505 * t476;
t334 = (t359 * t394 + t362 * t390) * t506 * t476;
t333 = (t358 * t393 + t361 * t389) * t507 * t476;
t332 = (t356 * t392 + t357 * t388) * t508 * t476;
t331 = -t335 * Ifges(3,3) + t355 * t399;
t330 = -t334 * Ifges(3,3) + t354 * t398;
t329 = -t333 * Ifges(3,3) + t353 * t397;
t328 = -t332 * Ifges(3,3) + t352 * t396;
t327 = -t335 * t399 + t355 * t387;
t326 = -t334 * t398 + t354 * t386;
t325 = -t333 * t397 + t353 * t385;
t324 = -t332 * t396 + t352 * t384;
t1 = [m(4) + (t337 * t495 + t339 * t490 + t341 * t486 + t343 * t482 + (-t345 * t503 - t347 * t499 - t349 * t498 - t351 * t497) * t459) * t460, (t337 * t496 + t339 * t491 + t341 * t487 + t343 * t483 + (-t345 * t504 - t347 * t502 - t349 * t501 - t351 * t500) * t459) * t460, 0, -t345 * t332 - t347 * t333 - t349 * t334 - t351 * t335 + t337 * t352 + t339 * t353 + t341 * t354 + t343 * t355 + t465; (t336 * t495 + t338 * t490 + t340 * t486 + t342 * t482 + (-t344 * t503 - t346 * t499 - t348 * t498 - t350 * t497) * t459) * t460, m(4) + (t336 * t496 + t338 * t491 + t340 * t487 + t342 * t483 + (-t344 * t504 - t346 * t502 - t348 * t501 - t350 * t500) * t459) * t460, 0, -t344 * t332 - t346 * t333 - t348 * t334 - t350 * t335 + t336 * t352 + t338 * t353 + t340 * t354 + t342 * t355 + t466; 0, 0, (4 * m(1)) + (4 * m(2)) + 0.4e1 * m(3) + m(4), 0; (t324 * t495 + t325 * t490 + t326 * t486 + t327 * t482 + (-t328 * t503 - t329 * t499 - t330 * t498 - t331 * t497) * t459) * t460 + t465, (t324 * t496 + t325 * t491 + t326 * t487 + t327 * t483 + (-t328 * t504 - t329 * t502 - t330 * t501 - t331 * t500) * t459) * t460 + t466, 0, t324 * t352 + t325 * t353 + t326 * t354 + t327 * t355 - t328 * t332 - t329 * t333 - t330 * t334 - t331 * t335 + Ifges(4,3);];
MX  = t1;
