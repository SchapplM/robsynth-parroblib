% Calculate inertia matrix for parallel robot
% P4PRRR1G1P1A0
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
% rSges [4x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [4x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% MX [4x4]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-02 20:15
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P4PRRR1G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,4),zeros(4,3),zeros(4,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [4 1]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: xP has to be [4x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 4]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x4] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [4 3]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [4x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [4 3]), ...
  'P4PRRR1G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [4x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-02 20:10:30
% EndTime: 2020-03-02 20:10:32
% DurationCPUTime: 1.66s
% Computational Cost: add. (6127->180), mult. (4266->331), div. (720->6), fcn. (3824->26), ass. (0->156)
t448 = pkin(7) + qJ(2,4);
t428 = qJ(3,4) + t448;
t416 = sin(t428);
t417 = cos(t428);
t426 = sin(t448);
t427 = cos(t448);
t507 = 0.1e1 / (t416 * t427 - t426 * t417);
t449 = pkin(7) + qJ(2,3);
t435 = qJ(3,3) + t449;
t418 = sin(t435);
t421 = cos(t435);
t429 = sin(t449);
t432 = cos(t449);
t506 = 0.1e1 / (t418 * t432 - t429 * t421);
t450 = pkin(7) + qJ(2,2);
t436 = qJ(3,2) + t450;
t419 = sin(t436);
t422 = cos(t436);
t430 = sin(t450);
t433 = cos(t450);
t505 = 0.1e1 / (t419 * t433 - t430 * t422);
t451 = pkin(7) + qJ(2,1);
t437 = qJ(3,1) + t451;
t420 = sin(t437);
t423 = cos(t437);
t431 = sin(t451);
t434 = cos(t451);
t504 = 0.1e1 / (t420 * t434 - t431 * t423);
t503 = m(3) * pkin(2);
t452 = legFrame(4,3);
t438 = sin(t452);
t442 = cos(t452);
t388 = t442 * t416 + t438 * t417;
t368 = pkin(2) * (t442 * t426 + t438 * t427) + t388 * pkin(3);
t502 = t368 * t507;
t389 = -t416 * t438 + t442 * t417;
t369 = -pkin(2) * (t426 * t438 - t442 * t427) + t389 * pkin(3);
t501 = t369 * t507;
t453 = legFrame(3,3);
t439 = sin(t453);
t443 = cos(t453);
t390 = t443 * t418 + t439 * t421;
t370 = pkin(2) * (t443 * t429 + t439 * t432) + t390 * pkin(3);
t500 = t370 * t506;
t454 = legFrame(2,3);
t440 = sin(t454);
t444 = cos(t454);
t392 = t444 * t419 + t440 * t422;
t371 = pkin(2) * (t444 * t430 + t440 * t433) + t392 * pkin(3);
t499 = t371 * t505;
t455 = legFrame(1,3);
t441 = sin(t455);
t445 = cos(t455);
t394 = t445 * t420 + t441 * t423;
t372 = pkin(2) * (t445 * t431 + t441 * t434) + t394 * pkin(3);
t498 = t372 * t504;
t391 = -t418 * t439 + t443 * t421;
t373 = -pkin(2) * (t429 * t439 - t443 * t432) + t391 * pkin(3);
t497 = t373 * t506;
t393 = -t419 * t440 + t444 * t422;
t374 = -pkin(2) * (t430 * t440 - t444 * t433) + t393 * pkin(3);
t496 = t374 * t505;
t395 = -t420 * t441 + t445 * t423;
t375 = -pkin(2) * (t431 * t441 - t445 * t434) + t395 * pkin(3);
t495 = t375 * t504;
t494 = t507 * t388;
t493 = t507 * t389;
t492 = t506 * t390;
t491 = t506 * t391;
t490 = t505 * t392;
t489 = t505 * t393;
t488 = t504 * t394;
t487 = t504 * t395;
t484 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t414 = t484 * m(3) + Icges(3,3);
t469 = 0.1e1 / pkin(3);
t486 = t414 * t469;
t470 = 0.1e1 / pkin(2);
t485 = t469 * t470;
t474 = ((t426 * rSges(3,1) - rSges(3,2) * t427) * t416 + (rSges(3,1) * t427 + t426 * rSges(3,2)) * t417) * t503;
t364 = t474 + t414;
t483 = t364 * t507 * t469;
t473 = ((t429 * rSges(3,1) - rSges(3,2) * t432) * t418 + (rSges(3,1) * t432 + t429 * rSges(3,2)) * t421) * t503;
t365 = t473 + t414;
t482 = t365 * t506 * t469;
t472 = ((t430 * rSges(3,1) - rSges(3,2) * t433) * t419 + (rSges(3,1) * t433 + t430 * rSges(3,2)) * t422) * t503;
t366 = t472 + t414;
t481 = t366 * t505 * t469;
t471 = ((t431 * rSges(3,1) - rSges(3,2) * t434) * t420 + (rSges(3,1) * t434 + t431 * rSges(3,2)) * t423) * t503;
t367 = t471 + t414;
t480 = t367 * t504 * t469;
t479 = t507 * t486;
t478 = t506 * t486;
t477 = t505 * t486;
t476 = t504 * t486;
t475 = Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t484) * m(3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2);
t466 = koppelP(1,1);
t465 = koppelP(2,1);
t464 = koppelP(3,1);
t463 = koppelP(4,1);
t462 = koppelP(1,2);
t461 = koppelP(2,2);
t460 = koppelP(3,2);
t459 = koppelP(4,2);
t458 = rSges(4,1);
t457 = rSges(4,2);
t456 = xP(4);
t447 = cos(t456);
t446 = sin(t456);
t413 = -t446 * t462 + t447 * t466;
t412 = -t446 * t461 + t447 * t465;
t411 = -t446 * t460 + t447 * t464;
t410 = -t446 * t459 + t447 * t463;
t409 = -t446 * t466 - t447 * t462;
t408 = -t446 * t465 - t447 * t461;
t407 = -t446 * t464 - t447 * t460;
t406 = -t446 * t463 - t447 * t459;
t397 = m(4) * (-t446 * t457 + t447 * t458);
t396 = m(4) * (-t446 * t458 - t447 * t457);
t363 = 0.2e1 * t471 + t475;
t362 = 0.2e1 * t472 + t475;
t361 = 0.2e1 * t473 + t475;
t360 = 0.2e1 * t474 + t475;
t359 = (t394 * t413 + t395 * t409) * t470 * t504;
t358 = (t392 * t412 + t393 * t408) * t470 * t505;
t357 = (t390 * t411 + t391 * t407) * t470 * t506;
t356 = (t388 * t410 + t389 * t406) * t470 * t507;
t355 = (t372 * t413 + t375 * t409) * t504 * t485;
t354 = (t371 * t412 + t374 * t408) * t505 * t485;
t353 = (t370 * t411 + t373 * t407) * t506 * t485;
t352 = (t368 * t410 + t369 * t406) * t507 * t485;
t351 = (t367 * t487 - t375 * t476) * t470;
t350 = (t367 * t488 - t372 * t476) * t470;
t349 = (t366 * t489 - t374 * t477) * t470;
t348 = (t366 * t490 - t371 * t477) * t470;
t347 = (t365 * t491 - t373 * t478) * t470;
t346 = (t365 * t492 - t370 * t478) * t470;
t345 = (t364 * t493 - t369 * t479) * t470;
t344 = (t364 * t494 - t368 * t479) * t470;
t343 = (t363 * t487 - t375 * t480) * t470;
t342 = (t363 * t488 - t372 * t480) * t470;
t341 = (t362 * t489 - t374 * t481) * t470;
t340 = (t362 * t490 - t371 * t481) * t470;
t339 = (t361 * t491 - t373 * t482) * t470;
t338 = (t361 * t492 - t370 * t482) * t470;
t337 = (t360 * t493 - t369 * t483) * t470;
t336 = (t360 * t494 - t368 * t483) * t470;
t335 = -t355 * t414 + t359 * t367;
t334 = -t354 * t414 + t358 * t366;
t333 = -t353 * t414 + t357 * t365;
t332 = -t352 * t414 + t356 * t364;
t331 = -t355 * t367 + t359 * t363;
t330 = -t354 * t366 + t358 * t362;
t329 = -t353 * t365 + t357 * t361;
t328 = -t352 * t364 + t356 * t360;
t1 = [m(4) + (t337 * t493 + t339 * t491 + t341 * t489 + t343 * t487 + (-t345 * t501 - t347 * t497 - t349 * t496 - t351 * t495) * t469) * t470, (t337 * t494 + t339 * t492 + t341 * t490 + t343 * t488 + (-t345 * t502 - t347 * t500 - t349 * t499 - t351 * t498) * t469) * t470, 0, t337 * t356 + t339 * t357 + t341 * t358 + t343 * t359 - t345 * t352 - t347 * t353 - t349 * t354 - t351 * t355 + t396; (t336 * t493 + t338 * t491 + t340 * t489 + t342 * t487 + (-t344 * t501 - t346 * t497 - t348 * t496 - t350 * t495) * t469) * t470, m(4) + (t336 * t494 + t338 * t492 + t340 * t490 + t342 * t488 + (-t344 * t502 - t346 * t500 - t348 * t499 - t350 * t498) * t469) * t470, 0, t336 * t356 + t338 * t357 + t340 * t358 + t342 * t359 - t344 * t352 - t346 * t353 - t348 * t354 - t350 * t355 + t397; 0, 0, (4 * m(1)) + 0.4e1 * m(2) + 0.4e1 * m(3) + m(4), 0; t396 + (t328 * t493 + t329 * t491 + t330 * t489 + t331 * t487 + (-t332 * t501 - t333 * t497 - t334 * t496 - t335 * t495) * t469) * t470, t397 + (t328 * t494 + t329 * t492 + t330 * t490 + t331 * t488 + (-t332 * t502 - t333 * t500 - t334 * t499 - t335 * t498) * t469) * t470, 0, t331 * t359 - t335 * t355 + t330 * t358 - t334 * t354 + t329 * t357 - t333 * t353 + t328 * t356 - t332 * t352 + Icges(4,3) + m(4) * (t457 ^ 2 + t458 ^ 2);];
MX  = t1;
