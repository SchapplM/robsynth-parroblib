% Calculate inertia matrix for parallel robot
% P3RPP1G1P1A0
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
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
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
% MX [3x3]
%   inertia matrix in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPP1G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:00
% EndTime: 2019-05-03 14:52:01
% DurationCPUTime: 0.87s
% Computational Cost: add. (5377->236), mult. (7134->341), div. (360->3), fcn. (2528->14), ass. (0->169)
t473 = 2 * pkin(1);
t472 = 2 * qJ(2,1);
t471 = 2 * qJ(3,1);
t470 = 2 * qJ(2,2);
t469 = 2 * qJ(3,2);
t468 = 2 * qJ(2,3);
t467 = 2 * qJ(3,3);
t443 = pkin(1) ^ 2;
t466 = 1 + t443;
t427 = qJ(3,3) ^ 2;
t385 = pkin(1) * t467 + t427 + t466;
t428 = qJ(2,3) ^ 2;
t382 = 1 / (t428 + t385);
t465 = m(3) * t382;
t429 = qJ(3,2) ^ 2;
t386 = pkin(1) * t469 + t429 + t466;
t430 = qJ(2,2) ^ 2;
t383 = 1 / (t430 + t386);
t464 = m(3) * t383;
t431 = qJ(3,1) ^ 2;
t387 = pkin(1) * t471 + t431 + t466;
t432 = qJ(2,1) ^ 2;
t384 = 1 / (t432 + t387);
t463 = m(3) * t384;
t407 = rSges(3,2) + qJ(2,3);
t462 = m(3) * t407;
t408 = rSges(3,2) + qJ(2,2);
t461 = m(3) * t408;
t409 = rSges(3,2) + qJ(2,1);
t460 = m(3) * t409;
t459 = (pkin(1) - rSges(2,2)) * m(2);
t458 = rSges(3,3) + qJ(3,1);
t457 = rSges(3,3) + qJ(3,2);
t456 = rSges(3,3) + qJ(3,3);
t413 = pkin(1) + qJ(3,3);
t416 = sin(qJ(1,3));
t455 = t416 * t413;
t414 = pkin(1) + qJ(3,2);
t417 = sin(qJ(1,2));
t454 = t417 * t414;
t415 = pkin(1) + qJ(3,1);
t418 = sin(qJ(1,1));
t453 = t418 * t415;
t419 = cos(qJ(1,3));
t452 = t419 * qJ(2,3);
t420 = cos(qJ(1,2));
t451 = t420 * qJ(2,2);
t421 = cos(qJ(1,1));
t450 = t421 * qJ(2,1);
t449 = t415 * t450;
t448 = t414 * t451;
t447 = t413 * t452;
t446 = rSges(3,2) ^ 2 + rSges(3,3) ^ 2 + t443;
t445 = ((rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,1) + Icges(3,1) + Icges(1,3));
t444 = rSges(2,3) ^ 2 + t443 + (-2 * pkin(1) + rSges(2,2)) * rSges(2,2);
t438 = koppelP(1,1);
t437 = koppelP(2,1);
t436 = koppelP(3,1);
t435 = koppelP(1,2);
t434 = koppelP(2,2);
t433 = koppelP(3,2);
t426 = rSges(4,1);
t425 = rSges(4,2);
t424 = xP(3);
t423 = m(2) + m(3);
t412 = legFrame(1,3);
t411 = legFrame(2,3);
t410 = legFrame(3,3);
t405 = 1 + t432;
t404 = 1 + t430;
t403 = 1 + t428;
t399 = cos(t424);
t398 = sin(t424);
t397 = cos(t412);
t396 = cos(t411);
t395 = cos(t410);
t394 = sin(t412);
t393 = sin(t411);
t392 = sin(t410);
t390 = qJ(2,1) * t453;
t389 = qJ(2,2) * t454;
t388 = qJ(2,3) * t455;
t381 = (-pkin(1) - t458) * m(3) - t459;
t380 = (-pkin(1) - t457) * m(3) - t459;
t379 = (-pkin(1) - t456) * m(3) - t459;
t378 = t418 * qJ(2,1) + t415 * t421;
t377 = t417 * qJ(2,2) + t414 * t420;
t376 = t416 * qJ(2,3) + t413 * t419;
t375 = t450 - t453;
t374 = t451 - t454;
t373 = t452 - t455;
t372 = -t398 * t435 + t399 * t438;
t371 = -t398 * t434 + t399 * t437;
t370 = -t398 * t433 + t399 * t436;
t369 = -t398 * t438 - t399 * t435;
t368 = -t398 * t437 - t399 * t434;
t367 = -t398 * t436 - t399 * t433;
t366 = -t405 * t421 + t390;
t365 = -t404 * t420 + t389;
t364 = -t403 * t419 + t388;
t363 = t418 * t405 + t449;
t362 = t417 * t404 + t448;
t361 = t416 * t403 + t447;
t360 = m(4) * (-t398 * t425 + t399 * t426);
t359 = m(4) * (-t398 * t426 - t399 * t425);
t358 = t387 * t421 + t390;
t357 = t386 * t420 + t389;
t356 = t385 * t419 + t388;
t355 = t418 * t387 - t449;
t354 = t417 * t386 - t448;
t353 = t416 * t385 - t447;
t352 = t394 * t375 + t378 * t397;
t351 = t393 * t374 + t377 * t396;
t350 = t392 * t373 + t376 * t395;
t349 = t375 * t397 - t394 * t378;
t348 = t374 * t396 - t393 * t377;
t347 = t373 * t395 - t392 * t376;
t346 = t394 * t363 + t366 * t397;
t345 = t393 * t362 + t365 * t396;
t344 = t392 * t361 + t364 * t395;
t343 = t363 * t397 - t394 * t366;
t342 = t362 * t396 - t393 * t365;
t341 = t361 * t395 - t392 * t364;
t340 = (rSges(3,2) * t472 + rSges(3,3) * t471 + t458 * t473 + t431 + t432 + t446) * m(3) + (rSges(2,3) * t472 + t432 + t444) * m(2) + t445;
t339 = (rSges(3,2) * t470 + rSges(3,3) * t469 + t457 * t473 + t429 + t430 + t446) * m(3) + (rSges(2,3) * t470 + t430 + t444) * m(2) + t445;
t338 = (rSges(3,2) * t468 + rSges(3,3) * t467 + t456 * t473 + t427 + t428 + t446) * m(3) + (rSges(2,3) * t468 + t428 + t444) * m(2) + t445;
t337 = -t394 * t355 + t358 * t397;
t336 = -t393 * t354 + t357 * t396;
t335 = -t392 * t353 + t356 * t395;
t334 = t355 * t397 + t394 * t358;
t333 = t354 * t396 + t393 * t357;
t332 = t353 * t395 + t392 * t356;
t331 = (t346 * t423 + t352 * t381) * t384;
t330 = (t345 * t423 + t351 * t380) * t383;
t329 = (t344 * t423 + t350 * t379) * t382;
t328 = (t343 * t423 + t349 * t381) * t384;
t327 = (t342 * t423 + t348 * t380) * t383;
t326 = (t341 * t423 + t347 * t379) * t382;
t325 = (t352 * t409 + t334) * t463;
t324 = (t351 * t408 + t333) * t464;
t323 = (t350 * t407 + t332) * t465;
t322 = (t349 * t409 + t337) * t463;
t321 = (t348 * t408 + t336) * t464;
t320 = (t347 * t407 + t335) * t465;
t319 = (t349 * t369 + t352 * t372) * t384;
t318 = (t348 * t368 + t351 * t371) * t383;
t317 = (t347 * t367 + t350 * t370) * t382;
t316 = (t343 * t369 + t346 * t372) * t384;
t315 = (t342 * t368 + t345 * t371) * t383;
t314 = (t341 * t367 + t344 * t370) * t382;
t313 = (t334 * t372 + t337 * t369) * t384;
t312 = (t333 * t371 + t336 * t368) * t383;
t311 = (t332 * t370 + t335 * t367) * t382;
t310 = (t334 * t460 + t340 * t352 + t346 * t381) * t384;
t309 = (t333 * t461 + t339 * t351 + t345 * t380) * t383;
t308 = (t332 * t462 + t338 * t350 + t344 * t379) * t382;
t307 = (t337 * t460 + t340 * t349 + t343 * t381) * t384;
t306 = (t336 * t461 + t339 * t348 + t342 * t380) * t383;
t305 = (t335 * t462 + t338 * t347 + t341 * t379) * t382;
t304 = t316 * t423 + t319 * t381;
t303 = t315 * t423 + t318 * t380;
t302 = t314 * t423 + t317 * t379;
t301 = (t319 * t409 + t313) * m(3);
t300 = (t318 * t408 + t312) * m(3);
t299 = (t317 * t407 + t311) * m(3);
t298 = t313 * t460 + t316 * t381 + t319 * t340;
t297 = t312 * t461 + t315 * t380 + t318 * t339;
t296 = t311 * t462 + t314 * t379 + t317 * t338;
t1 = [m(4) + (t307 * t349 + t322 * t337 + t328 * t343) * t384 + (t306 * t348 + t321 * t336 + t327 * t342) * t383 + (t305 * t347 + t320 * t335 + t326 * t341) * t382, (t307 * t352 + t322 * t334 + t328 * t346) * t384 + (t306 * t351 + t321 * t333 + t327 * t345) * t383 + (t305 * t350 + t320 * t332 + t326 * t344) * t382, t305 * t317 + t306 * t318 + t307 * t319 + t320 * t311 + t321 * t312 + t322 * t313 + t326 * t314 + t327 * t315 + t328 * t316 + t359; (t310 * t349 + t325 * t337 + t331 * t343) * t384 + (t309 * t348 + t324 * t336 + t330 * t342) * t383 + (t308 * t347 + t323 * t335 + t329 * t341) * t382, m(4) + (t310 * t352 + t325 * t334 + t331 * t346) * t384 + (t309 * t351 + t324 * t333 + t330 * t345) * t383 + (t308 * t350 + t323 * t332 + t329 * t344) * t382, t308 * t317 + t309 * t318 + t310 * t319 + t323 * t311 + t324 * t312 + t325 * t313 + t329 * t314 + t330 * t315 + t331 * t316 + t360; t359 + (t298 * t349 + t301 * t337 + t304 * t343) * t384 + (t297 * t348 + t300 * t336 + t303 * t342) * t383 + (t296 * t347 + t299 * t335 + t302 * t341) * t382, t360 + (t298 * t352 + t301 * t334 + t304 * t346) * t384 + (t297 * t351 + t300 * t333 + t303 * t345) * t383 + (t296 * t350 + t299 * t332 + t302 * t344) * t382, t298 * t319 + t304 * t316 + t301 * t313 + t297 * t318 + t303 * t315 + t300 * t312 + t296 * t317 + t302 * t314 + t299 * t311 + Icges(4,3) + m(4) * (t425 ^ 2 + t426 ^ 2);];
MX  = t1;
