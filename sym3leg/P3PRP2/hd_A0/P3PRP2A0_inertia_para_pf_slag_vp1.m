% Calculate inertia matrix for parallel robot
% P3PRP2A0
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
%   pkin=[a2,a3,d2]';
% m [4x1]
%   mass of all robot links (including platform)
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
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MX = P3PRP2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:38:50
% EndTime: 2018-12-20 17:38:51
% DurationCPUTime: 1.07s
% Computational Cost: add. (6088->252), mult. (11247->377), div. (432->3), fcn. (5066->14), ass. (0->156)
t437 = 2 * rSges(3,3);
t402 = (pkin(2) ^ 2);
t436 = -t402 - 1;
t391 = (qJ(3,3) ^ 2);
t361 = -t391 - t436;
t383 = cos(qJ(2,3));
t374 = t383 ^ 2;
t380 = sin(qJ(2,3));
t419 = t383 * qJ(3,3);
t413 = 0.2e1 * t419;
t336 = 0.1e1 / (t380 * pkin(2) * t413 + t361 * t374 - t391 + t436);
t435 = m(3) * t336;
t392 = (qJ(3,2) ^ 2);
t362 = -t392 - t436;
t384 = cos(qJ(2,2));
t375 = t384 ^ 2;
t381 = sin(qJ(2,2));
t418 = t384 * qJ(3,2);
t412 = 0.2e1 * t418;
t337 = 0.1e1 / (t381 * pkin(2) * t412 + t362 * t375 - t392 + t436);
t434 = m(3) * t337;
t393 = (qJ(3,1) ^ 2);
t363 = -t393 - t436;
t385 = cos(qJ(2,1));
t376 = t385 ^ 2;
t382 = sin(qJ(2,1));
t417 = t385 * qJ(3,1);
t411 = 0.2e1 * t417;
t338 = 0.1e1 / (t382 * pkin(2) * t411 + t363 * t376 - t393 + t436);
t433 = m(3) * t338;
t432 = m(3) * t383;
t431 = m(3) * t384;
t430 = m(3) * t385;
t386 = pkin(2) + rSges(3,1);
t429 = t386 * m(3);
t379 = legFrame(1,3);
t366 = sin(t379);
t428 = qJ(3,1) * t366;
t369 = cos(t379);
t427 = qJ(3,1) * t369;
t378 = legFrame(2,3);
t365 = sin(t378);
t426 = qJ(3,2) * t365;
t368 = cos(t378);
t425 = qJ(3,2) * t368;
t377 = legFrame(3,3);
t364 = sin(t377);
t424 = qJ(3,3) * t364;
t367 = cos(t377);
t423 = qJ(3,3) * t367;
t422 = t380 * t383;
t421 = t381 * t384;
t420 = t382 * t385;
t354 = pkin(2) * t424;
t416 = t391 * t367 + t354;
t355 = pkin(2) * t426;
t415 = t392 * t368 + t355;
t356 = pkin(2) * t428;
t414 = t393 * t369 + t356;
t410 = pkin(2) * t427;
t409 = pkin(2) * t425;
t408 = pkin(2) * t423;
t407 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,2) + Icges(2,3);
t406 = (rSges(3,3) ^ 2) + t402 + (0.2e1 * pkin(2) + rSges(3,1)) * rSges(3,1);
t405 = t366 * t393 - t410;
t404 = t365 * t392 - t409;
t403 = t364 * t391 - t408;
t399 = koppelP(1,1);
t398 = koppelP(2,1);
t397 = koppelP(3,1);
t396 = koppelP(1,2);
t395 = koppelP(2,2);
t394 = koppelP(3,2);
t390 = rSges(4,1);
t389 = rSges(4,2);
t388 = xP(3);
t387 = m(2) * rSges(2,2);
t372 = m(1) + m(2) + m(3);
t371 = cos(t388);
t370 = sin(t388);
t353 = m(2) * rSges(2,1) + t429;
t352 = -t370 * t396 + t371 * t399;
t351 = -t370 * t395 + t371 * t398;
t350 = -t370 * t394 + t371 * t397;
t349 = -t370 * t399 - t371 * t396;
t348 = -t370 * t398 - t371 * t395;
t347 = -t370 * t397 - t371 * t394;
t346 = m(4) * (-t370 * t389 + t371 * t390);
t345 = m(4) * (-t370 * t390 - t371 * t389);
t344 = t363 * t366 + 0.2e1 * t410;
t343 = t362 * t365 + 0.2e1 * t409;
t342 = t361 * t364 + 0.2e1 * t408;
t341 = t369 * t363 - 0.2e1 * t356;
t340 = t368 * t362 - 0.2e1 * t355;
t339 = t367 * t361 - 0.2e1 * t354;
t335 = t353 * t385 - ((-rSges(3,3) - qJ(3,1)) * m(3) + t387) * t382;
t334 = t353 * t384 - ((-rSges(3,3) - qJ(3,2)) * m(3) + t387) * t381;
t333 = t353 * t383 - ((-rSges(3,3) - qJ(3,3)) * m(3) + t387) * t380;
t332 = (qJ(3,1) * t437 + t393 + t406) * m(3) + t407;
t331 = (qJ(3,2) * t437 + t392 + t406) * m(3) + t407;
t330 = (qJ(3,3) * t437 + t391 + t406) * m(3) + t407;
t329 = -0.2e1 * t369 * t417 + t382 * (pkin(2) * t369 - t428);
t328 = t366 * t411 - t382 * (pkin(2) * t366 + t427);
t327 = -0.2e1 * t368 * t418 + t381 * (pkin(2) * t368 - t426);
t326 = t365 * t412 - t381 * (pkin(2) * t365 + t425);
t325 = -0.2e1 * t367 * t419 + t380 * (pkin(2) * t367 - t424);
t324 = t364 * t413 - t380 * (pkin(2) * t364 + t423);
t323 = t405 * t385 - t382 * (t369 + t414);
t322 = t404 * t384 - t381 * (t368 + t415);
t321 = t403 * t383 - t380 * (t367 + t416);
t320 = t414 * t385 - t382 * (-t366 - t405);
t319 = t415 * t384 - t381 * (-t365 - t404);
t318 = t416 * t383 - t380 * (-t364 - t403);
t317 = -t341 * t420 + t344 * t376 - t366 * t402 - t366 - t410;
t316 = -t340 * t421 + t343 * t375 - t365 * t402 - t365 - t409;
t315 = -t339 * t422 + t342 * t374 - t364 * t402 - t364 - t408;
t314 = t341 * t376 + t344 * t420 - t402 * t369 + t356 - t369;
t313 = t340 * t375 + t343 * t421 - t402 * t368 + t355 - t368;
t312 = t339 * t374 + t342 * t422 - t402 * t367 + t354 - t367;
t311 = (t328 * t349 + t329 * t352) * t338;
t310 = (t326 * t348 + t327 * t351) * t337;
t309 = (t324 * t347 + t325 * t350) * t336;
t308 = (t320 * t349 + t323 * t352) * t338;
t307 = (t319 * t348 + t322 * t351) * t337;
t306 = (t318 * t347 + t321 * t350) * t336;
t305 = (t314 * t349 + t317 * t352) * t338;
t304 = (t313 * t348 + t316 * t351) * t337;
t303 = (t312 * t347 + t315 * t350) * t336;
t302 = (-t317 * t385 - t329 * t386 + t323) * t433;
t301 = (-t316 * t384 - t327 * t386 + t322) * t434;
t300 = (-t315 * t383 - t325 * t386 + t321) * t435;
t299 = (-t314 * t385 - t328 * t386 + t320) * t433;
t298 = (-t313 * t384 - t326 * t386 + t319) * t434;
t297 = (-t312 * t383 - t324 * t386 + t318) * t435;
t296 = (t317 * t372 - t323 * t430 + t329 * t335) * t338;
t295 = (t316 * t372 - t322 * t431 + t327 * t334) * t337;
t294 = (t315 * t372 - t321 * t432 + t325 * t333) * t336;
t293 = (t314 * t372 - t320 * t430 + t328 * t335) * t338;
t292 = (t313 * t372 - t319 * t431 + t326 * t334) * t337;
t291 = (t312 * t372 - t318 * t432 + t324 * t333) * t336;
t290 = (t317 * t335 - t323 * t429 + t329 * t332) * t338;
t289 = (t316 * t334 - t322 * t429 + t327 * t331) * t337;
t288 = (t315 * t333 - t321 * t429 + t325 * t330) * t336;
t287 = (t314 * t335 - t320 * t429 + t328 * t332) * t338;
t286 = (t313 * t334 - t319 * t429 + t326 * t331) * t337;
t285 = (t312 * t333 - t318 * t429 + t324 * t330) * t336;
t284 = (-t305 * t385 - t311 * t386 + t308) * m(3);
t283 = (-t304 * t384 - t310 * t386 + t307) * m(3);
t282 = (-t303 * t383 - t309 * t386 + t306) * m(3);
t281 = t305 * t372 - t308 * t430 + t311 * t335;
t280 = t304 * t372 - t307 * t431 + t310 * t334;
t279 = t303 * t372 - t306 * t432 + t309 * t333;
t278 = t305 * t335 - t308 * t429 + t311 * t332;
t277 = t304 * t334 - t307 * t429 + t310 * t331;
t276 = t303 * t333 - t306 * t429 + t309 * t330;
t1 = [m(4) + (t287 * t328 + t293 * t314 + t299 * t320) * t338 + (t286 * t326 + t292 * t313 + t298 * t319) * t337 + (t285 * t324 + t291 * t312 + t297 * t318) * t336 (t287 * t329 + t293 * t317 + t299 * t323) * t338 + (t286 * t327 + t292 * t316 + t298 * t322) * t337 + (t285 * t325 + t291 * t315 + t297 * t321) * t336, t285 * t309 + t286 * t310 + t287 * t311 + t291 * t303 + t292 * t304 + t293 * t305 + t297 * t306 + t298 * t307 + t299 * t308 + t345; (t290 * t328 + t296 * t314 + t302 * t320) * t338 + (t289 * t326 + t295 * t313 + t301 * t319) * t337 + (t288 * t324 + t294 * t312 + t300 * t318) * t336, m(4) + (t290 * t329 + t296 * t317 + t302 * t323) * t338 + (t289 * t327 + t295 * t316 + t301 * t322) * t337 + (t288 * t325 + t294 * t315 + t300 * t321) * t336, t288 * t309 + t289 * t310 + t290 * t311 + t294 * t303 + t295 * t304 + t296 * t305 + t300 * t306 + t301 * t307 + t302 * t308 + t346; t345 + (t278 * t328 + t281 * t314 + t284 * t320) * t338 + (t277 * t326 + t280 * t313 + t283 * t319) * t337 + (t276 * t324 + t279 * t312 + t282 * t318) * t336, t346 + (t278 * t329 + t281 * t317 + t284 * t323) * t338 + (t277 * t327 + t280 * t316 + t283 * t322) * t337 + (t276 * t325 + t279 * t315 + t282 * t321) * t336, t281 * t305 + t278 * t311 + t284 * t308 + t280 * t304 + t277 * t310 + t283 * t307 + t279 * t303 + t276 * t309 + t282 * t306 + Icges(4,3) + m(4) * (t389 ^ 2 + t390 ^ 2);];
MX  = t1;
