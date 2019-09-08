% Calculate inertia matrix for parallel robot
% P3PRP1G1P1A0
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
% Datum: 2019-05-03 14:42
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRP1G1P1A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1G1P1A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:41:33
% EndTime: 2019-05-03 14:41:34
% DurationCPUTime: 1.20s
% Computational Cost: add. (6088->252), mult. (11247->374), div. (432->3), fcn. (5066->14), ass. (0->157)
t434 = 2 * pkin(2);
t433 = 2 * rSges(3,3);
t392 = pkin(2) ^ 2;
t432 = 1 + t392;
t381 = (qJ(3,3) ^ 2);
t351 = -t381 + t432;
t373 = cos(qJ(2,3));
t364 = t373 ^ 2;
t370 = sin(qJ(2,3));
t415 = t373 * qJ(3,3);
t335 = 0.1e1 / (t370 * t415 * t434 + t351 * t364 - t381 - t432);
t431 = m(3) * t335;
t382 = (qJ(3,2) ^ 2);
t352 = -t382 + t432;
t374 = cos(qJ(2,2));
t365 = t374 ^ 2;
t371 = sin(qJ(2,2));
t414 = t374 * qJ(3,2);
t336 = 0.1e1 / (t371 * t414 * t434 + t352 * t365 - t382 - t432);
t430 = m(3) * t336;
t383 = (qJ(3,1) ^ 2);
t353 = -t383 + t432;
t375 = cos(qJ(2,1));
t366 = t375 ^ 2;
t372 = sin(qJ(2,1));
t413 = t375 * qJ(3,1);
t337 = 0.1e1 / (t372 * t413 * t434 + t353 * t366 - t383 - t432);
t429 = m(3) * t337;
t428 = m(3) * t373;
t427 = m(3) * t374;
t426 = m(3) * t375;
t376 = pkin(2) + rSges(3,1);
t425 = t376 * m(3);
t369 = legFrame(1,3);
t356 = sin(t369);
t424 = qJ(3,1) * t356;
t359 = cos(t369);
t423 = qJ(3,1) * t359;
t368 = legFrame(2,3);
t355 = sin(t368);
t422 = qJ(3,2) * t355;
t358 = cos(t368);
t421 = qJ(3,2) * t358;
t367 = legFrame(3,3);
t354 = sin(t367);
t420 = qJ(3,3) * t354;
t357 = cos(t367);
t419 = qJ(3,3) * t357;
t418 = t370 * t373;
t417 = t371 * t374;
t416 = t372 * t375;
t412 = -0.2e1 * t415;
t411 = -0.2e1 * t414;
t410 = -0.2e1 * t413;
t409 = pkin(2) * t424;
t408 = pkin(2) * t422;
t407 = pkin(2) * t420;
t406 = pkin(2) * t419;
t405 = pkin(2) * t421;
t404 = pkin(2) * t423;
t403 = (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + Icges(3,2) + Icges(2,3);
t402 = rSges(3,3) ^ 2 + t392 + (t434 + rSges(3,1)) * rSges(3,1);
t401 = -t356 * t353 + 0.2e1 * t404;
t400 = -t383 * t356 - t404;
t399 = t359 * t383 - t409;
t398 = -t355 * t352 + 0.2e1 * t405;
t397 = -t382 * t355 - t405;
t396 = t358 * t382 - t408;
t395 = -t354 * t351 + 0.2e1 * t406;
t394 = -t381 * t354 - t406;
t393 = t357 * t381 - t407;
t389 = koppelP(1,1);
t388 = koppelP(2,1);
t387 = koppelP(3,1);
t386 = koppelP(1,2);
t385 = koppelP(2,2);
t384 = koppelP(3,2);
t380 = rSges(4,1);
t379 = rSges(4,2);
t378 = xP(3);
t377 = m(2) * rSges(2,2);
t362 = m(1) + m(2) + m(3);
t361 = cos(t378);
t360 = sin(t378);
t349 = m(2) * rSges(2,1) + t425;
t348 = -t360 * t386 + t361 * t389;
t347 = -t360 * t385 + t361 * t388;
t346 = -t360 * t384 + t361 * t387;
t345 = -t360 * t389 - t361 * t386;
t344 = -t360 * t388 - t361 * t385;
t343 = -t360 * t387 - t361 * t384;
t342 = m(4) * (-t360 * t379 + t361 * t380);
t341 = m(4) * (-t360 * t380 - t361 * t379);
t340 = t353 * t359 + 0.2e1 * t409;
t339 = t352 * t358 + 0.2e1 * t408;
t338 = t351 * t357 + 0.2e1 * t407;
t334 = t349 * t375 - t372 * ((-rSges(3,3) - qJ(3,1)) * m(3) + t377);
t333 = t349 * t374 - t371 * ((-rSges(3,3) - qJ(3,2)) * m(3) + t377);
t332 = t349 * t373 - t370 * ((-rSges(3,3) - qJ(3,3)) * m(3) + t377);
t331 = (qJ(3,1) * t433 + t383 + t402) * m(3) + t403;
t330 = (qJ(3,2) * t433 + t382 + t402) * m(3) + t403;
t329 = (qJ(3,3) * t433 + t381 + t402) * m(3) + t403;
t328 = t359 * t410 + t372 * (pkin(2) * t359 + t424);
t327 = t356 * t410 + t372 * (pkin(2) * t356 - t423);
t326 = t358 * t411 + t371 * (pkin(2) * t358 + t422);
t325 = t355 * t411 + t371 * (pkin(2) * t355 - t421);
t324 = t357 * t412 + t370 * (pkin(2) * t357 + t420);
t323 = t354 * t412 + t370 * (pkin(2) * t354 - t419);
t322 = t399 * t375 - t372 * (t356 - t400);
t321 = t396 * t374 - t371 * (t355 - t397);
t320 = t393 * t373 - t370 * (t354 - t394);
t319 = t400 * t375 + t372 * (-t359 - t399);
t318 = t397 * t374 + t371 * (-t358 - t396);
t317 = t394 * t373 + t370 * (-t357 - t393);
t316 = -t340 * t416 + t392 * t356 + t401 * t366 + t356 - t404;
t315 = t340 * t366 - t359 * t392 + t401 * t416 - t359 - t409;
t314 = -t339 * t417 + t392 * t355 + t398 * t365 + t355 - t405;
t313 = t339 * t365 - t358 * t392 + t398 * t417 - t358 - t408;
t312 = -t338 * t418 + t392 * t354 + t395 * t364 + t354 - t406;
t311 = t338 * t364 - t357 * t392 + t395 * t418 - t357 - t407;
t310 = (t327 * t348 + t328 * t345) * t337;
t309 = (t325 * t347 + t326 * t344) * t336;
t308 = (t323 * t346 + t324 * t343) * t335;
t307 = (t319 * t345 + t322 * t348) * t337;
t306 = (t318 * t344 + t321 * t347) * t336;
t305 = (t317 * t343 + t320 * t346) * t335;
t304 = (t315 * t348 + t316 * t345) * t337;
t303 = (t313 * t347 + t314 * t344) * t336;
t302 = (t311 * t346 + t312 * t343) * t335;
t301 = (-t316 * t375 - t328 * t376 + t319) * t429;
t300 = (-t315 * t375 - t327 * t376 + t322) * t429;
t299 = (-t314 * t374 - t326 * t376 + t318) * t430;
t298 = (-t313 * t374 - t325 * t376 + t321) * t430;
t297 = (-t312 * t373 - t324 * t376 + t317) * t431;
t296 = (-t311 * t373 - t323 * t376 + t320) * t431;
t295 = (t316 * t362 - t319 * t426 + t328 * t334) * t337;
t294 = (t315 * t362 - t322 * t426 + t327 * t334) * t337;
t293 = (t314 * t362 - t318 * t427 + t326 * t333) * t336;
t292 = (t313 * t362 - t321 * t427 + t325 * t333) * t336;
t291 = (t312 * t362 - t317 * t428 + t324 * t332) * t335;
t290 = (t311 * t362 - t320 * t428 + t323 * t332) * t335;
t289 = (t316 * t334 - t319 * t425 + t328 * t331) * t337;
t288 = (t315 * t334 - t322 * t425 + t327 * t331) * t337;
t287 = (t314 * t333 - t318 * t425 + t326 * t330) * t336;
t286 = (t313 * t333 - t321 * t425 + t325 * t330) * t336;
t285 = (t312 * t332 - t317 * t425 + t324 * t329) * t335;
t284 = (t311 * t332 - t320 * t425 + t323 * t329) * t335;
t283 = (-t304 * t375 - t310 * t376 + t307) * m(3);
t282 = (-t303 * t374 - t309 * t376 + t306) * m(3);
t281 = (-t302 * t373 - t308 * t376 + t305) * m(3);
t280 = t304 * t362 - t307 * t426 + t310 * t334;
t279 = t303 * t362 - t306 * t427 + t309 * t333;
t278 = t302 * t362 - t305 * t428 + t308 * t332;
t277 = t304 * t334 - t307 * t425 + t310 * t331;
t276 = t303 * t333 - t306 * t425 + t309 * t330;
t275 = t302 * t332 - t305 * t425 + t308 * t329;
t1 = [m(4) + (t289 * t328 + t295 * t316 + t301 * t319) * t337 + (t287 * t326 + t293 * t314 + t299 * t318) * t336 + (t285 * t324 + t291 * t312 + t297 * t317) * t335, (t289 * t327 + t295 * t315 + t301 * t322) * t337 + (t287 * t325 + t293 * t313 + t299 * t321) * t336 + (t285 * t323 + t291 * t311 + t297 * t320) * t335, t285 * t308 + t287 * t309 + t289 * t310 + t291 * t302 + t293 * t303 + t295 * t304 + t297 * t305 + t299 * t306 + t301 * t307 + t341; (t288 * t328 + t294 * t316 + t300 * t319) * t337 + (t286 * t326 + t292 * t314 + t298 * t318) * t336 + (t284 * t324 + t290 * t312 + t296 * t317) * t335, m(4) + (t288 * t327 + t294 * t315 + t300 * t322) * t337 + (t286 * t325 + t292 * t313 + t298 * t321) * t336 + (t284 * t323 + t290 * t311 + t296 * t320) * t335, t284 * t308 + t286 * t309 + t288 * t310 + t290 * t302 + t292 * t303 + t294 * t304 + t296 * t305 + t298 * t306 + t300 * t307 + t342; t341 + (t277 * t328 + t280 * t316 + t283 * t319) * t337 + (t276 * t326 + t279 * t314 + t282 * t318) * t336 + (t275 * t324 + t278 * t312 + t281 * t317) * t335, t342 + (t277 * t327 + t280 * t315 + t283 * t322) * t337 + (t276 * t325 + t279 * t313 + t282 * t321) * t336 + (t275 * t323 + t278 * t311 + t281 * t320) * t335, t280 * t304 + t277 * t310 + t283 * t307 + t279 * t303 + t276 * t309 + t282 * t306 + t278 * t302 + t275 * t308 + t281 * t305 + Icges(4,3) + m(4) * (t379 ^ 2 + t380 ^ 2);];
MX  = t1;
