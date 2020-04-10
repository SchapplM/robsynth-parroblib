% Calculate minimal parameter regressor of inertia matrix for parallel robot
% P3PRRR2G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
% koppelP [3x3]
%   coordinates of the platform coupling joints
%   (joints that link the end of legs with platform)
%   in platform coordinates
%   rows: number of leg
%   columns: x-, y-, z-coordinates

% Output:
% tau_reg [3*3x8]
%   minimal parameter regressor of inertia matrix for parallel robot
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:21
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR2G2P3A0_inertia_para_pf_regmin(xP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G2P3A0_inertia_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G2P3A0_inertia_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G2P3A0_inertia_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G2P3A0_inertia_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G2P3A0_inertia_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_MMreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:21:39
% EndTime: 2020-03-09 21:21:40
% DurationCPUTime: 1.02s
% Computational Cost: add. (929->175), mult. (1716->374), div. (918->9), fcn. (2226->21), ass. (0->176)
t356 = cos(qJ(3,3));
t357 = cos(qJ(2,3));
t350 = sin(qJ(3,3));
t351 = sin(qJ(2,3));
t427 = t351 * t350;
t311 = (pkin(2) * t356 + pkin(1)) * t357 - pkin(2) * t427;
t320 = t357 * t356 - t427;
t466 = t311 * t320;
t362 = 0.1e1 / pkin(2);
t465 = t311 * t362;
t358 = cos(qJ(3,2));
t359 = cos(qJ(2,2));
t352 = sin(qJ(3,2));
t353 = sin(qJ(2,2));
t426 = t353 * t352;
t312 = (pkin(2) * t358 + pkin(1)) * t359 - pkin(2) * t426;
t321 = t359 * t358 - t426;
t464 = t312 * t321;
t463 = t312 * t362;
t360 = cos(qJ(3,1));
t361 = cos(qJ(2,1));
t354 = sin(qJ(3,1));
t355 = sin(qJ(2,1));
t425 = t355 * t354;
t313 = (pkin(2) * t360 + pkin(1)) * t361 - pkin(2) * t425;
t322 = t361 * t360 - t425;
t462 = t313 * t322;
t461 = t313 * t362;
t342 = 0.1e1 / t350 ^ 2;
t460 = t320 ^ 2 * t342;
t344 = 0.1e1 / t352 ^ 2;
t459 = t321 ^ 2 * t344;
t346 = 0.1e1 / t354 ^ 2;
t458 = t322 ^ 2 * t346;
t347 = legFrame(3,2);
t335 = sin(t347);
t457 = t320 * t335;
t338 = cos(t347);
t456 = t320 * t338;
t348 = legFrame(2,2);
t336 = sin(t348);
t455 = t321 * t336;
t339 = cos(t348);
t454 = t321 * t339;
t349 = legFrame(1,2);
t337 = sin(t349);
t453 = t322 * t337;
t340 = cos(t349);
t452 = t322 * t340;
t332 = sin(qJ(2,3) + qJ(3,3));
t323 = t351 * pkin(1) + pkin(2) * t332;
t341 = 0.1e1 / t350;
t451 = t323 * t341;
t333 = sin(qJ(2,2) + qJ(3,2));
t324 = t353 * pkin(1) + pkin(2) * t333;
t343 = 0.1e1 / t352;
t450 = t324 * t343;
t334 = sin(qJ(2,1) + qJ(3,1));
t325 = t355 * pkin(1) + pkin(2) * t334;
t345 = 0.1e1 / t354;
t449 = t325 * t345;
t448 = t332 * t341;
t447 = t333 * t343;
t446 = t334 * t345;
t445 = t335 * t341;
t444 = t336 * t343;
t443 = t337 * t345;
t442 = t338 * t335;
t441 = t338 * t341;
t440 = t339 * t336;
t439 = t339 * t343;
t438 = t340 * t337;
t437 = t340 * t345;
t436 = t341 * t356;
t363 = 0.1e1 / pkin(1);
t435 = t341 * t363;
t434 = t342 * t356;
t433 = t343 * t358;
t432 = t343 * t363;
t431 = t344 * t358;
t430 = t345 * t360;
t429 = t345 * t363;
t428 = t346 * t360;
t424 = t362 * t363;
t423 = t341 * t466;
t422 = t311 * t445;
t421 = t311 * t441;
t420 = t343 * t464;
t419 = t312 * t444;
t418 = t312 * t439;
t417 = t345 * t462;
t416 = t313 * t443;
t415 = t313 * t437;
t414 = t320 * t332 * t342;
t413 = t320 * t445;
t412 = t320 * t441;
t411 = t320 * t436;
t410 = t320 * t434;
t409 = t321 * t333 * t344;
t408 = t321 * t444;
t407 = t321 * t439;
t406 = t321 * t433;
t405 = t321 * t431;
t404 = t322 * t334 * t346;
t403 = t322 * t443;
t402 = t322 * t437;
t401 = t322 * t430;
t400 = t322 * t428;
t399 = t332 * t436;
t398 = t332 * t434;
t397 = t333 * t433;
t396 = t333 * t431;
t395 = t334 * t430;
t394 = t334 * t428;
t393 = t341 * t424;
t392 = t343 * t424;
t391 = t345 * t424;
t390 = t332 * t435;
t389 = t333 * t432;
t388 = t334 * t429;
t387 = t391 * t438 * t462 + t392 * t440 * t464 + t393 * t442 * t466;
t386 = t311 * t410;
t385 = t311 * t398;
t384 = t312 * t405;
t383 = t312 * t396;
t382 = t313 * t400;
t381 = t313 * t394;
t380 = t335 * t411;
t379 = t338 * t411;
t378 = t338 * t410;
t377 = t336 * t406;
t376 = t339 * t406;
t375 = t339 * t405;
t374 = t337 * t401;
t373 = t340 * t401;
t372 = t340 * t400;
t371 = (t320 - t465) * t435;
t370 = (0.2e1 * t320 - t465) * t435;
t369 = (t321 - t463) * t432;
t368 = (0.2e1 * t321 - t463) * t432;
t367 = (t322 - t461) * t429;
t366 = (0.2e1 * t322 - t461) * t429;
t365 = (-t311 * t335 * t378 - t312 * t336 * t375 - t313 * t337 * t372) * t424;
t364 = 0.1e1 / pkin(1) ^ 2;
t331 = t340 ^ 2;
t330 = t339 ^ 2;
t329 = t338 ^ 2;
t328 = t337 ^ 2;
t327 = t336 ^ 2;
t326 = t335 ^ 2;
t316 = t325 * t391;
t315 = t324 * t392;
t314 = t323 * t393;
t310 = t316 - t388;
t309 = t315 - t389;
t308 = t314 - t390;
t307 = t316 - 0.2e1 * t388;
t306 = t315 - 0.2e1 * t389;
t305 = t314 - 0.2e1 * t390;
t304 = -t438 - t440 - t442;
t300 = t340 * t367;
t299 = t339 * t369;
t298 = t338 * t371;
t297 = t337 * t367;
t296 = t336 * t369;
t295 = t335 * t371;
t294 = t340 * t366;
t293 = t339 * t368;
t292 = t338 * t370;
t291 = t337 * t366;
t290 = t336 * t368;
t289 = t335 * t370;
t288 = (t438 * t458 + t440 * t459 + t442 * t460) * t364;
t287 = (-t338 * t414 - t339 * t409 - t340 * t404) * t364;
t286 = (-t335 * t414 - t336 * t409 - t337 * t404) * t364;
t1 = [t331 + t330 + t329, (t326 * t460 + t327 * t459 + t328 * t458) * t364, 0, 0, (t295 * t413 + t296 * t408 + t297 * t403 + (-t295 * t422 - t296 * t419 - t297 * t416) * t362) * t363, t289 * t380 + t290 * t377 + t291 * t374 + (-t326 * t386 - t327 * t384 - t328 * t382) * t424, -t289 * t457 - t290 * t455 - t291 * t453 + (t326 * t423 + t327 * t420 + t328 * t417) * t424, 1; t304, t288, 0, 0, (t298 * t413 + t299 * t408 + t300 * t403 + (-t298 * t422 - t299 * t419 - t300 * t416) * t362) * t363, t292 * t380 + t293 * t377 + t294 * t374 + t365, -t292 * t457 - t293 * t455 - t294 * t453 + t387, 0; 0, t286, 0, 0, (t308 * t413 + t309 * t408 + t310 * t403 + (-t308 * t422 - t309 * t419 - t310 * t416) * t362) * t363, t305 * t380 + t306 * t377 + t307 * t374 + (t335 * t385 + t336 * t383 + t337 * t381) * t424, -t305 * t457 - t306 * t455 - t307 * t453 + (-t332 * t422 - t333 * t419 - t334 * t416) * t424, 0; t304, t288, 0, 0, (t295 * t412 + t296 * t407 + t297 * t402 + (-t295 * t421 - t296 * t418 - t297 * t415) * t362) * t363, t289 * t379 + t290 * t376 + t291 * t373 + t365, -t289 * t456 - t290 * t454 - t291 * t452 + t387, 0; t328 + t327 + t326, (t329 * t460 + t330 * t459 + t331 * t458) * t364, 0, 0, (t298 * t412 + t299 * t407 + t300 * t402 + (-t298 * t421 - t299 * t418 - t300 * t415) * t362) * t363, t292 * t379 + t293 * t376 + t294 * t373 + (-t329 * t386 - t330 * t384 - t331 * t382) * t424, -t292 * t456 - t293 * t454 - t294 * t452 + (t329 * t423 + t330 * t420 + t331 * t417) * t424, 1; 0, t287, 0, 0, (t308 * t412 + t309 * t407 + t310 * t402 + (-t308 * t421 - t309 * t418 - t310 * t415) * t362) * t363, t305 * t379 + t306 * t376 + t307 * t373 + (t338 * t385 + t339 * t383 + t340 * t381) * t424, -t305 * t456 - t306 * t454 - t307 * t452 + (-t332 * t421 - t333 * t418 - t334 * t415) * t424, 0; 0, t286, 0, 0, (-t295 * t448 - t296 * t447 - t297 * t446 + (t295 * t451 + t296 * t450 + t297 * t449) * t362) * t363, -t289 * t399 - t290 * t397 - t291 * t395 + (t323 * t335 * t410 + t324 * t336 * t405 + t325 * t337 * t400) * t424, t332 * t289 + t333 * t290 + t334 * t291 + (-t323 * t413 - t324 * t408 - t325 * t403) * t424, 0; 0, t287, 0, 0, (-t298 * t448 - t299 * t447 - t300 * t446 + (t298 * t451 + t299 * t450 + t300 * t449) * t362) * t363, -t292 * t399 - t293 * t397 - t294 * t395 + (t323 * t378 + t324 * t375 + t325 * t372) * t424, t332 * t292 + t333 * t293 + t334 * t294 + (-t323 * t412 - t324 * t407 - t325 * t402) * t424, 0; 0, (t332 ^ 2 * t342 + t333 ^ 2 * t344 + t334 ^ 2 * t346) * t364, 0, 0, (-t308 * t448 - t309 * t447 - t310 * t446 + (t308 * t451 + t309 * t450 + t310 * t449) * t362) * t363, -t305 * t399 - t306 * t397 - t307 * t395 + (-t323 * t398 - t324 * t396 - t325 * t394) * t424, t332 * t305 + t333 * t306 + t334 * t307 + (t323 * t448 + t324 * t447 + t325 * t446) * t424, 1;];
tau_reg  = t1;
