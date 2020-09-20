% Calculate minimal parameter regressor of vector of centrifugal and coriolis load for parallel robot
% P3PRRR2G3A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
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
% tau_reg [3x8]
%   minimal parameter regressor of vector of coriolis and centrifugal joint torques
%   in task space

% Quelle: HybrDyn-Toolbox
% Datum: 2020-03-09 21:20
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = P3PRRR2G3A0_coriolisvec_para_pf_regmin(xP, xDP, qJ, legFrame, ...
  koppelP, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_regmin: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_regmin: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_regmin: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_regmin: pkin has to be [2x1] (double)');
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_regmin: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRR2G3A0_coriolisvec_para_pf_regmin: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From invdyn_para_plfcoord_tauCreg_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:20:07
% EndTime: 2020-03-09 21:20:07
% DurationCPUTime: 0.74s
% Computational Cost: add. (4396->121), mult. (3558->238), div. (1356->9), fcn. (2406->18), ass. (0->132)
t372 = sin(qJ(3,3));
t361 = 0.1e1 / t372 ^ 2;
t375 = cos(qJ(3,3));
t445 = t361 * t375;
t373 = sin(qJ(3,2));
t364 = 0.1e1 / t373 ^ 2;
t376 = cos(qJ(3,2));
t444 = t364 * t376;
t374 = sin(qJ(3,1));
t367 = 0.1e1 / t374 ^ 2;
t377 = cos(qJ(3,1));
t443 = t367 * t377;
t442 = 2 * pkin(2);
t360 = 0.1e1 / t372;
t363 = 0.1e1 / t373;
t366 = 0.1e1 / t374;
t383 = 0.1e1 / pkin(1) ^ 2;
t416 = t361 * t383;
t381 = 1 / pkin(2);
t382 = 1 / pkin(1);
t411 = t381 * t382;
t357 = -legFrame(3,2) + qJ(2,3);
t345 = sin(t357);
t348 = cos(t357);
t378 = xDP(2);
t379 = xDP(1);
t351 = qJ(3,3) + t357;
t339 = sin(t351);
t342 = cos(t351);
t392 = t339 * t379 - t342 * t378;
t321 = t392 * pkin(2) + (t345 * t379 - t348 * t378) * pkin(1);
t432 = t321 * t360;
t318 = t411 * t432;
t426 = t392 * t360;
t324 = t382 * t426;
t313 = -t324 + t318;
t435 = t313 * t321;
t306 = t416 * t435;
t417 = t360 * t375;
t309 = -pkin(2) * t313 + t392 * t417;
t397 = (t375 * pkin(1) + pkin(2)) * t381 * t435;
t380 = pkin(2) ^ 2;
t438 = (t313 * t380 + ((-t324 + t318 / 0.2e1) * t375 * t442 - t426) * pkin(1)) * t381;
t297 = t306 + (-t397 - (-t309 - t438) * t392) * t416;
t441 = t297 * t360;
t414 = t364 * t383;
t358 = -legFrame(2,2) + qJ(2,2);
t346 = sin(t358);
t349 = cos(t358);
t352 = qJ(3,2) + t358;
t340 = sin(t352);
t343 = cos(t352);
t391 = t340 * t379 - t343 * t378;
t322 = t391 * pkin(2) + (t346 * t379 - t349 * t378) * pkin(1);
t431 = t322 * t363;
t319 = t411 * t431;
t425 = t391 * t363;
t325 = t382 * t425;
t315 = -t325 + t319;
t434 = t315 * t322;
t307 = t414 * t434;
t415 = t363 * t376;
t310 = -pkin(2) * t315 + t391 * t415;
t395 = (t376 * pkin(1) + pkin(2)) * t381 * t434;
t437 = (t315 * t380 + ((-t325 + t319 / 0.2e1) * t376 * t442 - t425) * pkin(1)) * t381;
t298 = t307 + (-t395 - (-t310 - t437) * t391) * t414;
t440 = t298 * t363;
t412 = t367 * t383;
t359 = -legFrame(1,2) + qJ(2,1);
t347 = sin(t359);
t350 = cos(t359);
t353 = qJ(3,1) + t359;
t341 = sin(t353);
t344 = cos(t353);
t390 = t341 * t379 - t344 * t378;
t323 = t390 * pkin(2) + (t347 * t379 - t350 * t378) * pkin(1);
t430 = t323 * t366;
t320 = t411 * t430;
t424 = t390 * t366;
t326 = t382 * t424;
t317 = -t326 + t320;
t433 = t317 * t323;
t308 = t412 * t433;
t413 = t366 * t377;
t311 = -t317 * pkin(2) + t390 * t413;
t393 = (t377 * pkin(1) + pkin(2)) * t381 * t433;
t436 = (t317 * t380 + ((-t326 + t320 / 0.2e1) * t377 * t442 - t424) * pkin(1)) * t381;
t299 = t308 + (-t393 - (-t311 - t436) * t390) * t412;
t439 = t299 * t366;
t327 = t392 ^ 2;
t429 = t327 * t361;
t328 = t391 ^ 2;
t428 = t328 * t364;
t329 = t390 ^ 2;
t427 = t329 * t367;
t423 = t339 * t360;
t422 = t340 * t363;
t421 = t341 * t366;
t420 = t342 * t360;
t419 = t343 * t363;
t418 = t344 * t366;
t300 = t309 * t392 * t416 + t306;
t410 = t300 * t417;
t301 = t310 * t391 * t414 + t307;
t409 = t301 * t415;
t302 = t311 * t390 * t412 + t308;
t408 = t302 * t413;
t312 = -0.2e1 * t324 + t318;
t407 = t312 * t432;
t314 = -0.2e1 * t325 + t319;
t406 = t314 * t431;
t316 = -0.2e1 * t326 + t320;
t405 = t316 * t430;
t404 = t327 * t360 * t445;
t403 = t328 * t363 * t444;
t402 = t329 * t366 * t443;
t294 = 0.2e1 * t306 + (-t397 - (-0.2e1 * t309 - t438) * t392) * t416;
t401 = t294 * t417;
t295 = 0.2e1 * t307 + (-t395 - (-0.2e1 * t310 - t437) * t391) * t414;
t400 = t295 * t415;
t296 = 0.2e1 * t308 + (-t393 - (-0.2e1 * t311 - t436) * t390) * t412;
t399 = t296 * t413;
t398 = t312 * t321 * t445;
t396 = t314 * t322 * t444;
t394 = t316 * t323 * t443;
t338 = -pkin(1) * t350 - pkin(2) * t344;
t337 = -pkin(1) * t349 - pkin(2) * t343;
t336 = -pkin(1) * t348 - pkin(2) * t342;
t335 = pkin(1) * t347 + pkin(2) * t341;
t334 = pkin(1) * t346 + pkin(2) * t340;
t333 = pkin(1) * t345 + pkin(2) * t339;
t1 = [0, (-t300 * t423 - t301 * t422 - t302 * t421) * t382, 0, 0, (-t297 * t423 - t298 * t422 - t299 * t421 + (t333 * t441 + t334 * t440 + t335 * t439) * t381) * t382, -t339 * t401 - t340 * t400 - t341 * t399 + (t333 * t410 + t334 * t409 + t335 * t408 + (t333 * t429 + t334 * t428 + t335 * t427) * t383 + (t339 * t407 + t340 * t406 + t341 * t405) * t382) * t381, t339 * t294 + t340 * t295 + t341 * t296 + (-t333 * t300 - t334 * t301 - t335 * t302 + (t333 * t404 + t334 * t403 + t335 * t402) * t383 + (t339 * t398 + t340 * t396 + t341 * t394) * t382) * t381, 0; 0, (t300 * t420 + t301 * t419 + t302 * t418) * t382, 0, 0, (t297 * t420 + t298 * t419 + t299 * t418 + (t336 * t441 + t337 * t440 + t338 * t439) * t381) * t382, t342 * t401 + t343 * t400 + t344 * t399 + (t336 * t410 + t337 * t409 + t338 * t408 + (t336 * t429 + t337 * t428 + t338 * t427) * t383 + (-t342 * t407 - t343 * t406 - t344 * t405) * t382) * t381, -t342 * t294 - t343 * t295 - t344 * t296 + (-t336 * t300 - t337 * t301 - t338 * t302 + (t336 * t404 + t337 * t403 + t338 * t402) * t383 + (-t342 * t398 - t343 * t396 - t344 * t394) * t382) * t381, 0; 0, 0, 0, 0, 0, 0, 0, 0;];
tau_reg  = t1;
