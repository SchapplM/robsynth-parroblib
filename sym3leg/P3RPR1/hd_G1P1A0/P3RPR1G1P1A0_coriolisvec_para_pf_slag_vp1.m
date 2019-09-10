% Calculate vector of centrifugal and coriolis load on the joints for
% P3RPR1G1P1A0
% Use Code from Maple symbolic Code Generation
%
% Input:
% xP [3x1]
%   Generalized platform coordinates
% xDP [3x1]
%   Generalized platform velocities
% qJ [2x3]
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
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d3]';
% m [3x1]
%   mass of all robot links (leg links until cut joint, platform)
% rSges [3x3]
%   center of mass of all robot links (in body frames)
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: x-, y-, z-coordinates
% Icges [3x6]
%   inertia of all robot links about their respective center of mass, in body frames
%   rows: links of the robot (leg links until cut joint, platform)
%   columns: xx, yy, zz, xy, xz, yz (see inertiavector2matrix.m)
%
% Output:
% taucX [3x1]
%   forces required to compensate Coriolis and centrifugal joint torques
%   in platform coordinates

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-03 14:58
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taucX = P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1(xP, xDP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(2,3),zeros(3,3),zeros(3,3),zeros(4,1),zeros(2+1,1),zeros(2+1,3),zeros(2+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(xDP) && all(size(xDP) == [3 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: xDP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [2 3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: qJ has to be [2x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: m has to be [3x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [3,3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: rSges has to be [3x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [3 6]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: Icges has to be [3x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPR1G1P1A0_coriolisvec_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From coriolisvec_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:58:04
% EndTime: 2019-05-03 14:58:05
% DurationCPUTime: 1.04s
% Computational Cost: add. (3385->170), mult. (5014->314), div. (666->3), fcn. (3892->14), ass. (0->158)
t432 = 2 * rSges(2,3);
t367 = xDP(3);
t352 = t367 ^ 2;
t371 = xP(3);
t349 = sin(t371);
t350 = cos(t371);
t380 = koppelP(3,2);
t383 = koppelP(3,1);
t333 = -t349 * t380 + t350 * t383;
t368 = xDP(2);
t318 = t333 * t367 + t368;
t330 = t349 * t383 + t350 * t380;
t369 = xDP(1);
t321 = -t330 * t367 + t369;
t360 = sin(qJ(1,3));
t363 = cos(qJ(1,3));
t370 = pkin(1) + pkin(2);
t336 = -t363 * qJ(2,3) + t360 * t370;
t339 = t360 * qJ(2,3) + t370 * t363;
t357 = legFrame(3,3);
t342 = sin(t357);
t345 = cos(t357);
t312 = -t342 * t336 + t339 * t345;
t375 = 0.1e1 / qJ(2,3);
t410 = t312 * t375;
t309 = t336 * t345 + t342 * t339;
t413 = t309 * t375;
t276 = t318 * t413 + t321 * t410;
t431 = 0.2e1 * t276;
t381 = koppelP(2,2);
t384 = koppelP(2,1);
t334 = -t349 * t381 + t350 * t384;
t319 = t334 * t367 + t368;
t331 = t349 * t384 + t350 * t381;
t322 = -t331 * t367 + t369;
t361 = sin(qJ(1,2));
t364 = cos(qJ(1,2));
t337 = -t364 * qJ(2,2) + t361 * t370;
t340 = t361 * qJ(2,2) + t370 * t364;
t358 = legFrame(2,3);
t343 = sin(t358);
t346 = cos(t358);
t313 = -t343 * t337 + t340 * t346;
t377 = 0.1e1 / qJ(2,2);
t409 = t313 * t377;
t310 = t337 * t346 + t343 * t340;
t412 = t310 * t377;
t277 = t319 * t412 + t322 * t409;
t430 = 0.2e1 * t277;
t382 = koppelP(1,2);
t385 = koppelP(1,1);
t335 = -t349 * t382 + t350 * t385;
t320 = t335 * t367 + t368;
t332 = t349 * t385 + t350 * t382;
t323 = -t332 * t367 + t369;
t362 = sin(qJ(1,1));
t365 = cos(qJ(1,1));
t338 = -t365 * qJ(2,1) + t362 * t370;
t341 = t362 * qJ(2,1) + t370 * t365;
t359 = legFrame(1,3);
t344 = sin(t359);
t347 = cos(t359);
t314 = -t344 * t338 + t341 * t347;
t379 = 0.1e1 / qJ(2,1);
t408 = t314 * t379;
t311 = t338 * t347 + t344 * t341;
t411 = t311 * t379;
t278 = t320 * t411 + t323 * t408;
t429 = 0.2e1 * t278;
t428 = 0.2e1 * t370;
t427 = m(2) * (rSges(2,3) + qJ(2,3));
t426 = m(2) * (rSges(2,3) + qJ(2,2));
t425 = m(2) * (rSges(2,3) + qJ(2,1));
t366 = pkin(1) + rSges(2,1);
t424 = m(2) * t366;
t423 = m(2) * t375;
t422 = m(2) * t377;
t421 = m(2) * t379;
t420 = m(3) * t352;
t324 = t342 * t363 + t345 * t360;
t325 = -t342 * t360 + t345 * t363;
t291 = (t318 * t324 + t321 * t325) * t375;
t374 = qJ(2,3) ^ 2;
t389 = (pkin(1) ^ 2);
t401 = -t389 + (-2 * pkin(1) - pkin(2)) * pkin(2);
t416 = t291 * t375;
t264 = (t276 * t428 + (-t374 + t401) * t291) * t416;
t267 = (-t370 * t291 + t431) * t416;
t393 = (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + Icges(2,2) + Icges(1,3);
t397 = rSges(2,3) ^ 2 + t389 + (2 * pkin(1) + rSges(2,1)) * rSges(2,1);
t315 = (qJ(2,3) * t432 + t374 + t397) * m(2) + t393;
t258 = t264 * t424 - t315 * t267;
t419 = t258 * t375;
t326 = t343 * t364 + t346 * t361;
t327 = -t343 * t361 + t346 * t364;
t292 = (t319 * t326 + t322 * t327) * t377;
t376 = qJ(2,2) ^ 2;
t415 = t292 * t377;
t265 = (t277 * t428 + (-t376 + t401) * t292) * t415;
t268 = (-t370 * t292 + t430) * t415;
t316 = (qJ(2,2) * t432 + t376 + t397) * m(2) + t393;
t259 = t265 * t424 - t316 * t268;
t418 = t259 * t377;
t328 = t344 * t365 + t347 * t362;
t329 = -t344 * t362 + t347 * t365;
t293 = (t320 * t328 + t323 * t329) * t379;
t378 = qJ(2,1) ^ 2;
t414 = t293 * t379;
t266 = (t278 * t428 + (-t378 + t401) * t293) * t414;
t269 = (-t370 * t293 + t429) * t414;
t317 = (qJ(2,1) * t432 + t378 + t397) * m(2) + t393;
t260 = t266 * t424 - t317 * t269;
t417 = t260 * t379;
t407 = t375 * t352;
t406 = t377 * t352;
t405 = t379 * t352;
t404 = t291 ^ 2 * t427;
t403 = t292 ^ 2 * t426;
t402 = t293 ^ 2 * t425;
t400 = t375 * t404;
t399 = t377 * t403;
t398 = t379 * t402;
t396 = t291 * t427 * t431;
t395 = t292 * t426 * t430;
t394 = t293 * t425 * t429;
t392 = t375 * t396;
t391 = t377 * t395;
t390 = t379 * t394;
t373 = rSges(3,1);
t372 = rSges(3,2);
t302 = (-t329 * t366 + t314) * t421;
t301 = (-t328 * t366 + t311) * t421;
t300 = (-t327 * t366 + t313) * t422;
t299 = (-t326 * t366 + t310) * t422;
t298 = (-t325 * t366 + t312) * t423;
t297 = (-t324 * t366 + t309) * t423;
t296 = (t328 * t335 - t329 * t332) * t379;
t295 = (t326 * t334 - t327 * t331) * t377;
t294 = (t324 * t333 - t325 * t330) * t375;
t287 = (-t314 * t424 + t317 * t329) * t379;
t286 = (-t311 * t424 + t317 * t328) * t379;
t285 = (-t313 * t424 + t316 * t327) * t377;
t284 = (-t310 * t424 + t316 * t326) * t377;
t283 = (-t312 * t424 + t315 * t325) * t375;
t282 = (-t309 * t424 + t315 * t324) * t375;
t281 = (t311 * t335 - t314 * t332) * t379;
t280 = (t310 * t334 - t313 * t331) * t377;
t279 = (t309 * t333 - t312 * t330) * t375;
t275 = (-t296 * t366 + t281) * m(2);
t274 = (-t295 * t366 + t280) * m(2);
t273 = (-t294 * t366 + t279) * m(2);
t272 = -t281 * t424 + t296 * t317;
t271 = -t280 * t424 + t295 * t316;
t270 = -t279 * t424 + t294 * t315;
t263 = (t269 * t366 - t266) * m(2);
t262 = (t268 * t366 - t265) * m(2);
t261 = (t267 * t366 - t264) * m(2);
t1 = [(-(t287 * t329 + t302 * t314) * t335 - (t287 * t328 + t302 * t311) * t332) * t405 + t329 * t417 + t263 * t408 + t329 * t390 - t314 * t398 + (-(t285 * t327 + t300 * t313) * t334 - (t285 * t326 + t300 * t310) * t331) * t406 + t327 * t418 + t262 * t409 + t327 * t391 - t313 * t399 + (-(t283 * t325 + t298 * t312) * t333 - (t283 * t324 + t298 * t309) * t330) * t407 + t325 * t419 + t261 * t410 + t325 * t392 - t312 * t400 - (-t349 * t372 + t350 * t373) * t420; (-(t286 * t329 + t301 * t314) * t335 - (t286 * t328 + t301 * t311) * t332) * t405 + t328 * t417 + t263 * t411 + t328 * t390 - t311 * t398 + (-(t284 * t327 + t299 * t313) * t334 - (t284 * t326 + t299 * t310) * t331) * t406 + t326 * t418 + t262 * t412 + t326 * t391 - t310 * t399 + (-(t282 * t325 + t297 * t312) * t333 - (t282 * t324 + t297 * t309) * t330) * t407 + t324 * t419 + t261 * t413 + t324 * t392 - t309 * t400 - (t349 * t373 + t350 * t372) * t420; (-(t272 * t329 + t275 * t314) * t335 - (t272 * t328 + t275 * t311) * t332) * t405 + (-(t271 * t327 + t274 * t313) * t334 - (t271 * t326 + t274 * t310) * t331) * t406 + (-(t270 * t325 + t273 * t312) * t333 - (t270 * t324 + t273 * t309) * t330) * t407 + (t260 + t394) * t296 + (t259 + t395) * t295 + (t258 + t396) * t294 + (t263 - t402) * t281 + (t262 - t403) * t280 + (t261 - t404) * t279;];
taucX  = t1;
