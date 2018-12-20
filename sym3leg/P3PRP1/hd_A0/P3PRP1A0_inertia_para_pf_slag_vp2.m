% Calculate inertia matrix for parallel robot
% P3PRP1A0
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
% Datum: 2018-12-20 17:35
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MX = P3PRP1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:34:59
% EndTime: 2018-12-20 17:35:00
% DurationCPUTime: 1.00s
% Computational Cost: add. (5978->244), mult. (10997->365), div. (432->3), fcn. (5066->14), ass. (0->153)
t424 = 2 * qJ(3,1);
t423 = 2 * qJ(3,2);
t422 = 2 * qJ(3,3);
t383 = (pkin(2) ^ 2);
t421 = 1 + t383;
t368 = cos(qJ(2,3));
t420 = m(3) * t368;
t369 = cos(qJ(2,2));
t419 = m(3) * t369;
t370 = cos(qJ(2,1));
t418 = m(3) * t370;
t417 = mrSges(3,3) - mrSges(2,2);
t364 = legFrame(1,3);
t350 = sin(t364);
t416 = qJ(3,1) * t350;
t353 = cos(t364);
t415 = qJ(3,1) * t353;
t363 = legFrame(2,3);
t349 = sin(t363);
t414 = qJ(3,2) * t349;
t352 = cos(t363);
t413 = qJ(3,2) * t352;
t362 = legFrame(3,3);
t348 = sin(t362);
t412 = qJ(3,3) * t348;
t351 = cos(t362);
t411 = qJ(3,3) * t351;
t365 = sin(qJ(2,3));
t410 = t365 * t368;
t366 = sin(qJ(2,2));
t409 = t366 * t369;
t367 = sin(qJ(2,1));
t408 = t367 * t370;
t374 = qJ(3,3) ^ 2;
t407 = t374 + t383;
t375 = qJ(3,2) ^ 2;
t406 = t375 + t383;
t376 = qJ(3,1) ^ 2;
t405 = t376 + t383;
t404 = -0.2e1 * qJ(3,1) * t370;
t403 = -0.2e1 * qJ(3,2) * t369;
t402 = -0.2e1 * qJ(3,3) * t368;
t401 = pkin(2) * t416;
t400 = pkin(2) * t414;
t399 = pkin(2) * t412;
t398 = pkin(2) * t411;
t397 = pkin(2) * t413;
t396 = pkin(2) * t415;
t395 = 0.2e1 * pkin(2) * mrSges(3,1) + Ifges(3,2) + Ifges(2,3);
t354 = -m(3) * pkin(2) - mrSges(3,1);
t371 = xP(3);
t355 = sin(t371);
t356 = cos(t371);
t372 = mrSges(4,2);
t373 = mrSges(4,1);
t394 = -t355 * t372 + t356 * t373;
t393 = -t355 * t373 - t356 * t372;
t347 = -t376 + t421;
t392 = -t350 * t347 + 0.2e1 * t396;
t391 = -t376 * t350 - t396;
t390 = t353 * t376 - t401;
t346 = -t375 + t421;
t389 = -t349 * t346 + 0.2e1 * t397;
t388 = -t375 * t349 - t397;
t387 = t352 * t375 - t400;
t345 = -t374 + t421;
t386 = -t348 * t345 + 0.2e1 * t398;
t385 = -t374 * t348 - t398;
t384 = t351 * t374 - t399;
t382 = koppelP(1,1);
t381 = koppelP(2,1);
t380 = koppelP(3,1);
t379 = koppelP(1,2);
t378 = koppelP(2,2);
t377 = koppelP(3,2);
t361 = t370 ^ 2;
t360 = t369 ^ 2;
t359 = t368 ^ 2;
t357 = m(1) + m(2) + m(3);
t344 = mrSges(2,1) - t354;
t342 = -t355 * t379 + t356 * t382;
t341 = -t355 * t378 + t356 * t381;
t340 = -t355 * t377 + t356 * t380;
t339 = -t355 * t382 - t356 * t379;
t338 = -t355 * t381 - t356 * t378;
t337 = -t355 * t380 - t356 * t377;
t336 = t405 * m(3) + (mrSges(3,3) * t424) + t395;
t335 = t406 * m(3) + (mrSges(3,3) * t423) + t395;
t334 = t407 * m(3) + (mrSges(3,3) * t422) + t395;
t333 = t347 * t353 + 0.2e1 * t401;
t332 = t346 * t352 + 0.2e1 * t400;
t331 = t345 * t351 + 0.2e1 * t399;
t330 = (m(3) * qJ(3,1) + t417) * t367 + t344 * t370;
t329 = (m(3) * qJ(3,2) + t417) * t366 + t344 * t369;
t328 = (m(3) * qJ(3,3) + t417) * t365 + t344 * t368;
t327 = 0.1e1 / (pkin(2) * t408 * t424 + t347 * t361 - t405 - 0.1e1);
t326 = 0.1e1 / (pkin(2) * t409 * t423 + t346 * t360 - t406 - 0.1e1);
t325 = 0.1e1 / (pkin(2) * t410 * t422 + t345 * t359 - t407 - 0.1e1);
t324 = t353 * t404 + t367 * (pkin(2) * t353 + t416);
t323 = t350 * t404 + t367 * (pkin(2) * t350 - t415);
t322 = t352 * t403 + t366 * (pkin(2) * t352 + t414);
t321 = t349 * t403 + t366 * (pkin(2) * t349 - t413);
t320 = t351 * t402 + t365 * (pkin(2) * t351 + t412);
t319 = t348 * t402 + t365 * (pkin(2) * t348 - t411);
t318 = t390 * t370 - t367 * (t350 - t391);
t317 = t387 * t369 - t366 * (t349 - t388);
t316 = t384 * t368 - t365 * (t348 - t385);
t315 = t391 * t370 + t367 * (-t353 - t390);
t314 = t388 * t369 + t366 * (-t352 - t387);
t313 = t385 * t368 + t365 * (-t351 - t384);
t312 = -t333 * t408 + t383 * t350 + t392 * t361 + t350 - t396;
t311 = t333 * t361 - t353 * t383 + t392 * t408 - t353 - t401;
t310 = -t332 * t409 + t383 * t349 + t389 * t360 + t349 - t397;
t309 = t332 * t360 - t352 * t383 + t389 * t409 - t352 - t400;
t308 = -t331 * t410 + t383 * t348 + t386 * t359 + t348 - t398;
t307 = t331 * t359 - t351 * t383 + t386 * t410 - t351 - t399;
t306 = (t323 * t342 + t324 * t339) * t327;
t305 = (t321 * t341 + t322 * t338) * t326;
t304 = (t319 * t340 + t320 * t337) * t325;
t303 = (t315 * t339 + t318 * t342) * t327;
t302 = (t314 * t338 + t317 * t341) * t326;
t301 = (t313 * t337 + t316 * t340) * t325;
t300 = (t311 * t342 + t312 * t339) * t327;
t299 = (t309 * t341 + t310 * t338) * t326;
t298 = (t307 * t340 + t308 * t337) * t325;
t297 = (t324 * t354 + (-t312 * t370 + t315) * m(3)) * t327;
t296 = (t323 * t354 + (-t311 * t370 + t318) * m(3)) * t327;
t295 = (t322 * t354 + (-t310 * t369 + t314) * m(3)) * t326;
t294 = (t321 * t354 + (-t309 * t369 + t317) * m(3)) * t326;
t293 = (t320 * t354 + (-t308 * t368 + t313) * m(3)) * t325;
t292 = (t319 * t354 + (-t307 * t368 + t316) * m(3)) * t325;
t291 = (t312 * t357 - t315 * t418 + t324 * t330) * t327;
t290 = (t311 * t357 - t318 * t418 + t323 * t330) * t327;
t289 = (t310 * t357 - t314 * t419 + t322 * t329) * t326;
t288 = (t309 * t357 - t317 * t419 + t321 * t329) * t326;
t287 = (t308 * t357 - t313 * t420 + t320 * t328) * t325;
t286 = (t307 * t357 - t316 * t420 + t319 * t328) * t325;
t285 = (t312 * t330 + t315 * t354 + t324 * t336) * t327;
t284 = (t311 * t330 + t318 * t354 + t323 * t336) * t327;
t283 = (t310 * t329 + t314 * t354 + t322 * t335) * t326;
t282 = (t309 * t329 + t317 * t354 + t321 * t335) * t326;
t281 = (t308 * t328 + t313 * t354 + t320 * t334) * t325;
t280 = (t307 * t328 + t316 * t354 + t319 * t334) * t325;
t279 = t306 * t354 + (-t300 * t370 + t303) * m(3);
t278 = t305 * t354 + (-t299 * t369 + t302) * m(3);
t277 = t304 * t354 + (-t298 * t368 + t301) * m(3);
t276 = t300 * t357 - t303 * t418 + t306 * t330;
t275 = t299 * t357 - t302 * t419 + t305 * t329;
t274 = t298 * t357 - t301 * t420 + t304 * t328;
t273 = t300 * t330 + t303 * t354 + t306 * t336;
t272 = t299 * t329 + t302 * t354 + t305 * t335;
t271 = t298 * t328 + t301 * t354 + t304 * t334;
t1 = [m(4) + (t285 * t324 + t291 * t312 + t297 * t315) * t327 + (t283 * t322 + t289 * t310 + t295 * t314) * t326 + (t281 * t320 + t287 * t308 + t293 * t313) * t325 (t285 * t323 + t291 * t311 + t297 * t318) * t327 + (t283 * t321 + t289 * t309 + t295 * t317) * t326 + (t281 * t319 + t287 * t307 + t293 * t316) * t325, t281 * t304 + t283 * t305 + t285 * t306 + t287 * t298 + t289 * t299 + t291 * t300 + t293 * t301 + t295 * t302 + t297 * t303 + t393; (t284 * t324 + t290 * t312 + t296 * t315) * t327 + (t282 * t322 + t288 * t310 + t294 * t314) * t326 + (t280 * t320 + t286 * t308 + t292 * t313) * t325, m(4) + (t284 * t323 + t290 * t311 + t296 * t318) * t327 + (t282 * t321 + t288 * t309 + t294 * t317) * t326 + (t280 * t319 + t286 * t307 + t292 * t316) * t325, t280 * t304 + t282 * t305 + t284 * t306 + t286 * t298 + t288 * t299 + t290 * t300 + t292 * t301 + t294 * t302 + t296 * t303 + t394; (t273 * t324 + t276 * t312 + t279 * t315) * t327 + (t272 * t322 + t275 * t310 + t278 * t314) * t326 + (t271 * t320 + t274 * t308 + t277 * t313) * t325 + t393 (t273 * t323 + t276 * t311 + t279 * t318) * t327 + (t272 * t321 + t275 * t309 + t278 * t317) * t326 + (t271 * t319 + t274 * t307 + t277 * t316) * t325 + t394, t271 * t304 + t272 * t305 + t273 * t306 + t274 * t298 + t275 * t299 + t276 * t300 + t277 * t301 + t278 * t302 + t279 * t303 + Ifges(4,3);];
MX  = t1;
