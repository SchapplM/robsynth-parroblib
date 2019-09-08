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
% Datum: 2019-05-03 14:53
% Revision: abbb0d669c4fc7889a31e0cf750ab51a4f2eb1ce (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPP1G1P1A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPP1G1P1A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-03 14:52:11
% EndTime: 2019-05-03 14:52:11
% DurationCPUTime: 0.67s
% Computational Cost: add. (5213->215), mult. (6884->331), div. (360->3), fcn. (2528->14), ass. (0->157)
t434 = 2 * pkin(1);
t433 = 2 * qJ(3,1);
t432 = 2 * qJ(3,2);
t431 = 2 * qJ(3,3);
t430 = 2 * mrSges(3,2) + 2 * mrSges(2,3);
t414 = pkin(1) ^ 2;
t429 = 1 + t414;
t428 = (-mrSges(2,2) + mrSges(3,3));
t388 = -qJ(3,3) - pkin(1);
t392 = sin(qJ(1,3));
t427 = t392 * t388;
t389 = -qJ(3,2) - pkin(1);
t393 = sin(qJ(1,2));
t426 = t393 * t389;
t390 = -qJ(3,1) - pkin(1);
t394 = sin(qJ(1,1));
t425 = t394 * t390;
t395 = cos(qJ(1,3));
t424 = t395 * qJ(2,3);
t396 = cos(qJ(1,2));
t423 = t396 * qJ(2,2);
t397 = cos(qJ(1,1));
t422 = t397 * qJ(2,1);
t421 = t390 * t422;
t420 = t389 * t423;
t419 = t388 * t424;
t398 = (m(2) + m(3));
t418 = t398 * t414 + Ifges(2,1) + Ifges(3,1) + Ifges(1,3);
t406 = qJ(3,1) ^ 2;
t362 = pkin(1) * t433 + t406 + t429;
t402 = qJ(3,3) ^ 2;
t360 = pkin(1) * t431 + t402 + t429;
t404 = qJ(3,2) ^ 2;
t361 = pkin(1) * t432 + t404 + t429;
t399 = xP(3);
t377 = sin(t399);
t378 = cos(t399);
t400 = mrSges(4,2);
t401 = mrSges(4,1);
t417 = -t377 * t400 + t378 * t401;
t416 = -m(2) * pkin(1) - t428;
t415 = -t377 * t401 - t378 * t400;
t413 = koppelP(1,1);
t412 = koppelP(2,1);
t411 = koppelP(3,1);
t410 = koppelP(1,2);
t409 = koppelP(2,2);
t408 = koppelP(3,2);
t407 = qJ(2,1) ^ 2;
t405 = qJ(2,2) ^ 2;
t403 = qJ(2,3) ^ 2;
t387 = legFrame(1,3);
t386 = legFrame(2,3);
t385 = legFrame(3,3);
t384 = 0.1e1 + t407;
t383 = 0.1e1 + t405;
t382 = 0.1e1 + t403;
t375 = cos(t387);
t374 = cos(t386);
t373 = cos(t385);
t372 = sin(t387);
t371 = sin(t386);
t370 = sin(t385);
t369 = m(3) * qJ(2,1) + mrSges(3,2);
t368 = m(3) * qJ(2,2) + mrSges(3,2);
t367 = m(3) * qJ(2,3) + mrSges(3,2);
t365 = qJ(2,1) * t425;
t364 = qJ(2,2) * t426;
t363 = qJ(2,3) * t427;
t359 = t390 * m(3) + t416;
t358 = t389 * m(3) + t416;
t357 = t388 * m(3) + t416;
t356 = 0.1e1 / (t407 + t362);
t355 = 0.1e1 / (t405 + t361);
t354 = 0.1e1 / (t403 + t360);
t353 = t394 * qJ(2,1) - t390 * t397;
t352 = t393 * qJ(2,2) - t389 * t396;
t351 = t392 * qJ(2,3) - t388 * t395;
t350 = t422 + t425;
t349 = t423 + t426;
t348 = t424 + t427;
t347 = -t377 * t410 + t378 * t413;
t346 = -t377 * t409 + t378 * t412;
t345 = -t377 * t408 + t378 * t411;
t344 = -t377 * t413 - t378 * t410;
t343 = -t377 * t412 - t378 * t409;
t342 = -t377 * t411 - t378 * t408;
t341 = -t384 * t397 - t365;
t340 = -t383 * t396 - t364;
t339 = -t382 * t395 - t363;
t338 = t394 * t384 - t421;
t337 = t393 * t383 - t420;
t336 = t392 * t382 - t419;
t335 = t362 * t397 - t365;
t334 = t361 * t396 - t364;
t333 = t360 * t395 - t363;
t332 = t394 * t362 + t421;
t331 = t393 * t361 + t420;
t330 = t392 * t360 + t419;
t329 = t372 * t350 + t353 * t375;
t328 = t371 * t349 + t352 * t374;
t327 = t370 * t348 + t351 * t373;
t326 = t350 * t375 - t372 * t353;
t325 = t349 * t374 - t371 * t352;
t324 = t348 * t373 - t370 * t351;
t323 = m(3) * t406 + (mrSges(3,3) * t433) + t398 * t407 + (m(3) * qJ(3,1) + t428) * t434 + qJ(2,1) * t430 + t418;
t322 = m(3) * t404 + (mrSges(3,3) * t432) + t398 * t405 + (m(3) * qJ(3,2) + t428) * t434 + qJ(2,2) * t430 + t418;
t321 = m(3) * t402 + (mrSges(3,3) * t431) + t398 * t403 + (m(3) * qJ(3,3) + t428) * t434 + qJ(2,3) * t430 + t418;
t320 = t372 * t338 + t341 * t375;
t319 = t371 * t337 + t340 * t374;
t318 = t370 * t336 + t339 * t373;
t317 = t338 * t375 - t372 * t341;
t316 = t337 * t374 - t371 * t340;
t315 = t336 * t373 - t370 * t339;
t314 = -t372 * t332 + t335 * t375;
t313 = -t371 * t331 + t334 * t374;
t312 = -t370 * t330 + t333 * t373;
t311 = t332 * t375 + t372 * t335;
t310 = t331 * t374 + t371 * t334;
t309 = t330 * t373 + t370 * t333;
t308 = (t320 * t398 + t329 * t359) * t356;
t307 = (t319 * t398 + t328 * t358) * t355;
t306 = (t318 * t398 + t327 * t357) * t354;
t305 = (t317 * t398 + t326 * t359) * t356;
t304 = (t316 * t398 + t325 * t358) * t355;
t303 = (t315 * t398 + t324 * t357) * t354;
t302 = (m(3) * t311 + t329 * t369) * t356;
t301 = (m(3) * t310 + t328 * t368) * t355;
t300 = (m(3) * t309 + t327 * t367) * t354;
t299 = (m(3) * t314 + t326 * t369) * t356;
t298 = (m(3) * t313 + t325 * t368) * t355;
t297 = (m(3) * t312 + t324 * t367) * t354;
t296 = (t326 * t344 + t329 * t347) * t356;
t295 = (t325 * t343 + t328 * t346) * t355;
t294 = (t324 * t342 + t327 * t345) * t354;
t293 = (t317 * t344 + t320 * t347) * t356;
t292 = (t316 * t343 + t319 * t346) * t355;
t291 = (t315 * t342 + t318 * t345) * t354;
t290 = (t311 * t347 + t314 * t344) * t356;
t289 = (t310 * t346 + t313 * t343) * t355;
t288 = (t309 * t345 + t312 * t342) * t354;
t287 = (t311 * t369 + t320 * t359 + t323 * t329) * t356;
t286 = (t310 * t368 + t319 * t358 + t322 * t328) * t355;
t285 = (t309 * t367 + t318 * t357 + t321 * t327) * t354;
t284 = (t314 * t369 + t317 * t359 + t323 * t326) * t356;
t283 = (t313 * t368 + t316 * t358 + t322 * t325) * t355;
t282 = (t312 * t367 + t315 * t357 + t321 * t324) * t354;
t281 = t293 * t398 + t296 * t359;
t280 = t292 * t398 + t295 * t358;
t279 = t291 * t398 + t294 * t357;
t278 = t290 * m(3) + t296 * t369;
t277 = t289 * m(3) + t295 * t368;
t276 = t288 * m(3) + t294 * t367;
t275 = t290 * t369 + t293 * t359 + t296 * t323;
t274 = t289 * t368 + t292 * t358 + t295 * t322;
t273 = t288 * t367 + t291 * t357 + t294 * t321;
t1 = [m(4) + (t284 * t326 + t299 * t314 + t305 * t317) * t356 + (t283 * t325 + t298 * t313 + t304 * t316) * t355 + (t282 * t324 + t297 * t312 + t303 * t315) * t354, (t284 * t329 + t299 * t311 + t305 * t320) * t356 + (t283 * t328 + t298 * t310 + t304 * t319) * t355 + (t282 * t327 + t297 * t309 + t303 * t318) * t354, t282 * t294 + t283 * t295 + t284 * t296 + t297 * t288 + t298 * t289 + t299 * t290 + t303 * t291 + t304 * t292 + t305 * t293 + t415; (t287 * t326 + t302 * t314 + t308 * t317) * t356 + (t286 * t325 + t301 * t313 + t307 * t316) * t355 + (t285 * t324 + t300 * t312 + t306 * t315) * t354, m(4) + (t287 * t329 + t302 * t311 + t308 * t320) * t356 + (t286 * t328 + t301 * t310 + t307 * t319) * t355 + (t285 * t327 + t300 * t309 + t306 * t318) * t354, t285 * t294 + t286 * t295 + t287 * t296 + t300 * t288 + t301 * t289 + t302 * t290 + t306 * t291 + t307 * t292 + t308 * t293 + t417; (t275 * t326 + t278 * t314 + t281 * t317) * t356 + (t274 * t325 + t277 * t313 + t280 * t316) * t355 + (t273 * t324 + t276 * t312 + t279 * t315) * t354 + t415, (t275 * t329 + t278 * t311 + t281 * t320) * t356 + (t274 * t328 + t277 * t310 + t280 * t319) * t355 + (t273 * t327 + t276 * t309 + t279 * t318) * t354 + t417, t273 * t294 + t274 * t295 + t275 * t296 + t276 * t288 + t277 * t289 + t278 * t290 + t279 * t291 + t280 * t292 + t281 * t293 + Ifges(4,3);];
MX  = t1;
