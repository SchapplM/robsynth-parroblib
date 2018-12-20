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
% Datum: 2018-12-20 17:39
% Revision: f9720dcdc4676342702b46a014e894344751412a
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function MX = P3PRP2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRP2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-12-20 17:39:00
% EndTime: 2018-12-20 17:39:01
% DurationCPUTime: 0.91s
% Computational Cost: add. (5978->244), mult. (10997->369), div. (432->3), fcn. (5066->14), ass. (0->154)
t429 = 2 * mrSges(3,3);
t393 = (pkin(2) ^ 2);
t428 = -t393 - 1;
t378 = cos(qJ(2,3));
t427 = m(3) * t378;
t379 = cos(qJ(2,2));
t426 = m(3) * t379;
t380 = cos(qJ(2,1));
t425 = m(3) * t380;
t424 = mrSges(3,3) - mrSges(2,2);
t374 = legFrame(1,3);
t360 = sin(t374);
t423 = qJ(3,1) * t360;
t363 = cos(t374);
t422 = qJ(3,1) * t363;
t373 = legFrame(2,3);
t359 = sin(t373);
t421 = qJ(3,2) * t359;
t362 = cos(t373);
t420 = qJ(3,2) * t362;
t372 = legFrame(3,3);
t358 = sin(t372);
t419 = qJ(3,3) * t358;
t361 = cos(t372);
t418 = qJ(3,3) * t361;
t375 = sin(qJ(2,3));
t417 = t375 * t378;
t376 = sin(qJ(2,2));
t416 = t376 * t379;
t377 = sin(qJ(2,1));
t415 = t377 * t380;
t414 = t378 * qJ(3,3);
t413 = t379 * qJ(3,2);
t412 = t380 * qJ(3,1);
t347 = pkin(2) * t419;
t384 = qJ(3,3) ^ 2;
t411 = t384 * t361 + t347;
t348 = pkin(2) * t421;
t385 = qJ(3,2) ^ 2;
t410 = t385 * t362 + t348;
t349 = pkin(2) * t423;
t386 = qJ(3,1) ^ 2;
t409 = t386 * t363 + t349;
t408 = t384 + t393;
t407 = t385 + t393;
t406 = t386 + t393;
t405 = 0.2e1 * t414;
t404 = 0.2e1 * t413;
t403 = 0.2e1 * t412;
t402 = pkin(2) * t422;
t401 = pkin(2) * t420;
t400 = pkin(2) * t418;
t399 = 0.2e1 * pkin(2) * mrSges(3,1) + Ifges(3,2) + Ifges(2,3);
t364 = -m(3) * pkin(2) - mrSges(3,1);
t381 = xP(3);
t365 = sin(t381);
t366 = cos(t381);
t382 = mrSges(4,2);
t383 = mrSges(4,1);
t398 = -t365 * t382 + t366 * t383;
t397 = -t365 * t383 - t366 * t382;
t396 = t360 * t386 - t402;
t395 = t359 * t385 - t401;
t394 = t358 * t384 - t400;
t392 = koppelP(1,1);
t391 = koppelP(2,1);
t390 = koppelP(3,1);
t389 = koppelP(1,2);
t388 = koppelP(2,2);
t387 = koppelP(3,2);
t371 = t380 ^ 2;
t370 = t379 ^ 2;
t369 = t378 ^ 2;
t367 = m(1) + m(2) + m(3);
t357 = -t386 - t428;
t356 = -t385 - t428;
t355 = -t384 - t428;
t354 = mrSges(2,1) - t364;
t346 = -t365 * t389 + t366 * t392;
t345 = -t365 * t388 + t366 * t391;
t344 = -t365 * t387 + t366 * t390;
t343 = -t365 * t392 - t366 * t389;
t342 = -t365 * t391 - t366 * t388;
t341 = -t365 * t390 - t366 * t387;
t340 = t406 * m(3) + qJ(3,1) * t429 + t399;
t339 = t407 * m(3) + qJ(3,2) * t429 + t399;
t338 = t408 * m(3) + qJ(3,3) * t429 + t399;
t337 = t357 * t360 + 0.2e1 * t402;
t336 = t356 * t359 + 0.2e1 * t401;
t335 = t355 * t358 + 0.2e1 * t400;
t334 = t363 * t357 - 0.2e1 * t349;
t333 = t362 * t356 - 0.2e1 * t348;
t332 = t361 * t355 - 0.2e1 * t347;
t331 = (m(3) * qJ(3,1) + t424) * t377 + t354 * t380;
t330 = (m(3) * qJ(3,2) + t424) * t376 + t354 * t379;
t329 = (m(3) * qJ(3,3) + t424) * t375 + t354 * t378;
t328 = 0.1e1 / (t377 * pkin(2) * t403 + t357 * t371 - t406 - 0.1e1);
t327 = 0.1e1 / (t376 * pkin(2) * t404 + t356 * t370 - t407 - 0.1e1);
t326 = 0.1e1 / (t375 * pkin(2) * t405 + t355 * t369 - t408 - 0.1e1);
t325 = -0.2e1 * t363 * t412 + t377 * (pkin(2) * t363 - t423);
t324 = t360 * t403 - t377 * (pkin(2) * t360 + t422);
t323 = -0.2e1 * t362 * t413 + t376 * (pkin(2) * t362 - t421);
t322 = t359 * t404 - t376 * (pkin(2) * t359 + t420);
t321 = -0.2e1 * t361 * t414 + t375 * (pkin(2) * t361 - t419);
t320 = t358 * t405 - t375 * (pkin(2) * t358 + t418);
t319 = t396 * t380 - t377 * (t363 + t409);
t318 = t395 * t379 - t376 * (t362 + t410);
t317 = t394 * t378 - t375 * (t361 + t411);
t316 = t409 * t380 - t377 * (-t360 - t396);
t315 = t410 * t379 - t376 * (-t359 - t395);
t314 = t411 * t378 - t375 * (-t358 - t394);
t313 = -t334 * t415 + t337 * t371 - t360 * t393 - t360 - t402;
t312 = -t333 * t416 + t336 * t370 - t359 * t393 - t359 - t401;
t311 = -t332 * t417 + t335 * t369 - t358 * t393 - t358 - t400;
t310 = t334 * t371 + t337 * t415 - t393 * t363 + t349 - t363;
t309 = t333 * t370 + t336 * t416 - t393 * t362 + t348 - t362;
t308 = t332 * t369 + t335 * t417 - t393 * t361 + t347 - t361;
t307 = (t324 * t343 + t325 * t346) * t328;
t306 = (t322 * t342 + t323 * t345) * t327;
t305 = (t320 * t341 + t321 * t344) * t326;
t304 = (t316 * t343 + t319 * t346) * t328;
t303 = (t315 * t342 + t318 * t345) * t327;
t302 = (t314 * t341 + t317 * t344) * t326;
t301 = (t310 * t343 + t313 * t346) * t328;
t300 = (t309 * t342 + t312 * t345) * t327;
t299 = (t308 * t341 + t311 * t344) * t326;
t298 = (t325 * t364 + (-t313 * t380 + t319) * m(3)) * t328;
t297 = (t323 * t364 + (-t312 * t379 + t318) * m(3)) * t327;
t296 = (t321 * t364 + (-t311 * t378 + t317) * m(3)) * t326;
t295 = (t324 * t364 + (-t310 * t380 + t316) * m(3)) * t328;
t294 = (t322 * t364 + (-t309 * t379 + t315) * m(3)) * t327;
t293 = (t320 * t364 + (-t308 * t378 + t314) * m(3)) * t326;
t292 = (t313 * t367 - t319 * t425 + t325 * t331) * t328;
t291 = (t312 * t367 - t318 * t426 + t323 * t330) * t327;
t290 = (t311 * t367 - t317 * t427 + t321 * t329) * t326;
t289 = (t310 * t367 - t316 * t425 + t324 * t331) * t328;
t288 = (t309 * t367 - t315 * t426 + t322 * t330) * t327;
t287 = (t308 * t367 - t314 * t427 + t320 * t329) * t326;
t286 = (t313 * t331 + t319 * t364 + t325 * t340) * t328;
t285 = (t312 * t330 + t318 * t364 + t323 * t339) * t327;
t284 = (t311 * t329 + t317 * t364 + t321 * t338) * t326;
t283 = (t310 * t331 + t316 * t364 + t324 * t340) * t328;
t282 = (t309 * t330 + t315 * t364 + t322 * t339) * t327;
t281 = (t308 * t329 + t314 * t364 + t320 * t338) * t326;
t280 = t307 * t364 + (-t301 * t380 + t304) * m(3);
t279 = t306 * t364 + (-t300 * t379 + t303) * m(3);
t278 = t305 * t364 + (-t299 * t378 + t302) * m(3);
t277 = t301 * t367 - t304 * t425 + t307 * t331;
t276 = t300 * t367 - t303 * t426 + t306 * t330;
t275 = t299 * t367 - t302 * t427 + t305 * t329;
t274 = t301 * t331 + t304 * t364 + t307 * t340;
t273 = t300 * t330 + t303 * t364 + t306 * t339;
t272 = t299 * t329 + t302 * t364 + t305 * t338;
t1 = [m(4) + (t283 * t324 + t289 * t310 + t295 * t316) * t328 + (t282 * t322 + t288 * t309 + t294 * t315) * t327 + (t281 * t320 + t287 * t308 + t293 * t314) * t326 (t283 * t325 + t289 * t313 + t295 * t319) * t328 + (t282 * t323 + t288 * t312 + t294 * t318) * t327 + (t281 * t321 + t287 * t311 + t293 * t317) * t326, t281 * t305 + t282 * t306 + t283 * t307 + t287 * t299 + t288 * t300 + t289 * t301 + t293 * t302 + t294 * t303 + t295 * t304 + t397; (t286 * t324 + t292 * t310 + t298 * t316) * t328 + (t285 * t322 + t291 * t309 + t297 * t315) * t327 + (t284 * t320 + t290 * t308 + t296 * t314) * t326, m(4) + (t286 * t325 + t292 * t313 + t298 * t319) * t328 + (t285 * t323 + t291 * t312 + t297 * t318) * t327 + (t284 * t321 + t290 * t311 + t296 * t317) * t326, t284 * t305 + t285 * t306 + t286 * t307 + t290 * t299 + t291 * t300 + t292 * t301 + t296 * t302 + t297 * t303 + t298 * t304 + t398; (t274 * t324 + t277 * t310 + t280 * t316) * t328 + (t273 * t322 + t276 * t309 + t279 * t315) * t327 + (t272 * t320 + t275 * t308 + t278 * t314) * t326 + t397 (t274 * t325 + t277 * t313 + t280 * t319) * t328 + (t273 * t323 + t276 * t312 + t279 * t318) * t327 + (t272 * t321 + t275 * t311 + t278 * t317) * t326 + t398, t272 * t305 + t273 * t306 + t274 * t307 + t275 * t299 + t276 * t300 + t277 * t301 + t278 * t302 + t279 * t303 + t280 * t304 + Ifges(4,3);];
MX  = t1;
