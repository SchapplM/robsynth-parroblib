% Calculate inertia matrix for parallel robot
% P3RPRRR12V1G2A0
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4]';
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
% Datum: 2020-08-06 18:25
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G2A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:24:35
% EndTime: 2020-08-06 18:24:36
% DurationCPUTime: 0.73s
% Computational Cost: add. (1929->187), mult. (2418->359), div. (396->7), fcn. (1671->18), ass. (0->144)
t409 = 2 * mrSges(3,1);
t408 = 2 * mrSges(3,2);
t407 = 2 * mrSges(2,3);
t406 = -2 * Ifges(3,4);
t352 = (pkin(1) + pkin(5));
t346 = cos(qJ(3,3));
t405 = pkin(3) * t346;
t347 = cos(qJ(1,3));
t404 = pkin(3) * t347;
t348 = cos(qJ(3,2));
t403 = pkin(3) * t348;
t349 = cos(qJ(1,2));
t402 = pkin(3) * t349;
t350 = cos(qJ(3,1));
t401 = pkin(3) * t350;
t351 = cos(qJ(1,1));
t400 = pkin(3) * t351;
t399 = (mrSges(3,3) - mrSges(2,2));
t328 = pkin(6) + t352;
t341 = sin(qJ(1,3));
t305 = qJ(2,3) * t347 - t328 * t341;
t337 = legFrame(3,2);
t321 = sin(t337);
t324 = cos(t337);
t340 = sin(qJ(3,3));
t365 = t346 * qJ(2,3);
t287 = (-t305 * t324 + t321 * t405) * t340 + (t346 - 0.1e1) * (t346 + 0.1e1) * t324 * t404 + t321 * t365;
t330 = 0.1e1 / t340;
t398 = t287 * t330;
t343 = sin(qJ(1,2));
t306 = qJ(2,2) * t349 - t328 * t343;
t338 = legFrame(2,2);
t322 = sin(t338);
t325 = cos(t338);
t342 = sin(qJ(3,2));
t364 = t348 * qJ(2,2);
t288 = (-t306 * t325 + t322 * t403) * t342 + (t348 - 0.1e1) * (t348 + 0.1e1) * t325 * t402 + t322 * t364;
t331 = 0.1e1 / t342;
t397 = t288 * t331;
t345 = sin(qJ(1,1));
t307 = qJ(2,1) * t351 - t328 * t345;
t339 = legFrame(1,2);
t323 = sin(t339);
t326 = cos(t339);
t344 = sin(qJ(3,1));
t363 = t350 * qJ(2,1);
t289 = (-t307 * t326 + t323 * t401) * t344 + (t350 - 0.1e1) * (t350 + 0.1e1) * t326 * t400 + t323 * t363;
t332 = 0.1e1 / t344;
t396 = t289 * t332;
t333 = t346 ^ 2;
t293 = (t305 * t321 + t324 * t405) * t340 + (-t333 + 0.1e1) * t321 * t404 + t324 * t365;
t395 = t293 * t330;
t334 = t348 ^ 2;
t294 = (t306 * t322 + t325 * t403) * t342 + (-t334 + 0.1e1) * t322 * t402 + t325 * t364;
t394 = t294 * t331;
t335 = t350 ^ 2;
t295 = (t307 * t323 + t326 * t401) * t344 + (-t335 + 0.1e1) * t323 * t400 + t326 * t363;
t393 = t295 * t332;
t316 = t340 * pkin(3) + qJ(2,3);
t302 = t341 * t316 + t328 * t347;
t311 = -m(2) * pkin(1) - t352 * m(3) - t399;
t313 = 0.1e1 / t316;
t353 = m(2) + m(3);
t296 = (t302 * t353 + t311 * t347) * t313;
t392 = t296 * t330;
t317 = t342 * pkin(3) + qJ(2,2);
t303 = t343 * t317 + t328 * t349;
t314 = 0.1e1 / t317;
t297 = (t303 * t353 + t311 * t349) * t314;
t391 = t297 * t331;
t318 = t344 * pkin(3) + qJ(2,1);
t304 = t345 * t318 + t328 * t351;
t315 = 0.1e1 / t318;
t298 = (t304 * t353 + t311 * t351) * t315;
t390 = t298 * t332;
t308 = t346 * mrSges(3,1) - t340 * mrSges(3,2);
t389 = t308 * t330;
t309 = t348 * mrSges(3,1) - t342 * mrSges(3,2);
t388 = t309 * t331;
t310 = t350 * mrSges(3,1) - t344 * mrSges(3,2);
t387 = t310 * t332;
t386 = t311 * t330;
t385 = t311 * t331;
t384 = t311 * t332;
t383 = t321 * t330;
t382 = t321 * t341;
t381 = t322 * t331;
t380 = t322 * t343;
t379 = t323 * t332;
t378 = t323 * t345;
t377 = t324 * t330;
t376 = t324 * t341;
t375 = t325 * t331;
t374 = t325 * t343;
t373 = t326 * t332;
t372 = t326 * t345;
t371 = t330 * t353;
t355 = 0.1e1 / pkin(3);
t370 = t330 * t355;
t369 = t331 * t353;
t368 = t331 * t355;
t367 = t332 * t353;
t366 = t332 * t355;
t362 = t321 * t370;
t361 = t322 * t368;
t360 = t323 * t366;
t359 = t324 * t370;
t358 = t325 * t368;
t357 = t326 * t366;
t356 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + (2 * (m(3) * pkin(5) + t399) * pkin(1)) + t353 * (pkin(1) ^ 2) + (2 * mrSges(3,3) * pkin(5)) + (m(3) * pkin(5) ^ 2);
t336 = Ifges(3,1) - Ifges(3,2);
t320 = -t352 * mrSges(3,1) + Ifges(3,5);
t319 = t352 * mrSges(3,2) - Ifges(3,6);
t301 = t319 * t344 + t320 * t350;
t300 = t319 * t342 + t320 * t348;
t299 = t319 * t340 + t320 * t346;
t292 = t344 * t350 * t406 + t336 * t335 + (t353 * qJ(2,1) + t344 * t409 + t350 * t408 + t407) * qJ(2,1) + t356;
t291 = t342 * t348 * t406 + t336 * t334 + (t353 * qJ(2,2) + t342 * t409 + t348 * t408 + t407) * qJ(2,2) + t356;
t290 = t340 * t346 * t406 + t336 * t333 + (t353 * qJ(2,3) + t340 * t409 + t346 * t408 + t407) * qJ(2,3) + t356;
t286 = (t301 * t351 + t304 * t310) * t315;
t285 = (t300 * t349 + t303 * t309) * t314;
t284 = (t299 * t347 + t302 * t308) * t313;
t283 = (t292 * t351 + t304 * t311) * t315;
t282 = (t291 * t349 + t303 * t311) * t314;
t281 = (t290 * t347 + t302 * t311) * t313;
t280 = -t310 * t357 + (t295 * t367 - t311 * t378) * t315;
t279 = -t309 * t358 + (t294 * t369 - t311 * t380) * t314;
t278 = -t308 * t359 + (t293 * t371 - t311 * t382) * t313;
t277 = -t310 * t360 + (t289 * t367 + t311 * t372) * t315;
t276 = -t309 * t361 + (t288 * t369 + t311 * t374) * t314;
t275 = -t308 * t362 + (t287 * t371 + t311 * t376) * t313;
t274 = -Ifges(3,3) * t357 + (t295 * t387 - t301 * t378) * t315;
t273 = -Ifges(3,3) * t358 + (t294 * t388 - t300 * t380) * t314;
t272 = -Ifges(3,3) * t359 + (t293 * t389 - t299 * t382) * t313;
t271 = -Ifges(3,3) * t360 + (t289 * t387 + t301 * t372) * t315;
t270 = -Ifges(3,3) * t361 + (t288 * t388 + t300 * t374) * t314;
t269 = -Ifges(3,3) * t362 + (t287 * t389 + t299 * t376) * t313;
t268 = -t301 * t357 + (-t292 * t378 + t295 * t384) * t315;
t267 = -t300 * t358 + (-t291 * t380 + t294 * t385) * t314;
t266 = -t299 * t359 + (-t290 * t382 + t293 * t386) * t313;
t265 = -t301 * t360 + (t289 * t384 + t292 * t372) * t315;
t264 = -t300 * t361 + (t288 * t385 + t291 * t374) * t314;
t263 = -t299 * t362 + (t287 * t386 + t290 * t376) * t313;
t1 = [m(4) + (t265 * t372 + t277 * t396) * t315 + (t264 * t374 + t276 * t397) * t314 + (t263 * t376 + t275 * t398) * t313 + (-t269 * t383 - t270 * t381 - t271 * t379) * t355, (-t265 * t378 + t277 * t393) * t315 + (-t264 * t380 + t276 * t394) * t314 + (-t263 * t382 + t275 * t395) * t313 + (-t269 * t377 - t270 * t375 - t271 * t373) * t355, (t265 * t351 + t277 * t304) * t315 + (t264 * t349 + t276 * t303) * t314 + (t263 * t347 + t275 * t302) * t313; (t268 * t372 + t280 * t396) * t315 + (t267 * t374 + t279 * t397) * t314 + (t266 * t376 + t278 * t398) * t313 + (-t272 * t383 - t273 * t381 - t274 * t379) * t355, m(4) + (-t268 * t378 + t280 * t393) * t315 + (-t267 * t380 + t279 * t394) * t314 + (-t266 * t382 + t278 * t395) * t313 + (-t272 * t377 - t273 * t375 - t274 * t373) * t355, (t268 * t351 + t280 * t304) * t315 + (t267 * t349 + t279 * t303) * t314 + (t266 * t347 + t278 * t302) * t313; (t283 * t372 + t289 * t390) * t315 + (t282 * t374 + t288 * t391) * t314 + (t281 * t376 + t287 * t392) * t313 + (-t284 * t383 - t285 * t381 - t286 * t379) * t355, (-t283 * t378 + t295 * t390) * t315 + (-t282 * t380 + t294 * t391) * t314 + (-t281 * t382 + t293 * t392) * t313 + (-t284 * t377 - t285 * t375 - t286 * t373) * t355, m(4) + (t283 * t351 + t298 * t304) * t315 + (t282 * t349 + t297 * t303) * t314 + (t281 * t347 + t296 * t302) * t313;];
MX  = t1;
