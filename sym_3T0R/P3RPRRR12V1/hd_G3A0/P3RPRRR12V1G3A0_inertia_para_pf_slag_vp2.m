% Calculate inertia matrix for parallel robot
% P3RPRRR12V1G3A0
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
% Datum: 2020-08-06 18:28
% Revision: 8f4ff0ee124033641e65b154ac60823cef59ef1f (2020-07-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2(xP, qJ, legFrame, ...
  koppelP, pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(6,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: Ifges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRRR12V1G3A0_inertia_para_pf_slag_vp2: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-08-06 18:28:07
% EndTime: 2020-08-06 18:28:08
% DurationCPUTime: 0.80s
% Computational Cost: add. (1929->187), mult. (2451->362), div. (396->7), fcn. (1704->18), ass. (0->141)
t404 = 2 * mrSges(3,1);
t403 = 2 * mrSges(3,2);
t402 = 2 * mrSges(2,3);
t401 = -2 * Ifges(3,4);
t347 = (pkin(1) + pkin(5));
t400 = (mrSges(3,3) - mrSges(2,2));
t332 = legFrame(3,2);
t316 = sin(t332);
t319 = cos(t332);
t341 = cos(qJ(3,3));
t328 = t341 ^ 2;
t335 = sin(qJ(3,3));
t336 = sin(qJ(1,3));
t323 = pkin(6) + t347;
t342 = cos(qJ(1,3));
t352 = qJ(2,3) * t336 + t323 * t342;
t363 = t341 * qJ(2,3);
t366 = t335 * t341;
t285 = t352 * t319 * t335 + t316 * t363 + (t316 * t366 + (-t328 + 0.1e1) * t319 * t336) * pkin(3);
t325 = 0.1e1 / t335;
t399 = t285 * t325;
t333 = legFrame(2,2);
t317 = sin(t333);
t320 = cos(t333);
t343 = cos(qJ(3,2));
t329 = t343 ^ 2;
t337 = sin(qJ(3,2));
t338 = sin(qJ(1,2));
t344 = cos(qJ(1,2));
t353 = qJ(2,2) * t338 + t323 * t344;
t362 = t343 * qJ(2,2);
t365 = t337 * t343;
t286 = t353 * t320 * t337 + t317 * t362 + (t317 * t365 + (-t329 + 0.1e1) * t320 * t338) * pkin(3);
t326 = 0.1e1 / t337;
t398 = t286 * t326;
t334 = legFrame(1,2);
t318 = sin(t334);
t321 = cos(t334);
t345 = cos(qJ(3,1));
t330 = t345 ^ 2;
t339 = sin(qJ(3,1));
t340 = sin(qJ(1,1));
t346 = cos(qJ(1,1));
t354 = qJ(2,1) * t340 + t323 * t346;
t361 = t345 * qJ(2,1);
t364 = t339 * t345;
t287 = t354 * t321 * t339 + t318 * t361 + (t318 * t364 + (-t330 + 0.1e1) * t321 * t340) * pkin(3);
t327 = 0.1e1 / t339;
t397 = t287 * t327;
t288 = (t319 * pkin(3) * t341 - t352 * t316) * t335 + t336 * pkin(3) * (t341 - 0.1e1) * (t341 + 0.1e1) * t316 + t319 * t363;
t396 = t288 * t325;
t289 = (t320 * pkin(3) * t343 - t353 * t317) * t337 + t338 * pkin(3) * (t343 - 0.1e1) * (t343 + 0.1e1) * t317 + t320 * t362;
t395 = t289 * t326;
t290 = (t321 * pkin(3) * t345 - t354 * t318) * t339 + t340 * pkin(3) * (t345 - 0.1e1) * (t345 + 0.1e1) * t318 + t321 * t361;
t394 = t290 * t327;
t311 = t335 * pkin(3) + qJ(2,3);
t300 = t311 * t342 - t323 * t336;
t306 = -m(2) * pkin(1) - t347 * m(3) - t400;
t308 = 0.1e1 / t311;
t348 = m(2) + m(3);
t294 = (t300 * t348 - t306 * t336) * t308;
t393 = t294 * t325;
t312 = t337 * pkin(3) + qJ(2,2);
t301 = t312 * t344 - t323 * t338;
t309 = 0.1e1 / t312;
t295 = (t301 * t348 - t306 * t338) * t309;
t392 = t295 * t326;
t313 = t339 * pkin(3) + qJ(2,1);
t302 = t313 * t346 - t323 * t340;
t310 = 0.1e1 / t313;
t296 = (t302 * t348 - t306 * t340) * t310;
t391 = t296 * t327;
t303 = t341 * mrSges(3,1) - t335 * mrSges(3,2);
t390 = t303 * t325;
t304 = t343 * mrSges(3,1) - t337 * mrSges(3,2);
t389 = t304 * t326;
t305 = t345 * mrSges(3,1) - t339 * mrSges(3,2);
t388 = t305 * t327;
t387 = t306 * t325;
t386 = t306 * t326;
t385 = t306 * t327;
t384 = t316 * t325;
t383 = t316 * t342;
t382 = t317 * t326;
t381 = t317 * t344;
t380 = t318 * t327;
t379 = t318 * t346;
t378 = t319 * t325;
t377 = t319 * t342;
t376 = t320 * t326;
t375 = t320 * t344;
t374 = t321 * t327;
t373 = t321 * t346;
t372 = t325 * t348;
t350 = 0.1e1 / pkin(3);
t371 = t325 * t350;
t370 = t326 * t348;
t369 = t326 * t350;
t368 = t327 * t348;
t367 = t327 * t350;
t360 = t316 * t371;
t359 = t317 * t369;
t358 = t318 * t367;
t357 = t319 * t371;
t356 = t320 * t369;
t355 = t321 * t367;
t351 = Ifges(2,1) + Ifges(3,2) + Ifges(1,3) + (2 * (m(3) * pkin(5) + t400) * pkin(1)) + t348 * (pkin(1) ^ 2) + (2 * mrSges(3,3) * pkin(5)) + (m(3) * pkin(5) ^ 2);
t331 = Ifges(3,1) - Ifges(3,2);
t315 = -t347 * mrSges(3,1) + Ifges(3,5);
t314 = t347 * mrSges(3,2) - Ifges(3,6);
t299 = t339 * t314 + t315 * t345;
t298 = t337 * t314 + t315 * t343;
t297 = t335 * t314 + t315 * t341;
t293 = t364 * t401 + t331 * t330 + (t348 * qJ(2,1) + t339 * t404 + t345 * t403 + t402) * qJ(2,1) + t351;
t292 = t365 * t401 + t331 * t329 + (t348 * qJ(2,2) + t337 * t404 + t343 * t403 + t402) * qJ(2,2) + t351;
t291 = t366 * t401 + t331 * t328 + (t348 * qJ(2,3) + t335 * t404 + t341 * t403 + t402) * qJ(2,3) + t351;
t284 = (-t299 * t340 + t302 * t305) * t310;
t283 = (-t298 * t338 + t301 * t304) * t309;
t282 = (-t297 * t336 + t300 * t303) * t308;
t281 = (-t293 * t340 + t302 * t306) * t310;
t280 = (-t292 * t338 + t301 * t306) * t309;
t279 = (-t291 * t336 + t300 * t306) * t308;
t278 = -t305 * t355 + (t290 * t368 - t306 * t379) * t310;
t277 = -t304 * t356 + (t289 * t370 - t306 * t381) * t309;
t276 = -t303 * t357 + (t288 * t372 - t306 * t383) * t308;
t275 = -t305 * t358 + (t287 * t368 + t306 * t373) * t310;
t274 = -t304 * t359 + (t286 * t370 + t306 * t375) * t309;
t273 = -t303 * t360 + (t285 * t372 + t306 * t377) * t308;
t272 = -Ifges(3,3) * t355 + (t290 * t388 - t299 * t379) * t310;
t271 = -Ifges(3,3) * t356 + (t289 * t389 - t298 * t381) * t309;
t270 = -Ifges(3,3) * t357 + (t288 * t390 - t297 * t383) * t308;
t269 = -Ifges(3,3) * t358 + (t287 * t388 + t299 * t373) * t310;
t268 = -Ifges(3,3) * t359 + (t286 * t389 + t298 * t375) * t309;
t267 = -Ifges(3,3) * t360 + (t285 * t390 + t297 * t377) * t308;
t266 = -t299 * t355 + (t290 * t385 - t293 * t379) * t310;
t265 = -t298 * t356 + (t289 * t386 - t292 * t381) * t309;
t264 = -t297 * t357 + (t288 * t387 - t291 * t383) * t308;
t263 = -t299 * t358 + (t287 * t385 + t293 * t373) * t310;
t262 = -t298 * t359 + (t286 * t386 + t292 * t375) * t309;
t261 = -t297 * t360 + (t285 * t387 + t291 * t377) * t308;
t1 = [m(4) + (t263 * t373 + t275 * t397) * t310 + (t262 * t375 + t274 * t398) * t309 + (t261 * t377 + t273 * t399) * t308 + (-t267 * t384 - t268 * t382 - t269 * t380) * t350, (-t263 * t379 + t275 * t394) * t310 + (-t262 * t381 + t274 * t395) * t309 + (-t261 * t383 + t273 * t396) * t308 + (-t267 * t378 - t268 * t376 - t269 * t374) * t350, (-t263 * t340 + t275 * t302) * t310 + (-t262 * t338 + t274 * t301) * t309 + (-t261 * t336 + t273 * t300) * t308; (t266 * t373 + t278 * t397) * t310 + (t265 * t375 + t277 * t398) * t309 + (t264 * t377 + t276 * t399) * t308 + (-t270 * t384 - t271 * t382 - t272 * t380) * t350, m(4) + (-t266 * t379 + t278 * t394) * t310 + (-t265 * t381 + t277 * t395) * t309 + (-t264 * t383 + t276 * t396) * t308 + (-t270 * t378 - t271 * t376 - t272 * t374) * t350, (-t266 * t340 + t278 * t302) * t310 + (-t265 * t338 + t277 * t301) * t309 + (-t264 * t336 + t276 * t300) * t308; (t281 * t373 + t287 * t391) * t310 + (t280 * t375 + t286 * t392) * t309 + (t279 * t377 + t285 * t393) * t308 + (-t282 * t384 - t283 * t382 - t284 * t380) * t350, (-t281 * t379 + t290 * t391) * t310 + (-t280 * t381 + t289 * t392) * t309 + (-t279 * t383 + t288 * t393) * t308 + (-t282 * t378 - t283 * t376 - t284 * t374) * t350, m(4) + (-t281 * t340 + t296 * t302) * t310 + (-t280 * t338 + t295 * t301) * t309 + (-t279 * t336 + t294 * t300) * t308;];
MX  = t1;
