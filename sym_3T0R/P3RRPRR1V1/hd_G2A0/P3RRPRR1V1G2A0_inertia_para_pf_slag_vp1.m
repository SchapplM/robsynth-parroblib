% Calculate inertia matrix for parallel robot
% P3RRPRR1V1G2A0
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
%   pkin=[a3,a4,d4]';
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
% Datum: 2022-11-04 17:08
% Revision: e482436b586c4f286726c907c195760c5ac72455 (2022-11-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(3,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: pkin has to be [3x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RRPRR1V1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2022-11-04 17:07:58
% EndTime: 2022-11-04 17:07:59
% DurationCPUTime: 0.68s
% Computational Cost: add. (2184->174), mult. (2649->325), div. (342->7), fcn. (1164->24), ass. (0->150)
t418 = m(2) / 0.2e1;
t417 = -Icges(2,1) / 0.2e1 - Icges(3,1) / 0.2e1;
t416 = Icges(2,2) / 0.2e1 + Icges(3,2) / 0.2e1;
t415 = m(3) / 0.2e1;
t414 = m(3) * rSges(3,2);
t413 = rSges(2,3) * m(2);
t346 = pkin(1) + rSges(3,1);
t356 = rSges(2,2) ^ 2;
t358 = rSges(2,1) ^ 2;
t412 = -(rSges(3,2) + t346) * (rSges(3,2) - t346) * m(3) / 0.2e1 + (-t356 + t358) * t418 + t416 + t417;
t334 = sin(qJ(2,3));
t340 = cos(qJ(2,3));
t305 = -rSges(3,2) * t334 + t340 * t346;
t411 = m(3) * t305;
t336 = sin(qJ(2,2));
t342 = cos(qJ(2,2));
t306 = -rSges(3,2) * t336 + t342 * t346;
t410 = m(3) * t306;
t338 = sin(qJ(2,1));
t344 = cos(qJ(2,1));
t307 = -rSges(3,2) * t338 + t344 * t346;
t409 = m(3) * t307;
t328 = pkin(3) + qJ(3,3);
t317 = 0.1e1 / t328;
t408 = m(3) * t317;
t329 = pkin(3) + qJ(3,2);
t318 = 0.1e1 / t329;
t407 = m(3) * t318;
t330 = pkin(3) + qJ(3,1);
t319 = 0.1e1 / t330;
t406 = m(3) * t319;
t405 = m(3) * t346;
t325 = rSges(3,3) + qJ(3,3);
t364 = -rSges(2,1) * t413 + Icges(2,5) + Icges(3,5);
t366 = rSges(2,2) * t413 - Icges(2,6) - Icges(3,6);
t284 = (-t325 * t405 + t364) * t334 - t340 * (t325 * t414 + t366);
t404 = t284 * t317;
t351 = pkin(2) + pkin(1);
t324 = 0.1e1 / t351;
t403 = t284 * t324;
t326 = rSges(3,3) + qJ(3,2);
t285 = (-t326 * t405 + t364) * t336 - t342 * (t326 * t414 + t366);
t402 = t285 * t318;
t401 = t285 * t324;
t327 = rSges(3,3) + qJ(3,1);
t286 = (-t327 * t405 + t364) * t338 - t344 * (t327 * t414 + t366);
t400 = t286 * t319;
t399 = t286 * t324;
t365 = rSges(3,2) ^ 2 + (pkin(1) ^ 2) + ((2 * pkin(1) + rSges(3,1)) * rSges(3,1));
t373 = t356 + t358;
t398 = (m(2) * t373 + m(3) * t365 + Icges(2,3) + Icges(3,3)) * t324;
t331 = legFrame(3,2);
t310 = sin(t331);
t313 = cos(t331);
t335 = sin(qJ(1,3));
t381 = t335 * t340;
t295 = -t310 * t381 + t313 * t334;
t321 = 0.1e1 / t340;
t397 = t295 * t321;
t332 = legFrame(2,2);
t311 = sin(t332);
t314 = cos(t332);
t337 = sin(qJ(1,2));
t379 = t337 * t342;
t296 = -t311 * t379 + t314 * t336;
t322 = 0.1e1 / t342;
t396 = t296 * t322;
t333 = legFrame(1,2);
t312 = sin(t333);
t315 = cos(t333);
t339 = sin(qJ(1,1));
t377 = t339 * t344;
t297 = -t312 * t377 + t315 * t338;
t323 = 0.1e1 / t344;
t395 = t297 * t323;
t298 = t310 * t334 + t313 * t381;
t394 = t298 * t321;
t299 = t311 * t336 + t314 * t379;
t393 = t299 * t322;
t300 = t312 * t338 + t315 * t377;
t392 = t300 * t323;
t391 = t305 * t321;
t390 = t306 * t322;
t389 = t307 * t323;
t388 = t310 * t321;
t387 = t311 * t322;
t386 = t312 * t323;
t385 = t313 * t321;
t384 = t314 * t322;
t383 = t315 * t323;
t382 = t334 * t351;
t380 = t336 * t351;
t378 = t338 * t351;
t376 = t340 * t351;
t375 = t342 * t351;
t374 = t344 * t351;
t372 = t321 * t403;
t341 = cos(qJ(1,3));
t371 = t341 * t403;
t370 = t322 * t401;
t343 = cos(qJ(1,2));
t369 = t343 * t401;
t368 = t323 * t399;
t345 = cos(qJ(1,1));
t367 = t345 * t399;
t363 = -t328 * t341 + t335 * t376;
t362 = -t329 * t343 + t337 * t375;
t361 = -t330 * t345 + t339 * t374;
t360 = Icges(1,3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (0.2e1 * rSges(2,3) ^ 2 + t373) * t418 + t416 - t417;
t354 = 0.2e1 * qJ(2,1);
t353 = 0.2e1 * qJ(2,2);
t352 = 0.2e1 * qJ(2,3);
t304 = -m(2) * rSges(2,1) * rSges(2,2) - rSges(3,2) * t405 + Icges(2,4) + Icges(3,4);
t303 = t330 * t339 + t345 * t374;
t302 = t329 * t337 + t343 * t375;
t301 = t328 * t335 + t341 * t376;
t292 = t312 * t378 + t315 * t361;
t291 = t311 * t380 + t314 * t362;
t290 = t310 * t382 + t313 * t363;
t289 = -t312 * t361 + t315 * t378;
t288 = -t311 * t362 + t314 * t380;
t287 = -t310 * t363 + t313 * t382;
t283 = (-t307 * t345 + t303) * t406;
t282 = (-t306 * t343 + t302) * t407;
t281 = (-t305 * t341 + t301) * t408;
t280 = (-t300 * t389 + t292) * t406;
t279 = (-t299 * t390 + t291) * t407;
t278 = (-t298 * t391 + t290) * t408;
t277 = (-t297 * t389 + t289) * t406;
t276 = (-t296 * t390 + t288) * t407;
t275 = (-t295 * t391 + t287) * t408;
t274 = cos(t354) * t412 + t304 * sin(t354) + (0.2e1 * t327 ^ 2 + t365) * t415 + t360;
t273 = cos(t353) * t412 + t304 * sin(t353) + (0.2e1 * t326 ^ 2 + t365) * t415 + t360;
t272 = cos(t352) * t412 + t304 * sin(t352) + (0.2e1 * t325 ^ 2 + t365) * t415 + t360;
t271 = (t300 * t400 + t312 * t398) * t323;
t270 = (t299 * t402 + t311 * t398) * t322;
t269 = (t298 * t404 + t310 * t398) * t321;
t268 = (t297 * t400 + t315 * t398) * t323;
t267 = (t296 * t402 + t314 * t398) * t322;
t266 = (t295 * t404 + t313 * t398) * t321;
t265 = (t274 * t345 - t303 * t409) * t319;
t264 = (t273 * t343 - t302 * t410) * t318;
t263 = (t272 * t341 - t301 * t411) * t317;
t262 = t312 * t368 + (t274 * t392 - t292 * t409) * t319;
t261 = t311 * t370 + (t273 * t393 - t291 * t410) * t318;
t260 = t310 * t372 + (t272 * t394 - t290 * t411) * t317;
t259 = t315 * t368 + (t274 * t395 - t289 * t409) * t319;
t258 = t314 * t370 + (t273 * t396 - t288 * t410) * t318;
t257 = t313 * t372 + (t272 * t397 - t287 * t411) * t317;
t1 = [m(4) + (t262 * t392 + t280 * t292) * t319 + (t261 * t393 + t279 * t291) * t318 + (t260 * t394 + t278 * t290) * t317 + (t269 * t388 + t270 * t387 + t271 * t386) * t324, (t262 * t395 + t280 * t289) * t319 + (t261 * t396 + t279 * t288) * t318 + (t260 * t397 + t278 * t287) * t317 + (t269 * t385 + t270 * t384 + t271 * t383) * t324, (t262 * t345 + t280 * t303) * t319 + (t261 * t343 + t279 * t302) * t318 + (t260 * t341 + t278 * t301) * t317; (t259 * t392 + t277 * t292) * t319 + (t258 * t393 + t276 * t291) * t318 + (t257 * t394 + t275 * t290) * t317 + (t266 * t388 + t267 * t387 + t268 * t386) * t324, m(4) + (t259 * t395 + t277 * t289) * t319 + (t258 * t396 + t276 * t288) * t318 + (t257 * t397 + t275 * t287) * t317 + (t266 * t385 + t267 * t384 + t268 * t383) * t324, (t259 * t345 + t277 * t303) * t319 + (t258 * t343 + t276 * t302) * t318 + (t257 * t341 + t275 * t301) * t317; (t283 * t292 + (t265 * t300 + t312 * t367) * t323) * t319 + (t282 * t291 + (t264 * t299 + t311 * t369) * t322) * t318 + (t281 * t290 + (t263 * t298 + t310 * t371) * t321) * t317, (t283 * t289 + (t265 * t297 + t315 * t367) * t323) * t319 + (t282 * t288 + (t264 * t296 + t314 * t369) * t322) * t318 + (t281 * t287 + (t263 * t295 + t313 * t371) * t321) * t317, m(4) + (t265 * t345 + t283 * t303) * t319 + (t264 * t343 + t282 * t302) * t318 + (t263 * t341 + t281 * t301) * t317;];
MX  = t1;
