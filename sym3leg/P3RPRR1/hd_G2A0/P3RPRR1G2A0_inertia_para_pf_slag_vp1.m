% Calculate inertia matrix for parallel robot
% P3RPRR1G2A0
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2020-03-09 21:25
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3RPRR1G2A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(7,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3RPRR1G2A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:24:49
% EndTime: 2020-03-09 21:24:50
% DurationCPUTime: 0.73s
% Computational Cost: add. (2961->170), mult. (3048->225), div. (243->4), fcn. (1260->65), ass. (0->118)
t359 = pkin(7) + qJ(3,3);
t387 = pkin(1) * sin(t359) + sin(qJ(3,3)) * pkin(2);
t301 = 0.1e1 / t387;
t398 = t301 / 0.2e1;
t360 = pkin(7) + qJ(3,2);
t386 = pkin(1) * sin(t360) + sin(qJ(3,2)) * pkin(2);
t302 = 0.1e1 / t386;
t396 = t302 / 0.2e1;
t361 = pkin(7) + qJ(3,1);
t385 = pkin(1) * sin(t361) + sin(qJ(3,1)) * pkin(2);
t303 = 0.1e1 / t385;
t394 = t303 / 0.2e1;
t376 = 0.1e1 / pkin(3);
t427 = t376 / 0.2e1;
t364 = legFrame(1,2);
t391 = qJ(1,1) + pkin(7);
t333 = t364 + t391;
t326 = qJ(3,1) + t333;
t334 = -t364 + t391;
t327 = qJ(3,1) + t334;
t300 = -cos(t327) + cos(t326);
t426 = t300 * t394;
t363 = legFrame(2,2);
t390 = qJ(1,2) + pkin(7);
t331 = t363 + t390;
t324 = qJ(3,2) + t331;
t332 = -t363 + t390;
t325 = qJ(3,2) + t332;
t299 = -cos(t325) + cos(t324);
t425 = t299 * t396;
t362 = legFrame(3,2);
t389 = qJ(1,3) + pkin(7);
t329 = t362 + t389;
t322 = qJ(3,3) + t329;
t330 = -t362 + t389;
t323 = qJ(3,3) + t330;
t298 = -cos(t323) + cos(t322);
t424 = t298 * t398;
t297 = sin(t326) + sin(t327);
t423 = t297 * t394;
t296 = sin(t324) + sin(t325);
t422 = t296 * t396;
t295 = sin(t322) + sin(t323);
t421 = t295 * t398;
t420 = (-t385 * rSges(3,2) + (pkin(1) * cos(t361) + pkin(2) * cos(qJ(3,1))) * rSges(3,1)) * m(3);
t419 = (-t386 * rSges(3,2) + (pkin(1) * cos(t360) + cos(qJ(3,2)) * pkin(2)) * rSges(3,1)) * m(3);
t418 = (-t387 * rSges(3,2) + (pkin(1) * cos(t359) + pkin(2) * cos(qJ(3,3))) * rSges(3,1)) * m(3);
t347 = qJ(1,3) + t362;
t348 = qJ(1,3) - t362;
t283 = t295 * pkin(3) + (sin(t329) + sin(t330)) * pkin(2) + (sin(t347) + sin(t348)) * pkin(1);
t417 = t283 * t301;
t349 = qJ(1,2) + t363;
t350 = qJ(1,2) - t363;
t284 = t296 * pkin(3) + (sin(t331) + sin(t332)) * pkin(2) + (sin(t349) + sin(t350)) * pkin(1);
t416 = t284 * t302;
t351 = qJ(1,1) + t364;
t352 = qJ(1,1) - t364;
t285 = t297 * pkin(3) + (sin(t333) + sin(t334)) * pkin(2) + (sin(t351) + sin(t352)) * pkin(1);
t415 = t285 * t303;
t286 = -t298 * pkin(3) + (cos(t330) - cos(t329)) * pkin(2) + (cos(t348) - cos(t347)) * pkin(1);
t414 = t286 * t301;
t287 = -t299 * pkin(3) + (cos(t332) - cos(t331)) * pkin(2) + (cos(t350) - cos(t349)) * pkin(1);
t413 = t287 * t302;
t288 = -t300 * pkin(3) + (cos(t334) - cos(t333)) * pkin(2) + (cos(t352) - cos(t351)) * pkin(1);
t412 = t288 * t303;
t392 = rSges(3,1) ^ 2 + rSges(3,2) ^ 2;
t321 = t392 * m(3) + Icges(3,3);
t289 = t418 + t321;
t411 = t289 * t376;
t290 = t419 + t321;
t410 = t290 * t376;
t291 = t420 + t321;
t409 = t291 * t376;
t335 = cos(qJ(1,3) + t359);
t292 = -pkin(3) * t335 - pkin(2) * cos(t389) - cos(qJ(1,3)) * pkin(1);
t408 = t292 * t301;
t336 = cos(qJ(1,2) + t360);
t293 = -pkin(3) * t336 - pkin(2) * cos(t390) - cos(qJ(1,2)) * pkin(1);
t407 = t293 * t302;
t337 = cos(qJ(1,1) + t361);
t294 = -pkin(3) * t337 - pkin(2) * cos(t391) - cos(qJ(1,1)) * pkin(1);
t406 = t294 * t303;
t399 = t301 * t335;
t397 = t302 * t336;
t395 = t303 * t337;
t393 = t321 * t376;
t353 = sin(t362);
t354 = sin(t363);
t355 = sin(t364);
t356 = cos(t362);
t357 = cos(t363);
t358 = cos(t364);
t372 = m(2) + m(3);
t388 = (t353 * t356 + t354 * t357 + t355 * t358) * t372;
t377 = pkin(1) ^ 2;
t378 = Icges(1,3) + Icges(2,3) + Icges(3,3) + (pkin(2) ^ 2 + t377 + t392) * m(3) + (rSges(1,1) ^ 2 + rSges(1,2) ^ 2) * m(1) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2 + t377) * m(2) + 0.2e1 * ((m(2) * rSges(2,1) + m(3) * pkin(2)) * cos(pkin(7)) - m(2) * sin(pkin(7)) * rSges(2,2)) * pkin(1);
t282 = t378 + 0.2e1 * t420;
t281 = t378 + 0.2e1 * t419;
t280 = t378 + 0.2e1 * t418;
t279 = (t291 * t337 + t294 * t393) * t303;
t278 = (t290 * t336 + t293 * t393) * t302;
t277 = (t289 * t335 + t292 * t393) * t301;
t276 = (t288 * t393 + t291 * t300) * t394;
t275 = (t287 * t393 + t290 * t299) * t396;
t274 = (t286 * t393 + t289 * t298) * t398;
t273 = (-t285 * t393 + t291 * t297) * t394;
t272 = (-t284 * t393 + t290 * t296) * t396;
t271 = (-t283 * t393 + t289 * t295) * t398;
t270 = (t282 * t337 + t294 * t409) * t303;
t269 = (t281 * t336 + t293 * t410) * t302;
t268 = (t280 * t335 + t292 * t411) * t301;
t267 = (t282 * t300 + t288 * t409) * t394;
t266 = (t281 * t299 + t287 * t410) * t396;
t265 = (t280 * t298 + t286 * t411) * t398;
t264 = (t282 * t297 - t285 * t409) * t394;
t263 = (t281 * t296 - t284 * t410) * t396;
t262 = (t280 * t295 - t283 * t411) * t398;
t1 = [m(4) + (t353 ^ 2 + t354 ^ 2 + t355 ^ 2) * t372 + t262 * t421 + t263 * t422 + t264 * t423 + (-t271 * t417 - t272 * t416 - t273 * t415) * t427, t262 * t424 + t263 * t425 + t264 * t426 + (t271 * t414 + t272 * t413 + t273 * t412) * t427 + t388, t262 * t399 + t263 * t397 + t264 * t395 + (t271 * t408 + t272 * t407 + t273 * t406) * t376; t265 * t421 + t266 * t422 + t267 * t423 + (-t274 * t417 - t275 * t416 - t276 * t415) * t427 + t388, m(4) + (t356 ^ 2 + t357 ^ 2 + t358 ^ 2) * t372 + t265 * t424 + t266 * t425 + t267 * t426 + (t274 * t414 + t275 * t413 + t276 * t412) * t427, t265 * t399 + t266 * t397 + t267 * t395 + (t274 * t408 + t275 * t407 + t276 * t406) * t376; t268 * t421 + t269 * t422 + t270 * t423 + (-t277 * t417 - t278 * t416 - t279 * t415) * t427, t268 * t424 + t269 * t425 + t270 * t426 + (t277 * t414 + t278 * t413 + t279 * t412) * t427, t268 * t399 + t269 * t397 + t270 * t395 + m(4) + (t277 * t408 + t278 * t407 + t279 * t406) * t376;];
MX  = t1;
