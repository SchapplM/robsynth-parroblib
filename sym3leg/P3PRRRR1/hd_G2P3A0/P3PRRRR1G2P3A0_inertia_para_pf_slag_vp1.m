% Calculate inertia matrix for parallel robot
% P3PRRRR1G2P3A0
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
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2020-03-09 21:16
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G2P3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:16:23
% EndTime: 2020-03-09 21:16:24
% DurationCPUTime: 0.72s
% Computational Cost: add. (1041->136), mult. (2220->270), div. (576->10), fcn. (1470->24), ass. (0->134)
t395 = m(3) / 0.2e1;
t394 = Icges(3,2) / 0.2e1;
t393 = rSges(3,3) * m(3);
t333 = rSges(3,2) ^ 2;
t334 = rSges(3,1) ^ 2;
t392 = (-t333 + t334) * t395 - Icges(3,1) / 0.2e1 + t394;
t315 = sin(qJ(3,3));
t321 = cos(qJ(3,3));
t391 = m(3) * (t315 * rSges(3,1) + t321 * rSges(3,2));
t317 = sin(qJ(3,2));
t323 = cos(qJ(3,2));
t390 = m(3) * (t317 * rSges(3,1) + t323 * rSges(3,2));
t319 = sin(qJ(3,1));
t325 = cos(qJ(3,1));
t389 = m(3) * (t319 * rSges(3,1) + t325 * rSges(3,2));
t312 = legFrame(3,2);
t296 = sin(t312);
t335 = 0.1e1 / pkin(2);
t388 = t296 * t335;
t313 = legFrame(2,2);
t297 = sin(t313);
t387 = t297 * t335;
t314 = legFrame(1,2);
t298 = sin(t314);
t386 = t298 * t335;
t299 = cos(t312);
t385 = t299 * t335;
t300 = cos(t313);
t384 = t300 * t335;
t301 = cos(t314);
t383 = t301 * t335;
t316 = sin(qJ(2,3));
t303 = 0.1e1 / t316;
t306 = 0.1e1 / t321;
t382 = t303 * t306;
t322 = cos(qJ(2,3));
t381 = t303 * t322;
t318 = sin(qJ(2,2));
t304 = 0.1e1 / t318;
t308 = 0.1e1 / t323;
t380 = t304 * t308;
t324 = cos(qJ(2,2));
t379 = t304 * t324;
t320 = sin(qJ(2,1));
t305 = 0.1e1 / t320;
t310 = 0.1e1 / t325;
t378 = t305 * t310;
t326 = cos(qJ(2,1));
t377 = t305 * t326;
t376 = t306 * t335;
t375 = t308 * t335;
t374 = t310 * t335;
t373 = t316 * t321;
t372 = t318 * t323;
t371 = t320 * t325;
t370 = t333 + t334;
t369 = t306 * t391;
t368 = t308 * t390;
t367 = t310 * t389;
t279 = t296 * t373 - t299 * t315;
t366 = t279 * t382;
t280 = t297 * t372 - t300 * t317;
t365 = t280 * t380;
t281 = t298 * t371 - t301 * t319;
t364 = t281 * t378;
t282 = t296 * t315 + t299 * t373;
t363 = t282 * t382;
t283 = t297 * t317 + t300 * t372;
t362 = t283 * t380;
t284 = t298 * t319 + t301 * t371;
t361 = t284 * t378;
t302 = m(1) + m(2) + m(3);
t360 = t302 * t382;
t359 = t302 * t380;
t358 = t302 * t378;
t357 = 0.1e1 / t321 ^ 2 * t315 * t381;
t356 = 0.1e1 / t323 ^ 2 * t317 * t379;
t355 = 0.1e1 / t325 ^ 2 * t319 * t377;
t354 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t370) * t395 + t394 + Icges(3,1) / 0.2e1;
t292 = m(2) * rSges(2,2) - t393;
t329 = m(2) * rSges(2,1);
t273 = (t329 + (rSges(3,1) * t321 - rSges(3,2) * t315) * m(3)) * t322 - t316 * t292;
t295 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t330 = 0.2e1 * qJ(3,3);
t264 = cos(t330) * t392 + t295 * sin(t330) + t354;
t293 = -rSges(3,2) * t393 + Icges(3,6);
t294 = rSges(3,1) * t393 - Icges(3,5);
t276 = t293 * t321 - t294 * t315;
t344 = t264 * t357 - t276 * t306;
t243 = t273 * t366 + t344 * t385;
t291 = t370 * m(3) + Icges(3,3);
t341 = t276 * t357 - t291 * t306;
t353 = t243 * t357 - (-t279 * t369 + t341 * t385) * t306;
t274 = (t329 + (rSges(3,1) * t323 - rSges(3,2) * t317) * m(3)) * t324 - t318 * t292;
t331 = 0.2e1 * qJ(3,2);
t265 = cos(t331) * t392 + t295 * sin(t331) + t354;
t277 = t293 * t323 - t294 * t317;
t343 = t265 * t356 - t277 * t308;
t244 = t274 * t365 + t343 * t384;
t340 = t277 * t356 - t291 * t308;
t352 = t244 * t356 - (-t280 * t368 + t340 * t384) * t308;
t275 = (t329 + (rSges(3,1) * t325 - rSges(3,2) * t319) * m(3)) * t326 - t320 * t292;
t332 = 0.2e1 * qJ(3,1);
t266 = cos(t332) * t392 + t295 * sin(t332) + t354;
t278 = t293 * t325 - t294 * t319;
t342 = t266 * t355 - t278 * t310;
t245 = t275 * t364 + t342 * t383;
t339 = t278 * t355 - t291 * t310;
t351 = t245 * t355 - (-t281 * t367 + t339 * t383) * t310;
t246 = t273 * t363 - t344 * t388;
t350 = t246 * t357 - (-t282 * t369 - t341 * t388) * t306;
t247 = t274 * t362 - t343 * t387;
t349 = t247 * t356 - (-t283 * t368 - t340 * t387) * t308;
t248 = t275 * t361 - t342 * t386;
t348 = t248 * t355 - (-t284 * t367 - t339 * t386) * t310;
t255 = (-t264 * t376 + t273 * t322) * t303;
t347 = t255 * t357 - (-t303 * t276 * t376 - t322 * t391) * t306;
t256 = (-t265 * t375 + t274 * t324) * t304;
t346 = t256 * t356 - (-t304 * t277 * t375 - t324 * t390) * t308;
t257 = (-t266 * t374 + t275 * t326) * t305;
t345 = t257 * t355 - (-t305 * t278 * t374 - t326 * t389) * t310;
t338 = t273 * t357 + t316 * t369;
t337 = t274 * t356 + t318 * t368;
t336 = t275 * t355 + t320 * t367;
t269 = (-t275 * t374 + t302 * t326) * t305;
t268 = (-t274 * t375 + t302 * t324) * t304;
t267 = (-t273 * t376 + t302 * t322) * t303;
t254 = t284 * t358 - t336 * t386;
t253 = t283 * t359 - t337 * t387;
t252 = t282 * t360 - t338 * t388;
t251 = t281 * t358 + t336 * t383;
t250 = t280 * t359 + t337 * t384;
t249 = t279 * t360 + t338 * t385;
t1 = [t249 * t366 + t250 * t365 + t251 * t364 + m(4) + (t353 * t299 + t352 * t300 + t351 * t301) * t335, t249 * t363 + t250 * t362 + t251 * t361 + (-t353 * t296 - t352 * t297 - t351 * t298) * t335, t249 * t381 + t250 * t379 + t251 * t377 + (-t243 * t382 - t244 * t380 - t245 * t378) * t335; t252 * t366 + t253 * t365 + t254 * t364 + (t350 * t299 + t349 * t300 + t348 * t301) * t335, t252 * t363 + t253 * t362 + t254 * t361 + m(4) + (-t350 * t296 - t349 * t297 - t348 * t298) * t335, t252 * t381 + t253 * t379 + t254 * t377 + (-t246 * t382 - t247 * t380 - t248 * t378) * t335; t267 * t366 + t268 * t365 + t269 * t364 + (t347 * t299 + t346 * t300 + t345 * t301) * t335, t267 * t363 + t268 * t362 + t269 * t361 + (-t347 * t296 - t346 * t297 - t345 * t298) * t335, t267 * t381 + t268 * t379 + t269 * t377 + m(4) + (-t255 * t382 - t256 * t380 - t257 * t378) * t335;];
MX  = t1;
