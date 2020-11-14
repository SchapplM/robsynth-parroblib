% Calculate inertia matrix for parallel robot
% P3PRRRR1G3A0
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
% Datum: 2020-03-09 21:02
% Revision: 0f11fd83bca0a8cdff505979e09e2c4d81033460 (2020-02-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MX = P3PRRRR1G3A0_inertia_para_pf_slag_vp1(xP, qJ, legFrame, ...
  koppelP, pkin, m, rSges, Icges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,3),zeros(3,3),zeros(3,3),zeros(2,1),zeros(3+1,1),zeros(3+1,3),zeros(3+1,6)}
assert(isreal(xP) && all(size(xP) == [3 1]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: xP has to be [3x1] (double)');
assert(isreal(qJ) && all(size(qJ) == [3 3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: qJ has to be [3x3] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [4 1]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: m has to be [4x1] (double)'); 
assert(isreal(rSges) && all(size(rSges) == [4,3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: rSges has to be [4x3] (double)');
assert(isreal(Icges) && all(size(Icges) == [4 6]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: Icges has to be [4x6] (double)'); 
assert(isreal(legFrame) && all(size(legFrame) == [3 3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: legFrame has to be [3x1] (double)');
assert(isreal(koppelP) && all(size(koppelP) == [3 3]), ...
  'P3PRRRR1G3A0_inertia_para_pf_slag_vp1: Koppelpunkt has to be [3x3] (double)');

%% Symbolic Calculation
% From inertia_para_plfcoord_par1_matlab.m
% OptimizationMode: 2
% StartTime: 2020-03-09 21:02:17
% EndTime: 2020-03-09 21:02:17
% DurationCPUTime: 0.55s
% Computational Cost: add. (900->136), mult. (1755->279), div. (441->10), fcn. (1062->24), ass. (0->131)
t368 = m(3) / 0.2e1;
t367 = Icges(3,2) / 0.2e1;
t300 = cos(qJ(3,3));
t285 = 0.1e1 / t300;
t366 = t285 ^ 2;
t302 = cos(qJ(3,2));
t287 = 0.1e1 / t302;
t365 = t287 ^ 2;
t304 = cos(qJ(3,1));
t289 = 0.1e1 / t304;
t364 = t289 ^ 2;
t363 = rSges(3,3) * m(3);
t312 = rSges(3,2) ^ 2;
t313 = rSges(3,1) ^ 2;
t362 = (-t312 + t313) * t368 - Icges(3,1) / 0.2e1 + t367;
t291 = legFrame(3,2);
t275 = sin(t291);
t278 = cos(t291);
t295 = sin(qJ(2,3));
t301 = cos(qJ(2,3));
t258 = -t275 * t301 + t278 * t295;
t282 = 0.1e1 / t295;
t361 = t258 * t282;
t259 = t295 * t275 + t278 * t301;
t360 = t259 * t282;
t292 = legFrame(2,2);
t276 = sin(t292);
t279 = cos(t292);
t297 = sin(qJ(2,2));
t303 = cos(qJ(2,2));
t260 = -t276 * t303 + t279 * t297;
t283 = 0.1e1 / t297;
t359 = t260 * t283;
t261 = t297 * t276 + t279 * t303;
t358 = t261 * t283;
t293 = legFrame(1,2);
t277 = sin(t293);
t280 = cos(t293);
t299 = sin(qJ(2,1));
t305 = cos(qJ(2,1));
t262 = -t277 * t305 + t280 * t299;
t284 = 0.1e1 / t299;
t357 = t262 * t284;
t263 = t299 * t277 + t280 * t305;
t356 = t263 * t284;
t294 = sin(qJ(3,3));
t264 = t294 * rSges(3,1) + t300 * rSges(3,2);
t355 = t264 * t285;
t296 = sin(qJ(3,2));
t265 = t296 * rSges(3,1) + t302 * rSges(3,2);
t354 = t265 * t287;
t298 = sin(qJ(3,1));
t266 = t298 * rSges(3,1) + t304 * rSges(3,2);
t353 = t266 * t289;
t352 = t282 * t285;
t351 = t282 * t294;
t350 = t283 * t287;
t349 = t283 * t296;
t348 = t284 * t289;
t347 = t284 * t298;
t314 = 0.1e1 / pkin(2);
t346 = t285 * t314;
t345 = 0.1e1 / t300 ^ 2 * t301;
t344 = t287 * t314;
t343 = 0.1e1 / t302 ^ 2 * t303;
t342 = t289 * t314;
t341 = 0.1e1 / t304 ^ 2 * t305;
t340 = t312 + t313;
t272 = -rSges(3,2) * t363 + Icges(3,6);
t273 = rSges(3,1) * t363 - Icges(3,5);
t255 = t272 * t300 - t273 * t294;
t339 = t255 * t282 * t366;
t256 = t272 * t302 - t273 * t296;
t338 = t256 * t283 * t365;
t257 = t272 * t304 - t273 * t298;
t337 = t257 * t284 * t364;
t336 = t275 * t352;
t335 = t275 * t346;
t334 = t276 * t350;
t333 = t276 * t344;
t332 = t277 * t348;
t331 = t277 * t342;
t330 = t278 * t352;
t329 = t278 * t346;
t328 = t279 * t350;
t327 = t279 * t344;
t326 = t280 * t348;
t325 = t280 * t342;
t324 = t285 * t351;
t323 = t287 * t349;
t322 = t289 * t347;
t321 = t314 * t345;
t320 = t314 * t343;
t319 = t314 * t341;
t318 = t345 * t351;
t317 = t343 * t349;
t316 = t341 * t347;
t315 = Icges(2,3) + (rSges(2,1) ^ 2 + rSges(2,2) ^ 2) * m(2) + (0.2e1 * rSges(3,3) ^ 2 + t340) * t368 + t367 + Icges(3,1) / 0.2e1;
t311 = 0.2e1 * qJ(3,1);
t310 = 0.2e1 * qJ(3,2);
t309 = 0.2e1 * qJ(3,3);
t308 = m(2) * rSges(2,1);
t281 = m(1) + m(2) + m(3);
t274 = -m(3) * rSges(3,1) * rSges(3,2) + Icges(3,4);
t271 = m(2) * rSges(2,2) - t363;
t270 = t340 * m(3) + Icges(3,3);
t254 = (t308 + (rSges(3,1) * t304 - rSges(3,2) * t298) * m(3)) * t305 - t299 * t271;
t253 = (t308 + (rSges(3,1) * t302 - rSges(3,2) * t296) * m(3)) * t303 - t297 * t271;
t252 = (t308 + (rSges(3,1) * t300 - rSges(3,2) * t294) * m(3)) * t301 - t295 * t271;
t251 = cos(t311) * t362 + t274 * sin(t311) + t315;
t250 = cos(t310) * t362 + t274 * sin(t310) + t315;
t249 = cos(t309) * t362 + t274 * sin(t309) + t315;
t248 = (-t254 * t325 + t263 * t281) * t284;
t247 = (t254 * t331 + t262 * t281) * t284;
t246 = (-t253 * t327 + t261 * t281) * t283;
t245 = (t253 * t333 + t260 * t281) * t283;
t244 = (-t252 * t329 + t259 * t281) * t282;
t243 = (t252 * t335 + t258 * t281) * t282;
t242 = -t299 * m(3) * t266 * t342 + (-t254 * t319 + t281 * t289) * t347;
t241 = -t297 * m(3) * t265 * t344 + (-t253 * t320 + t281 * t287) * t349;
t240 = -t295 * m(3) * t264 * t346 + (-t252 * t321 + t281 * t285) * t351;
t239 = (-t251 * t325 + t254 * t263) * t284;
t238 = (t251 * t331 + t254 * t262) * t284;
t237 = (-t250 * t327 + t253 * t261) * t283;
t236 = (t250 * t333 + t253 * t260) * t283;
t235 = (-t249 * t329 + t252 * t259) * t282;
t234 = (t249 * t335 + t252 * t258) * t282;
t233 = t257 * t342 + (-t251 * t319 + t254 * t289) * t347;
t232 = t256 * t344 + (-t250 * t320 + t253 * t287) * t349;
t231 = t255 * t346 + (-t249 * t321 + t252 * t285) * t351;
t1 = [t244 * t360 + t246 * t358 + t248 * t356 + m(4) + (-t235 * t330 - t237 * t328 - t239 * t326) * t314, t244 * t361 + t246 * t359 + t248 * t357 + (t235 * t336 + t237 * t334 + t239 * t332) * t314, t244 * t324 + t246 * t323 + t248 * t322 + (-t235 * t318 - t237 * t317 - t239 * t316 + (-t278 * t339 - t279 * t338 - t280 * t337) * t314 + (-t259 * t355 - t261 * t354 - t263 * t353) * m(3)) * t314; t243 * t360 + t245 * t358 + t247 * t356 + (-t234 * t330 - t236 * t328 - t238 * t326) * t314, t243 * t361 + t245 * t359 + t247 * t357 + m(4) + (t234 * t336 + t236 * t334 + t238 * t332) * t314, t243 * t324 + t245 * t323 + t247 * t322 + (-t234 * t318 - t236 * t317 - t238 * t316 + (t275 * t339 + t276 * t338 + t277 * t337) * t314 + (-t258 * t355 - t260 * t354 - t262 * t353) * m(3)) * t314; t240 * t360 + t241 * t358 + t242 * t356 + (-t231 * t330 - t232 * t328 - t233 * t326) * t314, t240 * t361 + t241 * t359 + t242 * t357 + (t231 * t336 + t232 * t334 + t233 * t332) * t314, t240 * t324 + t241 * t323 + t242 * t322 + m(4) + (-t231 * t318 - t232 * t317 - t233 * t316 + ((-t257 * t316 + t270 * t289) * t289 + (-t256 * t317 + t270 * t287) * t287 + (-t255 * t318 + t270 * t285) * t285) * t314 + (-t264 * t366 * t294 - t265 * t365 * t296 - t266 * t364 * t298) * m(3)) * t314;];
MX  = t1;
